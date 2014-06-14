# LICENSE_GOES_HERE

import collections
import logging
import pysam

from oncotator.datasources.Datasource import Datasource


class Cosmic(Datasource):
    """
    Cosmic Datasource expecting a tabix compressed tsv file from Cosmic.

    Retrieves annotations based on specific genomic position and entire transcript protein sequence.

    Preferred input annotations:
    transcript_protein_position_start -- used to get COSMIC annotations by transcript protein seq
    transcript_protein_position_end -- used to get COSMIC annotations by transcript protein seq

    Output columns are:
        COSMIC_n_overlapping_mutations
            Number of COSMIC entries that overlap by genomic position.

        COSMIC_overlapping_mutation_AAs
            Protein change of overlapping COSMIC entries.  Number of entries with protein change is in parentheses.  Multiple protein
            changes are delimited with "|" character.

        COSMIC_overlapping_mutation_descriptions
            COSMIC mutation description of overlapping COSMIC entries.  Number of entries with mutation description is in parentheses.  Multiple mutation
            descriptions are delimited with "|" character.

        COSMIC_overlapping_primary_sites
            Primary site of overlapping COSMIC entries.  Number of entries with primary site is in parentheses.  Multiple primary
            sites are delimited with "|" character.


        NOTE: fusion genes are handled in a different datasource.

    """
    def __init__(self, src_file, title='COSMIC', version=None, gpp_tabix_file=None):
        self.title = title

        super(Cosmic, self).__init__(src_file, title=self.title, version=version)

        if gpp_tabix_file is None:
            raise ValueError("A second index by gene protein position must be specified.")

        self.db_genomePos = pysam.Tabixfile(src_file)
        header = self.db_genomePos.header.next()
        self.src_headers = header.lstrip('#').strip().split('\t')

        self.db_geneProteinPos = pysam.Tabixfile(gpp_tabix_file)
        gppHeader = self.db_geneProteinPos.header.next()
        self.gpp_headers = gppHeader.lstrip('#').strip().split('\t')

        self.output_headers = ['COSMIC_n_overlapping_mutations', 'COSMIC_overlapping_mutation_AAs',
            'COSMIC_overlapping_mutation_descriptions', 'COSMIC_overlapping_primary_sites']
        self.logger = logging.getLogger(__name__)

    def annotate_mutation(self, mutation):
        # transcriptPos = [mutation['tx_start'],mutation['tx_end']]
        overlapping_cosmic_entries = self.fetch_overlapping_cosmic_records(mutation.chr, mutation.start, mutation.end)

        overlappingGeneProteinPositionEntries = []
        if "transcript_protein_position_start" in mutation.keys() and "transcript_protein_position_end" in mutation.keys():
            overlappingGeneProteinPositionEntries = self.fetch_overlapping_gene_proteinPos_records(mutation['gene'],mutation["transcript_protein_position_start"], mutation["transcript_protein_position_end"])

        overlapping_cosmic_entries.extend(overlappingGeneProteinPositionEntries)

        n_overlapping_mutations = len(overlapping_cosmic_entries)

        mutation_AAs = collections.Counter([entry['Mutation AA'] for entry in overlapping_cosmic_entries])
        mutation_descriptions = collections.Counter([entry['Mutation Description'] for entry in overlapping_cosmic_entries])
        primary_sites = collections.Counter([entry['Primary site'] for entry in overlapping_cosmic_entries])

        get_output_str = lambda counter_obj: '|'.join(['%s(%d)' % s for s in sorted(counter_obj.items(), key=lambda x: x[1], reverse=True)])

        overlapping_mutations_AA_str = get_output_str(mutation_AAs)
        overlapping_mutation_descriptions_str = get_output_str(mutation_descriptions)
        overlapping_primary_sites_str = get_output_str(primary_sites)

        mutation.createAnnotation('COSMIC_n_overlapping_mutations', str(n_overlapping_mutations), self.title)
        mutation.createAnnotation('COSMIC_overlapping_mutation_AAs', overlapping_mutations_AA_str, self.title)
        mutation.createAnnotation('COSMIC_overlapping_mutation_descriptions', overlapping_mutation_descriptions_str, self.title)
        mutation.createAnnotation('COSMIC_overlapping_primary_sites', overlapping_primary_sites_str, self.title)

        return mutation

    def _fetchRecords(self, db, headers, chromosome, start, end):
        finalResult = []
        results = None
        try:
            results = self._fetch_from_tabix_file(db, chromosome, start, end)
        except ValueError as ve:
            self.logger.debug(
                "Exception when looking for COSMIC records.  Empty set of records being returned: " + repr(ve))
        if results is not None:
            for res in results:
                finalResult.append(dict(zip(headers, res.strip().split('\t'))))
        return finalResult

    def fetch_overlapping_cosmic_records(self, chromosome, start, end):
        headers = self.src_headers
        db = self.db_genomePos
        return self._fetchRecords(db, headers, chromosome, start, end)

    def fetch_overlapping_gene_proteinPos_records(self, gene, startAA, endAA):
        headers = self.gpp_headers
        db = self.db_geneProteinPos


        if (startAA.strip() == "") or (endAA.strip() == ""):
            return []
        startAA = str(int(startAA) - 1)
        return self._fetchRecords(db, headers, gene, startAA, endAA)

    def _fetch_from_tabix_file(self, tabix_file, chromosome, start, end):
        #TODO: Low priority: This only supports human, if that matters.
        if chromosome == 'X':
            chromosome = "23"
        if chromosome == 'Y':
            chromosome = "24"

        return tabix_file.fetch(region='%s:%s-%s' % (chromosome, str(start), str(end)))