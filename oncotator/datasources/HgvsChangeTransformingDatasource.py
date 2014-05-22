import re
import logging

from Bio import Seq

from oncotator.datasources.ChangeTransformingDatasource import ChangeTransformingDatasource
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.utils.VariantClassifier import VariantClassifier
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils

AA_NAME_MAPPING = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val',
    '*': '*',
}

PROT_REGEXP = re.compile('p\.([A-Z*])(\d+)([A-Z*])')
PROT_FRAMESHIFT_REGEXP = re.compile('p\.([A-Z*])(\d+)fs')
PROT_INFRAME_DEL_REGEXP = re.compile('p\.([A-Z*])(\d+)del')
PROT_INFRAME_INS_REGEXP = re.compile('p\.(\d+)_(\d+)[A-Z*]>([A-Z*])')

class HgvsChangeTransformingDatasource(ChangeTransformingDatasource):

    def __init__(self, src_file, title='', version=None):
        super(ChangeTransformingDatasource, self).__init__(src_file, title=title, version=version)
        self.output_headers = ['HGVS_genomic_change', 'HGVS_coding_DNA_change', 'HGVS_protein_change']

        self.vcer = VariantClassifier()

        #### HACK #####
        self.gencode_ds = EnsemblTranscriptDatasource(src_file)

        #### HACK #####
        ensembl_id_mapping_fname = '/Users/aramos/Desktop/ot_hgvs/ensembl_id_mappingsGRCh37.p13.txt'
        import pandas as pd
        self.ensembl_id_mapping = dict()
        with open(ensembl_id_mapping_fname) as in_fh:
            for line in in_fh:
                tx_id, prot_id = line.split('\t')[1:]
                self.ensembl_id_mapping[tx_id] = prot_id.replace('\n', '')
        
    def annotate_mutation(self, mutation):
        # TODO
        # TODO: When annotating, you should probably give the annotation a name like:  self.title + "_" + output_annotation_name
        # You can assume that the annotations from the transript datasource have already been populated (such as protein_change and genome_change)

#        if mutation['start'] in ['67969313', '67970596']:
#            from IPython import embed; embed()

        hgvs_genomic_change = self._adjust_genome_change(mutation)
        hgvs_coding_dna_change = self._adjust_coding_DNA_change(mutation)
        hgvs_protein_change = self._adjust_protein_change(mutation)

        mutation.createAnnotation('HGVS_genomic_change', hgvs_genomic_change, self.title)
        mutation.createAnnotation('HGVS_coding_DNA_change', hgvs_coding_dna_change, self.title)
        mutation.createAnnotation('HGVS_protein_change', hgvs_protein_change, self.title)

        return mutation

    def _adjust_genome_change(self, mutation):
        if mutation['variant_type'] == 'DNP':
            #from IPython import embed; embed()
            change = '%d_%ddelins%s' % (mutation['start'], mutation['end'], mutation['alt_allele'])
        elif mutation['variant_type'] == 'INS':
            genome_pos_adjust = self._genome_pos_adjust_if_duplication(mutation['alt_allele'], mutation['ref_context'])
            if genome_pos_adjust > 0:
                start_pos = mutation['end'] + genome_pos_adjust
                end_pos = start_pos + len(mutation['alt_allele']) - 1
                change = '%d_%ddup%s' % (start_pos, end_pos, mutation['alt_allele'])
            else:
                change = '%d_%dins%s' % (mutation['start'], mutation['end'], mutation['alt_allele'])
        elif mutation['variant_type'] == 'DEL':
            if mutation['start'] == mutation['end']:
                change = str(mutation['start']) if mutation['start'] == mutation['end'] else '%d_%d' % (mutation['start'], mutation['end'])
                change = '%sdel%s' % (change, mutation['ref_allele'])
            else:
                change = '%d_%ddel%s' % (mutation['start'], mutation['end'], mutation['ref_allele'])
        else:
            change = '%d%s>%s' % (mutation['start'], mutation['ref_allele'], mutation['alt_allele'])

        chrn = mutation['chr']
        if mutation['build'] == '':
            adjusted_genome_change = '%s:g.%s' % (chrn, change)
        else:
            if mutation['build'] == 'hg19' and not chrn.startswith('chr'):
                chrn = 'chr' + chrn
            adjusted_genome_change = '%s.%s:g.%s' % (chrn, mutation['build'], change)

        return adjusted_genome_change

    def _adjust_coding_DNA_change(self, mutation):
        if mutation['variant_type'] == 'DNP':
            return self._get_cdna_change_for_ONP(mutation)
        else:
            vc = mutation['variant_classification']
            if vc in ['Intron', 'Splice_Site']:
                return self._get_cdna_change_for_intron(mutation)
            elif vc in ["5'UTR", 'De_novo_Start_OutOfFrame']:
                return self._get_cdna_change_for_5_utr(mutation)
            elif vc == "3'UTR":
                return self._get_cdna_change_for_3_utr(mutation)
            elif vc == 'IGR':
                return ''
            elif vc.endswith('Ins'):
                return self._get_cdna_change_for_ins(mutation)
            else:
                return '%s:%s' % (mutation['annotation_transcript'], mutation['transcript_change'])

    def _adjust_protein_change(self, mutation):
        if mutation['variant_classification'] in ['Intron', 'IGR', "3'UTR", "5'UTR", 'RNA',
            'lincRNA', 'Silent', 'Splice_Site', 'De_novo_Start_OutOfFrame']:
            return ''
        elif mutation['variant_classification'] in ['Frame_Shift_Ins', 'Frame_Shift_Del']:
            regx_res = PROT_FRAMESHIFT_REGEXP.match(mutation['protein_change'])
            ref_aa, aa_pos = [regx_res.group(i) for i in range(1, 3)]
            adjusted_prot_change = 'p.%s%sfs' % (AA_NAME_MAPPING[ref_aa], aa_pos)
        elif mutation['variant_classification'] == 'In_Frame_Del':
            regx_res = PROT_INFRAME_DEL_REGEXP.match(mutation['protein_change'])
            ref_aa, aa_pos = [regx_res.group(i) for i in range(1, 3)]
            adjusted_prot_change = 'p.%s%sdel' % (AA_NAME_MAPPING[ref_aa], aa_pos)
        elif mutation['variant_classification'] == 'In_Frame_Ins':
            regx_res = PROT_INFRAME_INS_REGEXP.match(mutation['protein_change'])
            aa_pos_1, aa_pos_2, alt_allele = [regx_res.group(i) for i in range(1, 4)]
            adjusted_prot_change = 'p.%s_%sins%s' % (aa_pos_1, aa_pos_2, alt_allele)
        else:
            regx_res = PROT_REGEXP.match(mutation['protein_change'])
            ref_aa, aa_pos, alt_aa = [regx_res.group(i) for i in range(1, 4)]
            adjusted_prot_change = 'p.%s%s%s' % (AA_NAME_MAPPING[ref_aa], aa_pos, AA_NAME_MAPPING[alt_aa])

        prot_id = self._get_ensembl_prot_id_from_tx_id(mutation['annotation_transcript'])
        adjusted_prot_change = '%s:%s' % (prot_id, adjusted_prot_change)
        return adjusted_prot_change

    def _get_ensembl_prot_id_from_tx_id(self, tx_id):
        if '.' in tx_id:
            tx_id = tx_id[:tx_id.index('.')]
        return self.ensembl_id_mapping[tx_id]

    def _get_cdna_change_for_intron(self, mutation):
        #### HACK #####
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        nearest_exon = self.vcer._determine_closest_exon(tx, int(mutation['start']),
            int(mutation['end']))
        dist_to_exon = self.vcer._get_splice_site_coordinates(tx, int(mutation['start']),
            int(mutation['end']), nearest_exon)
        tx_exons = tx.get_exons()
    
        if dist_to_exon < 0 and mutation['transcript_strand'] == '+':
            cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(
                tx_exons[nearest_exon][0], tx_exons[nearest_exon][0], tx)
            if mutation['variant_classification'] == 'Intron':
                #do not do for splice sites
                dist_to_exon = dist_to_exon - 1 #why substract 1? why only for positive strand transcript
        elif dist_to_exon < 0 and mutation['transcript_strand'] == '-':
            cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(
                tx_exons[nearest_exon][1], tx_exons[nearest_exon][1], tx)
        elif dist_to_exon > 0 and mutation['transcript_strand'] == '+':
            cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(
                tx_exons[nearest_exon][1], tx_exons[nearest_exon][1], tx)
        elif dist_to_exon > 0 and mutation['transcript_strand'] == '-':
            cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(
                tx_exons[nearest_exon][0], tx_exons[nearest_exon][0], tx)

        if cds_position_of_nearest_exon[0] == 0:
            #this means intron occurs before start codon and thus cds position needs to be adjusted to a negative number
            exon_coords = TranscriptProviderUtils.convert_genomic_space_to_exon_space(mutation['start'],
                tx.determine_cds_start() + 1, tx)
            cds_position_of_nearest_exon = exon_coords[0] - exon_coords[1]
        else:
            cds_position_of_nearest_exon = cds_position_of_nearest_exon[0]

        if dist_to_exon < 0:   
            cds_position_of_nearest_exon += 1 #why add 1?

#        if mutation['start'] == 80529551:
#            from IPython import embed; embed()
        
        sign = '-' if dist_to_exon < 0 else '+'
        dist_to_exon = abs(dist_to_exon)

        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)

        adjusted_tx_change = 'c.%d%s%d%s>%s' % (cds_position_of_nearest_exon, sign, dist_to_exon,
            tx_ref_allele, tx_alt_allele)

        return '%s:%s' % (mutation['annotation_transcript'], adjusted_tx_change)


#TranscriptProviderUtils.convert_genomic_space_to_cds_space(80529551, 80529551, tx)
#TranscriptProviderUtils.convert_genomic_space_to_cds_space(tx_exons[nearest_exon][0], tx_exons[nearest_exon][0], tx)
#self.vcer._get_splice_site_coordinates(tx, 484634, 484634, nearest_exon)
#self.vcer._get_splice_site_coordinates(tx, 484814, 484814, nearest_exon)

    def _get_cdna_change_for_5_utr(self, mutation):
        #### HACK #####
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        dist_from_cds_start = abs(mutation['start'] - tx.determine_cds_start())
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        adjusted_tx_change = 'c.-%d%s>%s' % (dist_from_cds_start, tx_ref_allele, tx_alt_allele)
        return '%s:%s' % (mutation['annotation_transcript'], adjusted_tx_change)

    def _get_cdna_change_for_3_utr(self, mutation):
        #### HACK #####
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        dist_from_cds_stop = abs(mutation['start'] - tx.determine_cds_stop() + 2)
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        adjusted_tx_change = 'c.*%d%s>%s' % (dist_from_cds_stop, tx_ref_allele, tx_alt_allele)
        return '%s:%s' % (mutation['annotation_transcript'], adjusted_tx_change)

    def _get_tx_alleles(self, mutation):
        tx_ref_allele, tx_alt_allele = mutation['ref_allele'], mutation['alt_allele']
        if mutation['transcript_strand'] == '-':
            tx_ref_allele, tx_alt_allele = Seq.reverse_complement(tx_ref_allele), Seq.reverse_complement(tx_alt_allele)
        return tx_ref_allele, tx_alt_allele

    def _get_cdna_change_for_ONP(self, mutation):
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        tx_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(mutation['start'], mutation['end'], tx)
        tx_pos = tuple(t + 1 for t in tx_pos)
        adjusted_tx_change = 'c.%d_%ddelins%s' % (tx_pos[0], tx_pos[1], mutation['alt_allele'])
        return '%s:%s' % (mutation['annotation_transcript'], adjusted_tx_change)

    def _genome_pos_adjust_if_duplication(self, alt_allele, ref_context):
        half_len = len(ref_context) / 2
        downstream_seq = ref_context[half_len:]
        return self._adjust_pos_if_repeated_in_seq(alt_allele, downstream_seq)

    def _get_cdna_change_for_ins(self, mutation):
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        tx_pos = TranscriptProviderUtils.convert_genomic_space_to_exon_space(mutation['start'], mutation['end'], tx)
        downstream_seq = tx.get_seq()[tx_pos[0]:]
        coding_pos_adjust = self._adjust_pos_if_repeated_in_seq(tx_alt_allele, downstream_seq)
        coding_pos_adjust = coding_pos_adjust - len(mutation['alt_allele'])

        coding_tx_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(mutation['start'], mutation['end'], tx)
        if coding_pos_adjust > 0:
            coding_tx_pos = tuple(t + 1 for t in coding_tx_pos)
            start_pos = coding_tx_pos[0] + coding_pos_adjust
            end_pos = start_pos + len(mutation['alt_allele']) - 1
            change = 'c.%d_%ddup%s' % (start_pos, end_pos, mutation['alt_allele'])
        else:
            change = 'c.%d_%dins%s' % (coding_tx_pos[0], coding_tx_pos[1], mutation['alt_allele'])
            #from IPython import embed; embed()

        return '%s:%s' % (mutation['annotation_transcript'], change)

    def _adjust_pos_if_repeated_in_seq(self, allele, downstream_seq):
        allele_len = len(allele)
        adjust_amt = 0
        for i in range(0, len(downstream_seq), allele_len):
            if allele == downstream_seq[i:i + allele_len]:
                adjust_amt += allele_len
            else:
                break

        #from IPython import embed; embed()
        return adjust_amt


