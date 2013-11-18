from copy import deepcopy
import logging
from shove.core import Shove
from oncotator.Transcript import Transcript
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.index.gaf import region2bin
from oncotator.utils.MutUtils import MutUtils
from oncotator.utils.install.GenomeBuildInstallUtils import GenomeBuildInstallUtils
from BCBio import GFF
from Bio import SeqIO
from Bio import Seq

class GenomeBuildFactory(object):
    """ Responsible for creating indices for genome builds (through ENSEMBL) and creating a set of datasource files.
        The methods in this class would typically be run in datasource creation, not during annotation.

        NOTE: This class only supports ENSEMBL/GENCODE, but there is only a minimal amount of ENSEMBL/GENCODE specific
        code
    """

    QUALS_TO_CHECK = ['gene_status', 'level', 'source', 'tag', 'ccdsid', 'transcript_name', 'transcript_type', 'havana_gene', 'havana_transcript', 'transcript_status']

    def __init__(self):
        self._transcript_index = dict()

    def _determine_protein_seq(self, tx):
        cds_start, cds_stop = tx.determine_cds_footprint()
        if cds_start == -1 or cds_stop == -1:
            return ""
        protein_seq = self._calculate_protein_sequence(tx.get_exons(), tx.get_seq(), cds_start, cds_stop, tx.get_strand())
        protein_seq = ''.join([protein_seq, '*'])
        return protein_seq

    def _calculate_protein_sequence(self, exons, seq, cds_start_genomic_space, cds_stop_genomic_space, strand):
        cds_start_exon_space, cds_stop_exon_space = TranscriptProviderUtils._convert_genomic_space_to_feature_space(int(cds_start_genomic_space), int(cds_stop_genomic_space), exons, strand)

        prot_seq = Seq.translate(seq[int(cds_start_exon_space):int(cds_stop_exon_space)])
        if len(prot_seq) > 0 and prot_seq[-1] == '*':
            prot_seq = prot_seq[:-1]

        return prot_seq

    def _convertGFFRecordToTranscript(self, gff_record, seq_dict, seq_dict_keys):
        """

        :param gff_record:
        :param seq_dict:
        :return: None if the record is a gene record or otherwise does not represent a transcript, CDS, *_codon, or exon
        """
        types_of_interest = ["exon", "CDS", "start_codon", "stop_codon"]
        if gff_record['type'] not in types_of_interest:
            return None

        quals = gff_record['quals']
        transcript_id = quals['transcript_id'][0]


        if transcript_id not in self._transcript_index.keys():
            contig = MutUtils.convertChromosomeStringToMutationDataFormat(gff_record['rec_id'])
            tx = Transcript(transcript_id, gene=quals['gene_name'][0], gene_id=quals['gene_id'][0], contig=contig)
            self._transcript_index[transcript_id] = tx

            # Set the gene_type based on gene_type or gene_biotype
            key = "gene_biotype"
            if key not in quals.keys():
                key = "gene_type"
            self._transcript_index[transcript_id].set_gene_type(quals.get(key, [""])[0])

            if gff_record['strand'] == 1:
                self._transcript_index[transcript_id].set_strand("+")
            else:
                self._transcript_index[transcript_id].set_strand("-")
            qual_keys = quals.keys()
            for attribute in GenomeBuildFactory.QUALS_TO_CHECK:
                if attribute in qual_keys:
                    self._transcript_index[transcript_id].add_other_attribute(attribute, quals[attribute])

            seq = seq_dict.get(transcript_id, None)

            if seq is not None:
                genome_seq_as_str = str(seq.seq)
            else:
                genome_seq_as_str = ""

            self._transcript_index[transcript_id].set_seq(genome_seq_as_str)

        gff_type = gff_record['type']
        if gff_type == 'exon':
            self._transcript_index[transcript_id].add_exon(gff_record['location'][0], gff_record['location'][1], quals['exon_number'][0])
        elif gff_type == 'CDS':
            self._transcript_index[transcript_id].add_cds(gff_record['location'][0], gff_record['location'][1])
        elif gff_type == 'start_codon':
            self._transcript_index[transcript_id].set_start_codon(gff_record['location'][0], gff_record['location'][1])
        elif gff_type == 'stop_codon':
            self._transcript_index[transcript_id].set_stop_codon(gff_record['location'][0], gff_record['location'][1])

    def _create_seq_dict(self, seq_fasta_fp):
        """Create a dictionary with keys to the sequenced bases.  Note that this includes strand.
        So if the ref is CGAT and strand is negative, this dict will have ATCG for the sequence

        This method assumes that no transcripts have "|" or ">" in the ID.

        This method will extract the transcript id from the GENCODE fasta.  Example:
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|

        :param seq_fasta_fp: open readable filepointer to a fasta file (.fa)
        :returns seq_dict: keys (usually transcript id, but can be more)
        """
        seq_dict = SeqIO.to_dict(SeqIO.parse(seq_fasta_fp, "fasta"))

        # Check for faulty parsing and create new entries that are simply the transcript id
        # For example: >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        seq_dict_keys = seq_dict.keys()
        for k in seq_dict_keys:
            key_list = k.split("|")
            if len(key_list) != 1:
                for id in key_list:
                    new_id = id.replace(">", "")
                    if new_id.startswith("ENST"):
                        seq_dict[new_id] = seq_dict[k]

        return seq_dict

    def build_ensembl_transcript_index(self, ensembl_input_gtfs, ensembl_input_fastas, output_filename, protocol="file"):
        """Create the transcript index (using shove) for ensembl.  Key is transcript ID.

        Note:  This method will hold the entire transcript index in RAM.

        :param ensembl_input_gtfs: (list)
        :param ensembl_input_fastas: (list) sequence data for transcripts corresponding to what is in the gtfs
        :param output_filename:
        :param protocol: shove protocol.  Usually "file" or "sqlite"
        """

        # Example code taken from http://biopython.org/wiki/GFF_Parsing
        shove = Shove(protocol + "://" + output_filename, "memory://")
        logging.getLogger(__name__).info("Transcript index being created: " + protocol + "://" + output_filename)

        seq_dict = {}
        for in_seq_file in ensembl_input_fastas:
            in_seq_handle = open(in_seq_file)
            seq_dict.update(self._create_seq_dict(in_seq_handle))
            in_seq_handle.close()
            logging.getLogger(__name__).info("Parsed fasta file: " + in_seq_file)

        for file_ctr, in_file in enumerate(ensembl_input_gtfs):
            in_handle = open(in_file)
            seq_dict_keys = seq_dict.keys()
            ctr = 0
            for rec in GFF.parse_simple(in_file): #(in_handle, base_dict=seq_dict):

                # transcript id seems to always be a list of length 1
                if len(rec['quals']['transcript_id']) > 1:
                    logging.getLogger(__name__).warn("ensembl records had more than one transcript id: " + str(rec['quals']['transcript_id']))

                self._convertGFFRecordToTranscript(rec, seq_dict, seq_dict_keys)
                ctr += 1
                if (ctr % 10000) == 0:
                    logging.getLogger(__name__).info("Added " + str(ctr) + " lines of gtf " + str(file_ctr+1) + " of " + str(len(ensembl_input_gtfs)) + " (" + in_file + ") into internal transcript index.")
            in_handle.close()
            logging.getLogger(__name__).info("Finished " + str(ctr) + " lines of gtf (" + in_file + ")")

        logging.getLogger(__name__).info("Populating final db with internal transcript index.")
        transcript_index_keys = self._transcript_index.keys()
        for k in transcript_index_keys:
            protein_sequence = self._determine_protein_seq(self._transcript_index[k])
            self._transcript_index[k].set_protein_seq(protein_sequence)
            shove[k] = self._transcript_index[k]

        logging.getLogger(__name__).info("Transcript index created " + str(len(shove.keys())) + " transcripts.")
        shove.close()


    def build_ensembl_transcripts_by_gene_index(self, ensembl_transcript_index_fname, output_filename, protocol="file"):
        """ Create an index for gene --> transcripts using a transcript index created in build_ensembl_transcript_index
        :param ensembl_transcript_index_fname: file/dir location for ensembl transcript db
        :return:
        """

        #TODO: This may need to be moved to the init of the transcript datasource as that may be faster.

        transcript_db = Shove(protocol + "://" + ensembl_transcript_index_fname, "memory://")
        output_db = Shove(protocol + "://" + output_filename, "memory://", optimize=False)

        transcript_keys = transcript_db.keys()

        for tx_id in transcript_keys:
            tx = transcript_db[tx_id]
            gene = tx.get_gene()
            try:
                tmpList = output_db[gene]
            except KeyError:
                output_db[gene] = []
                tmpList = output_db[gene]
            tmpList.append(tx)
            output_db[gene] = tmpList

        output_db.close()
        transcript_db.close()

    def build_ensembl_transcripts_by_genomic_location_index(self, ensembl_transcript_index_fname, output_filename, protocol="file"):
        """Create an index for genomic position to transcripts index, using a transcript index created in
            build_ensembl_transcript_index
        """
        transcript_db = Shove(protocol + "://" + ensembl_transcript_index_fname)
        output_db = Shove(protocol + "://" + output_filename, optimize=False)

        transcript_keys = transcript_db.keys()

        for tx_id in transcript_keys:
            tx = transcript_db[tx_id]
            start = tx.get_start()
            end = tx.get_end()
            genomic_location_bin = region2bin(start, end)
            key = tx.get_contig() + "_" + str(genomic_location_bin)
            try:
                tmpList = output_db[key]
            except KeyError:
                output_db[key] = []
                tmpList = output_db[key]

            tmpList.append(tx)
            output_db[key] = tmpList

        output_db.close()
        transcript_db.close()

    def construct_ensembl_indices(self, ensembl_input_gtfs, ensembl_input_fastas, base_output_filename):
        """

        :param ensembl_input_gtfs: (list) gtf input files
        :param ensembl_input_fastas: (list) fasta input files
        :param base_output_filename: Just the base output filename, such as "my_ensembl" without any extensions.
        :return:
        """
        ensembl_transcript_index_filename = base_output_filename + ".transcript.idx"
        self.build_ensembl_transcript_index(ensembl_input_gtfs, ensembl_input_fastas, ensembl_transcript_index_filename)
        self.build_ensembl_transcripts_by_gene_index(ensembl_transcript_index_filename, base_output_filename + ".transcript_by_gene.idx")
        self.build_ensembl_transcripts_by_genomic_location_index(ensembl_transcript_index_filename, base_output_filename + ".transcript_by_gp_bin.idx")
