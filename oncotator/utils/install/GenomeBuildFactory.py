import logging
from shove.core import Shove
from oncotator.Transcript import Transcript
from oncotator.index.gaf import region2bin
from oncotator.utils.install.GenomeBuildInstallUtils import GenomeBuildInstallUtils
from BCBio import GFF
from Bio import SeqIO

class GenomeBuildFactory(object):
    """ Responsible for creating indices for genome builds (through ENSEMBL) and creating a set of datasource files.
        The methods in this class would typically be run in datasource creation, not during annotation.
    """
    def __init__(self):
        self._transcript_index = dict()

    def _convertGFFRecordToTranscript(self, gff_record, seq_dict):
        quals = gff_record['quals']
        transcript_id = quals['transcript_id'][0]

        if transcript_id not in self._transcript_index.keys():
            self._transcript_index[transcript_id] = Transcript(transcript_id, gene=quals['gene_name'][0], gene_id=quals['gene_id'][0], contig=gff_record['rec_id'])

        if gff_record['type'] == 'exon':
            self._transcript_index[transcript_id].add_exon(gff_record['location'][0], gff_record['location'][1])
        elif gff_record['type'] == 'CDS':
            self._transcript_index[transcript_id].add_cds(gff_record['location'][0], gff_record['location'][1])

        seq = seq_dict.get(transcript_id, None)
        if seq is not None:
            genome_seq_as_str = str(seq.seq)
        else:
            genome_seq_as_str = ""

        self._transcript_index[transcript_id].set_seq(genome_seq_as_str)
        return self._transcript_index[transcript_id]

    def build_ensembl_transcript_index(self, ensembl_input_gtf, ensembl_input_fasta, output_filename, protocol="file"):
        """Create the transcript index (using shove) for ensembl.  Key is transcript ID
        :param ensembl_input_gtf:
        :param ensembl_input_fasta: sequence data for transcripts corresponding to what is in the gtf
        :param output_filename:
        :param protocol: shove protocol.  Usually "file" or "sqlite"
        """

        # Example code taken from http://biopython.org/wiki/GFF_Parsing
        shove = Shove(protocol + "://" + output_filename, "memory://")

        in_seq_file = ensembl_input_fasta
        in_seq_handle = open(in_seq_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
        in_seq_handle.close()

        in_file = ensembl_input_gtf
        in_handle = open(in_file)
        for rec in GFF.parse_simple(in_file): #(in_handle, base_dict=seq_dict):

            # transcript id seems to always be a list of length 1
            if len(rec['quals']['transcript_id']) > 1:
                logging.getLogger(__name__).warn("ensembl records had more than one transcript id: " + str(rec['quals']['transcript_id']))
            transcript_id = rec['quals']['transcript_id'][0]
            shove[transcript_id] = self._convertGFFRecordToTranscript(rec, seq_dict)
        shove.close()
        in_handle.close()

    def build_ensembl_transcripts_by_gene_index(self, ensembl_transcript_index_fname, output_filename, protocol="file"):
        """ Create an index for gene --> transcripts using a transcript index created in build_ensembl_transcript_index
        :param ensembl_transcript_index_fname: file/dir location for ensembl transcript db
        :return:
        """

        #TODO: This may need to be moved to the init of the transcript datasource as that may be faster.

        transcript_db = Shove(protocol + "://" + ensembl_transcript_index_fname, "memory://")
        output_db = Shove(protocol + "://" + output_filename, "memory://")

        transcript_keys = transcript_db.keys()

        for tx_id in transcript_keys:
            tx = transcript_db[tx_id]
            if tx.get_gene() not in output_db:
                output_db[tx.get_gene()] = [tx]
            else:
                output_db[tx.get_gene()].append(tx)

        output_db.close()
        transcript_db.close()

    def build_ensembl_transcripts_by_genomic_location_index(self, ensembl_transcript_index_fname, output_filename, protocol="file"):
        """Create an index for genomic position to transcripts index, using a transcript index created in
            build_ensembl_transcript_index
        """
        transcript_db = Shove(protocol + "://" + ensembl_transcript_index_fname, "memory://")
        output_db = Shove(protocol + "://" + output_filename, "memory://")

        transcript_keys = transcript_db.keys()

        for tx_id in transcript_keys:
            tx = transcript_db[tx_id]
            start = tx.get_start()
            end = tx.get_end()
            genomic_location_bin = region2bin(start, end)
            key = tx.get_contig() + "_" + str(genomic_location_bin)
            if key not in output_db:
                output_db[key] = [tx]
            else:
                output_db[key].append([tx])
        output_db.close()
        transcript_db.close()

    def build_ensembl_sequences_by_transcript_id_index(self, ensembl_transcript_index_fname):
        """ Create an index that returns the sequences from the fasta file given the transcript id as a key.
        :return:
        """
        #TODO: This method may need different parameters.
        raise NotImplementedError("This method is not finished.")


