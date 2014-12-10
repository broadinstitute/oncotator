"""
By downloading the PROGRAM you agree to the following terms of use:

BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY

This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").

WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and

WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.

NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:

1. DEFINITIONS
1.1 "PROGRAM" shall mean the object code and source code known as Oncotator 1.0 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.  BROAD acknowledges that the PROGRAM employs one or more public domain code(s) that are freely available for public use.

2. LICENSE
2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.  LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.

3. OWNERSHIP OF INTELLECTUAL PROPERTY
LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.

Copyright 2014 Broad Institute, Inc.
Notice of attribution:  The Oncotator 1.0 program was made available through the generosity of the Broad Institute, Inc.

LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.

4. INDEMNIFICATION
LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorney fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.

5. NO REPRESENTATIONS OR WARRANTIES
THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.

6. ASSIGNMENT
This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.

7. MISCELLANEOUS
7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
"""

from copy import deepcopy
import logging
from shove.core import Shove
from oncotator.Transcript import Transcript
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.index.gaf import region2bin
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.utils.MutUtils import MutUtils
from oncotator.utils.install.GenomeBuildInstallUtils import GenomeBuildInstallUtils
from BCBio import GFF
from Bio import SeqIO

class GenomeBuildFactory(object):
    """ Responsible for creating indices for genome builds (through ENSEMBL) and creating a set of datasource files.
        The methods in this class would typically be run in datasource creation, not during annotation.

        NOTE: This class only supports ENSEMBL/GENCODE, but there is only a minimal amount of ENSEMBL/GENCODE specific
        code
    """

    QUALS_TO_CHECK = ['gene_status', 'level', 'source', 'tag', 'ccdsid', 'transcript_name', 'transcript_type', 'havana_gene', 'havana_transcript', 'transcript_status']

    def __init__(self):
        self._transcript_index = dict()

    def _create_tx_id_to_protein_id_mapping(self, mapping_file):
        """
        Mapping file is assumed to have three columns and be a tsv:
        Ensembl Gene ID, Ensembl Transcript ID, and Ensembl Protein ID

        """
        result = dict()
        if mapping_file is None or mapping_file.strip() == "":
            return result
        tsv_reader = GenericTsvReader(mapping_file)
        for line_dict in tsv_reader:
            result[line_dict['Ensembl Transcript ID']] = line_dict['Ensembl Protein ID']
        return result

    def _determine_protein_seq(self, tx):
        cds_start, cds_stop = tx.determine_cds_footprint()
        if cds_start == -1 or cds_stop == -1:
            return ""
        protein_seq = self._calculate_protein_sequence(tx.get_exons(), tx.get_seq(), cds_start, cds_stop, tx.get_strand())
        protein_seq = ''.join([protein_seq, '*'])
        return protein_seq

    def _calculate_protein_sequence(self, exons, seq, cds_start_genomic_space, cds_stop_genomic_space, strand):
        cds_start_exon_space, cds_stop_exon_space = TranscriptProviderUtils._convert_genomic_space_to_feature_space(int(cds_start_genomic_space), int(cds_stop_genomic_space), exons, strand)

        prot_seq = MutUtils.translate_sequence(seq[int(cds_start_exon_space):int(cds_stop_exon_space)])
        if len(prot_seq) > 0 and prot_seq[-1] == '*':
            prot_seq = prot_seq[:-1]

        return prot_seq

    def _convertGFFRecordToTranscript(self, gff_record, seq_dict, seq_dict_keys, tx_to_protein_mapping):
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

        try:
            tx = self._transcript_index[transcript_id]
        except KeyError:

            # Create the initial record for this transcript.
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
                    self._transcript_index[transcript_id].add_other_attribute(attribute, "|".join(quals[attribute]))

            seq = seq_dict.get(transcript_id, None)

            if seq is not None:
                genome_seq_as_str = str(seq.seq)
            else:
                genome_seq_as_str = ""

            self._transcript_index[transcript_id].set_seq(genome_seq_as_str)

            tx_id_for_protein_lookup = transcript_id
            if '.' in transcript_id:
                tx_id_for_protein_lookup = tx_id_for_protein_lookup[:tx_id_for_protein_lookup.index('.')]
            self._transcript_index[transcript_id].set_protein_id(tx_to_protein_mapping.get(tx_id_for_protein_lookup, ""))

            tx = self._transcript_index[transcript_id]

        gff_type = gff_record['type']
        if gff_type == 'exon':
            tx.add_exon(gff_record['location'][0], gff_record['location'][1], quals['exon_number'][0])
        elif gff_type == 'CDS':
            tx.add_cds(gff_record['location'][0], gff_record['location'][1])
        elif gff_type == 'start_codon':
            tx.set_start_codon(gff_record['location'][0], gff_record['location'][1])
        elif gff_type == 'stop_codon':
            tx.set_stop_codon(gff_record['location'][0], gff_record['location'][1])

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

    def build_ensembl_transcript_index(self, ensembl_input_gtfs, ensembl_input_fastas, output_filename, protocol="file", protein_id_mapping_file=None):
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

        # Get the transcript ID to protein ID mapping
        tx_to_protein_mapping = self._create_tx_id_to_protein_id_mapping(protein_id_mapping_file)

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

                self._convertGFFRecordToTranscript(rec, seq_dict, seq_dict_keys, tx_to_protein_mapping)
                ctr += 1
                if (ctr % 10000) == 0:
                    logging.getLogger(__name__).info("Added " + str(ctr) + " lines of gtf " + str(file_ctr+1) + " of " + str(len(ensembl_input_gtfs)) + " (" + in_file + ") into internal transcript index.")
            in_handle.close()
            logging.getLogger(__name__).info("Finished " + str(ctr) + " lines of gtf (" + in_file + ")")

        logging.getLogger(__name__).info("Populating final db with internal transcript index.")
        transcript_index_keys = self._transcript_index.keys()
        for i,k in enumerate(transcript_index_keys):

            # Populate the protein sequence
            protein_sequence = self._determine_protein_seq(self._transcript_index[k])
            self._transcript_index[k].set_protein_seq(protein_sequence)

            shove[k] = self._transcript_index[k]
            if i % 10000 == 0:
                logging.getLogger(__name__).info("Saved %0.1f%% of transcript index to disk with protein sequence." % (float(i*100)/float(len(transcript_index_keys))))

        logging.getLogger(__name__).info("Transcript index created " + str(len(shove.keys())) + " transcripts: " + protocol + "://" + output_filename)
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

        for i,tx_id in enumerate(transcript_keys):
            tx = transcript_db[tx_id]
            gene = tx.get_gene()
            try:
                tmpList = output_db[gene]
            except KeyError:
                output_db[gene] = []
                tmpList = output_db[gene]
            tmpList.append(tx)
            output_db[gene] = tmpList
            if (i+1) % 10000 == 0:
                logging.getLogger(__name__).info("Gene index added " + str(i) + " transcripts so far.")
        logging.getLogger(__name__).info("Finished gene index with " + str(len(output_db.keys())) + " genes.")
        output_db.close()
        transcript_db.close()

    def build_ensembl_transcripts_by_genomic_location_index(self, ensembl_transcript_index_fname, output_filename, protocol="file"):
        """Create an index for genomic position to transcripts index, using a transcript index created in
            build_ensembl_transcript_index
        """
        transcript_db = Shove(protocol + "://" + ensembl_transcript_index_fname)
        output_db = Shove(protocol + "://" + output_filename, optimize=False)

        transcript_keys = transcript_db.keys()

        for i,tx_id in enumerate(transcript_keys):
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
            if (i+1) % 10000 == 0:
                logging.getLogger(__name__).info("Genomic position index added " + str(i) + " transcripts so far.")

        output_db.close()
        transcript_db.close()

    def construct_ensembl_indices(self, ensembl_input_gtfs, ensembl_input_fastas, base_output_filename, protein_id_mapping_file=None):
        """

        :param ensembl_input_gtfs: (list) gtf input files
        :param ensembl_input_fastas: (list) fasta input files
        :param base_output_filename: Just the base output filename, such as "my_ensembl" without any extensions.
        :return:
        """
        ensembl_transcript_index_filename = base_output_filename + ".transcript.idx"
        self.build_ensembl_transcript_index(ensembl_input_gtfs, ensembl_input_fastas, ensembl_transcript_index_filename, protein_id_mapping_file=protein_id_mapping_file)
        self.build_ensembl_transcripts_by_gene_index(ensembl_transcript_index_filename, base_output_filename + ".transcript_by_gene.idx")
        self.build_ensembl_transcripts_by_genomic_location_index(ensembl_transcript_index_filename, base_output_filename + ".transcript_by_gp_bin.idx")
