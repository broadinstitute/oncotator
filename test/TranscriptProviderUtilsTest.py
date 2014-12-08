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

import shutil
import Bio
from oncotator.DatasourceFactory import DatasourceFactory
import unittest
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.utils.VariantClassification import VariantClassification
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from test.TestUtils import TestUtils, data_provider_decorator

TestUtils.setupLogging(__file__, __name__)
class TranscriptProviderUtilsTest(unittest.TestCase):

    _multiprocess_can_split_ = True

    def test_convert_genomic_space_to_transcript_space(self):
        base_config_location = "testdata/ensembl/saccer/"
        ensembl_ds = DatasourceFactory.createDatasource(base_config_location + "ensembl.config", base_config_location)

        tx = ensembl_ds.get_overlapping_transcripts("I", "350", "350") # transcript starts at 335.
        start, end = TranscriptProviderUtils.convert_genomic_space_to_transcript_space("350", "350", tx[0])
        self.assertTrue(start == end)
        self.assertTrue(start == 16)

        tx = ensembl_ds.get_overlapping_transcripts("II", "764690", "764690") # transcript starts at 764697 (strand is '-').
        start, end = TranscriptProviderUtils.convert_genomic_space_to_transcript_space("764690", "764690", tx[0])
        self.assertTrue(start == end)
        self.assertTrue(start == 7)

        start, end = TranscriptProviderUtils.convert_genomic_space_to_transcript_space("764680", "764690", tx[0])
        self.assertTrue(start == (end - 10))
        self.assertTrue(start == 7)

    locs = lambda: (
        (("22108790", "22108790"), 11020), (("22108800", "22108890"), 10920)
    )

    @data_provider_decorator(locs)
    def test_convert_genomic_space_to_exon_space(self, loc, gt_d):
        """Test genomic --> exon transform on real data. """
        gencode_input_gtf = "testdata/gencode/MAPK1.gencode.v19.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/MAPK1.gencode.v19.pc_transcripts.fa"
        base_output_filename = "out/test_variant_classification"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([gencode_input_gtf], [gencode_input_fasta], base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, version="TEST")
        tx = ensembl_ds.get_overlapping_transcripts("22", "22108790", "22108790")

        start, end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(loc[0], loc[1], tx[0])
        loc_length = (int(loc[1]) - int(loc[0]))
        self.assertTrue((end - start) == loc_length, str(end) + " - " + str(start) + " was not correct length: " + str(loc_length))
        self.assertTrue(start == gt_d, "start position (" + str(start) + ") did not match gt (" + str(end) + ")" + "   exons: " + str(tx[0].get_exons()))

    simple_locs = lambda: (
        ([(5, 10, 1), (15, 20, 2), (25, 40, 3)], 6, 1, "+"),
        ([(25, 40, 1), (15, 20, 2), (5, 10, 3)], 6, 24, "-"),
        ([(22221611, 22221919, '1'), (22161952, 22162135, '2'), (22160138, 22160328, '3'), (22153300, 22153417, '4'), (22142982, 22143097, '5'), (22142545, 22142677, '6'), (22127161, 22127271, '7'), (22123483, 22123609, '8'), (22108788, 22118529, '9')],
            22162133, 308 + 2, "-"),
        ([(22221611, 22221919, '1'), (22161952, 22162135, '2'), (22160138, 22160328, '3'), (22153300, 22153417, '4'), (22142982, 22143097, '5'), (22142545, 22142677, '6'), (22127161, 22127271, '7'), (22123483, 22123609, '8'), (22108788, 22118529, '9')],
            22108800, 11022-12, "-"),
        ([(22221611, 22221919, '1'), (22161952, 22162135, '2'), (22160138, 22160328, '3'), (22153300, 22153417, '4'), (22142982, 22143097, '5'), (22142545, 22142677, '6'), (22127161, 22127271, '7'), (22123483, 22123609, '8'), (22108788, 22118529, '9')],
            22108890, 11022-102, "-"),
    )

    @data_provider_decorator(simple_locs)
    def test_transform_to_feature_space(self, exons, s, gt, strand):
        """Run some basic tests transforming genomic coordinates to exon coordinates, taking strand into account. """

        guess = TranscriptProviderUtils._transform_to_feature_space(exons, s, strand)
        self.assertTrue(guess == gt, "Did not transform genomic to exon space properly: " + str(exons) +  "   pos: " + str(s) + "  strand: " + strand + "  guess/gt: " + str(guess) + "/" + str(gt))

    def _create_ensembl_ds_from_testdata(self, gene):
        gencode_input_gtf = "testdata/gencode/" + gene + ".gencode.v19.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/" + gene + ".gencode.v19.pc_transcripts.fa"
        base_output_filename = "out/test_variant_classification"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([gencode_input_gtf], [gencode_input_fasta], base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, title="GENCODE", version="v18")
        return ensembl_ds

    def retrieve_test_transcript_MUC16(self):
        ensembl_ds = self._create_ensembl_ds_from_testdata("MUC16")
        tx = ensembl_ds.transcript_db['ENST00000397910.4']
        self.assertTrue(tx is not None, "Unit test appears to be misconfigured or a bug exists in the ensembl datasource code.")
        return tx

    def retrieve_test_transcript_MAPK1(self):
        ensembl_ds = self._create_ensembl_ds_from_testdata("MAPK1")
        tx = ensembl_ds.transcript_db['ENST00000215832.6']
        self.assertTrue(tx is not None, "Unit test appears to be misconfigured or a bug exists in the ensembl datasource code.")
        return tx

    # These are stranded, so since the test transcript is "-", these are reverse complement of what you would
    #  see in genome browser.
    seq_testdata = lambda: (
        ("22143048", "22143050", "AGA"),
        ("22143050", "22143050", "A"),
        ("22143044", "22143046", "ATG"),
        ("22108789", "22108795", "CTATAAA")
    )
    @data_provider_decorator(seq_testdata)
    def test_seq(self, start, end, gt):
        """Test that we can successfully determine the codon at an arbitrary location on test transcript"""
        tx = self.retrieve_test_transcript_MAPK1()

        transcript_position_start, transcript_position_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(start, end, tx)
        transcript_seq = tx.get_seq()
        seq = transcript_seq[transcript_position_start:transcript_position_end+1]
        self.assertTrue(seq == gt, "Incorrect seq found guess,gt (%s, %s)" %(seq, gt))

    # start, end, base (stranded), gt_codon (stranded)
    codon_tests_single_base = lambda: (

        ("22127164", "22127164", "G", "GAG"),
        ("22127165", "22127165", "C", "GAC")
    )

    @data_provider_decorator(codon_tests_single_base)
    def test_codon_single_base(self, start, end, ref_base_stranded, gt_codon):
        """Test that we can grab the proper three bases of a codon for an arbitrary single base """
        tx = self.retrieve_test_transcript_MAPK1()
        transcript_position_start, transcript_position_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(start, end, tx)
        cds_start, cds_stop = TranscriptProviderUtils.determine_cds_in_exon_space(tx)
        protein_position_start, protein_position_end = TranscriptProviderUtils.get_protein_positions(transcript_position_start, transcript_position_end, cds_start)
        cds_codon_start, cds_codon_end = TranscriptProviderUtils.get_cds_codon_positions(protein_position_start, protein_position_end, cds_start)

        codon_seq = tx.get_seq()[cds_codon_start:cds_codon_end+1]
        self.assertTrue(codon_seq == gt_codon, "Did not get correct codon (%s): %s    loc: %s-%s" %(gt_codon, codon_seq, start, end))

    transcript_change_testdata = lambda: (
        ("SNP", VariantClassification.MISSENSE, 6353, 6353, "C", "T", "c.6353C>T", ""),
        ("SNP", VariantClassification.NONSENSE, 6037, 6037, "C", "T", "c.6037C>T", ""),
        ("SNP", VariantClassification.MISSENSE, 192, 192, "G", "A", "c.192G>A", ""),
        ("SNP", VariantClassification.SPLICE_SITE, 316, 316, "G", "A", "c.316_splice", VariantClassification.INTRON),
        ("SNP", VariantClassification.SPLICE_SITE, 316, 316, "G", "A", "c.316G>A", VariantClassification.MISSENSE),
        ("SNP", VariantClassification.SPLICE_SITE, 2784, 2784, "C", "A", "c.2784C>A", VariantClassification.MISSENSE),
        ("DEL", VariantClassification.IN_FRAME_DEL, 1358, 1360, "AGA", "-", "c.1358_1360delAGA", ""),
        ("DEL", VariantClassification.IN_FRAME_DEL, 1358, 1360, "AGA", "", "c.1358_1360delAGA", ""),
        ("SNP", VariantClassification.RNA,	2543, 2543, "A", "G", "c.2543A>G", ""),
        ("INS", VariantClassification.FRAME_SHIFT_INS, 3990, 3991, "-", "TTCTTAAG", "c.3990_3991insTTCTTAAG", ""),
        ("INS", VariantClassification.FRAME_SHIFT_INS, 3990, 3991, "", "TTCTTAAG", "c.3990_3991insTTCTTAAG", ""),
        ("INS", VariantClassification.FRAME_SHIFT_INS, 3990, 3997, "-", "TTCTTAAG", "c.3990_3997insTTCTTAAG", ""),
        ("INS", VariantClassification.FRAME_SHIFT_INS, 3990, 3997, "", "TTCTTAAG", "c.3990_3997insTTCTTAAG", "")
    )
    @data_provider_decorator(transcript_change_testdata)
    def test_render_transcript_change(self, variant_type, vc, exon_position_start, exon_position_end, ref_allele_stranded, alt_allele_stranded, gt, secondary_vc):
        """Simple test of transcript change, once parameters have been rendered. """
        guess = TranscriptProviderUtils.render_transcript_change(variant_type, vc, exon_position_start, exon_position_end, ref_allele_stranded, alt_allele_stranded, secondary_vc)
        self.assertTrue(guess == gt, "Incorrect guess gt <> guess: %s <> %s" % (gt, guess))

    protein_change_testdata = lambda: (
        ("SNP", "Missense_Mutation", "", 2118, 2118, "S", "L", "-", "p.S2118L"),
        ("SNP", "Nonsense_Mutation",  "", 2013, 2013, "Q", "*", "-", "p.Q2013*"),
        ("SNP", "Splice_Site",  VariantClassification.MISSENSE, 106, 106, "V", "A", "-", "p.V106A"),
        ("SNP", "Splice_Site", VariantClassification.INTRON, 0, 0, "V", "A", "-", ""),
        ("DEL", "In_Frame_Del", "", 454, 454, "K", "-",	"+", "p.K454del"),
        ("SNP", "Nonstop_Mutation", "", 246, 246, "*", "S", "+", "p.*246S")
    )
    @data_provider_decorator(protein_change_testdata)
    def test_render_protein_change(self, variant_type, variant_classification, secondary_vc, prot_position_start, prot_position_end, ref_prot_allele, alt_prot_allele, strand, gt):
        """Simple test of protein change, once parameters have been rendered. """
        guess = TranscriptProviderUtils.render_protein_change(variant_type, variant_classification, prot_position_start, prot_position_end, ref_prot_allele, alt_prot_allele, secondary_vc)
        self.assertTrue(guess == gt, "Incorrect guess gt <> guess: %s <> %s" % (gt, guess))

    mutate_ref_sequence_testdata = lambda: (
        ("DEL", 22221919, 22221919, "T", "-", 0, 0, ""),
        ("INS", 22221919, 22221919, "-", "G", 0, 2, "ACGG"),
        ("INS", 22221919, 22221920, "-", "T", 0, 1, "AAG"),
        ("SNP", 22221919, 22221919, "T", "G", 0, 0, "C"),
        ("DEL", 22221919, 22221919, "T", "-", 0, 0, ""),
    )
    @data_provider_decorator(mutate_ref_sequence_testdata)
    def test_mutate_reference_seqeunce(self, vt, start, end, ref, alt, start_exon_space, end_exon_space, mutated_seq_gt):
        """ Test that we can render a mutated sequence with SNP, INS, and DEL
        """
        # mutated_seq_gt is stranded and this is a "-" transcript
        tx = self.retrieve_test_transcript_MAPK1()
        observed_allele = Bio.Seq.reverse_complement(alt)
        mutated_allele = TranscriptProviderUtils.mutate_reference_sequence(tx.get_seq()[start_exon_space : end_exon_space+1], start_exon_space, start_exon_space, end_exon_space, observed_allele, vt)
        self.assertTrue(mutated_seq_gt == mutated_allele, "No match (gt/guess)  %s/%s for %s." % (mutated_seq_gt, mutated_allele, str([vt, start, end, ref, alt, start_exon_space, end_exon_space, mutated_seq_gt])))

    def test_determine_closest_distance_from_exon_in_exon(self):
        tx = self.retrieve_test_transcript_MAPK1()

        # Right in exon 1
        left_diff, right_diff = TranscriptProviderUtils.determine_closest_distance_from_exon(22162000, 22162005, 1,  tx)
        self.assertTrue(left_diff < 0 and right_diff > 0, "left distance should be negative while right distance should be positive.")

    variant_type_examples = lambda: (
        ("A","G", VariantClassification.VT_SNP),
        ("TG","-", VariantClassification.VT_DEL),
        ("TC","CT", VariantClassification.VT_DNP),
        ("","A", VariantClassification.VT_INS),
        ("ACA","TCG", VariantClassification.VT_TNP),
        ("TTTTT","AAAAA", VariantClassification.VT_ONP)
    )
    @data_provider_decorator(variant_type_examples)
    def test_infer_variant_type(self,ref,alt,vt_gt):
        """test that we can tell a snp from an indel"""
        self.assertEqual(TranscriptProviderUtils.infer_variant_type(ref, alt), vt_gt)

    def test_variant_type_of_non_variant_is_Exception(self):
        """check that empty variants raise an exception"""
        self.assertRaises(Exception, TranscriptProviderUtils.infer_variant_type, "-","-")
        self.assertRaises(Exception, TranscriptProviderUtils.infer_variant_type, "","")

if __name__ == '__main__':
    unittest.main()
