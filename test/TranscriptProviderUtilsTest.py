import shutil
from oncotator.DatasourceCreator import DatasourceCreator
import unittest
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources import EnsemblTranscriptDatasource
from oncotator.utils.VariantClassification import VariantClassification
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from test.TestUtils import TestUtils, data_provider_decorator

TestUtils.setupLogging(__file__, __name__)
class TranscriptProviderUtilsTest(unittest.TestCase):

    _multiprocess_can_split_ = True

    def test_convert_genomic_space_to_transcript_space(self):
        base_config_location = "testdata/ensembl/saccer/"
        ensembl_ds = DatasourceCreator.createDatasource(base_config_location + "ensembl.config", base_config_location)

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
        gencode_input_gtf = "testdata/gencode/MAPK1.gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/MAPK1.gencode.v18.pc_transcripts.fa"
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
        gencode_input_gtf = "testdata/gencode/" + gene + ".gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/" + gene + ".gencode.v18.pc_transcripts.fa"
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
        ("SNP", VariantClassification.MISSENSE, 6353, 6353, "C", "T", "c.6353C>T"),
        ("SNP", VariantClassification.NONSENSE, 6037, 6037, "C", "T", "c.6037C>T"),
        ("SNP", VariantClassification.MISSENSE, 192, 192, "G", "A", "c.192G>A"),
        ("SNP", VariantClassification.SPLICE_SITE, 316, 316, "G", "A", "c.316_splice"),
        ("DEL", VariantClassification.IN_FRAME_DEL, 1358, 1360, "AGA", "-", "c.1358_1360delAGA"),
        ("DEL", VariantClassification.IN_FRAME_DEL, 1358, 1360, "AGA", "", "c.1358_1360delAGA"),
        ("SNP", VariantClassification.RNA,	2543, 2543, "A", "G", "c.2543A>G"),
        ("INS", VariantClassification.FRAME_SHIFT_INS, 3990, 3991, "-", "TTCTTAAG", "c.3990_3991insTTCTTAAG"),
        ("INS", VariantClassification.FRAME_SHIFT_INS, 3990, 3991, "", "TTCTTAAG", "c.3990_3991insTTCTTAAG"),
        ("INS", VariantClassification.FRAME_SHIFT_INS, 3990, 3997, "-", "TTCTTAAG", "c.3990_3997insTTCTTAAG"),
        ("INS", VariantClassification.FRAME_SHIFT_INS, 3990, 3997, "", "TTCTTAAG", "c.3990_3997insTTCTTAAG")
    )
    @data_provider_decorator(transcript_change_testdata)
    def test_render_transcript_change(self, variant_type, vc, exon_position_start, exon_position_end, ref_allele_stranded, alt_allele_stranded, gt):
        """Simple test of transcript change, once parameters have been rendered. """
        guess = TranscriptProviderUtils.render_transcript_change(variant_type, vc, exon_position_start, exon_position_end, ref_allele_stranded, alt_allele_stranded)
        self.assertTrue(guess == gt, "Incorrect guess gt <> guess: %s <> %s" % (gt, guess))

    protein_change_testdata = lambda: (
        ("SNP", "Missense_Mutation", 2118, 2118, "S", "L", "-", "p.S2118L"),
        ("SNP", "Nonsense_Mutation", 2013, 2013, "Q", "*", "-", "p.Q2013*"),
        ("SNP", "Splice_Site", 106, 106, "V", "A", "-", "p.V106_splice"),
        ("DEL", "In_Frame_Del", 454, 454, "K", "-",	"+", "p.K454del"),
        ("SNP", "Nonstop_Mutation", 246, 246, "*", "S", "+", "p.*246S")
    )
    @data_provider_decorator(protein_change_testdata)
    def test_render_protein_change(self, variant_type, variant_classification, prot_position_start, prot_position_end, ref_allele, alt_allele, strand, gt):
        """Simple test of protein change, once parameters have been rendered. """
        guess = TranscriptProviderUtils.render_protein_change(variant_type, variant_classification, prot_position_start, prot_position_end, ref_allele, alt_allele)
        self.assertTrue(guess == gt, "Incorrect guess gt <> guess: %s <> %s" % (gt, guess))

if __name__ == '__main__':
    unittest.main()
