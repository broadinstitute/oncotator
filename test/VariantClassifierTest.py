import shutil
from oncotator.DatasourceCreator import DatasourceCreator
from oncotator.MutationData import MutationData
from oncotator.datasources import EnsemblTranscriptDatasource
from oncotator.index.gaf import region2bin, region2bins
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from test.TestUtils import TestUtils
from MUC16Testdata import muc16testdata
__author__ = 'lichtens'

import unittest
from unittest_data_provider import data_provider
from oncotator.utils.VariantClassifier import VariantClassifier

TestUtils.setupLogging(__file__, __name__)
class VariantClassifierTest(unittest.TestCase):
    _multiprocess_can_split_ = True
    variants_indels = lambda: (
        ("22", "22221645", "22221647", "In_Frame_Del", "DEL", "GAG", "-"),
        ("22", "22221645", "22221645", "Frame_Shift_Del", "DEL", "G", "-"),
        ("22",	"22221645", "22221645", "Frame_Shift_Ins", "INS", "-", "A")
    )
    # TODO: Get recently downloaded test data and use that.

    def _create_ensembl_ds_from_testdata(self, gene):
        gencode_input_gtf = "testdata/gencode/" + gene + ".gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/" + gene + ".gencode.v18.pc_transcripts.fa"
        base_output_filename = "out/test_variant_classification"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(gencode_input_gtf, gencode_input_fasta, base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, title="GENCODE", version="v18")
        return ensembl_ds

    def _test_variant_classification(self, alt, chr, end, gt_vc, ref, start, vt, gene="MAPK1"):
        ensembl_ds = self._create_ensembl_ds_from_testdata(gene)
        recs = ensembl_ds.get_overlapping_transcripts(chr, start, end)
        tx = ensembl_ds._choose_transcript(recs, EnsemblTranscriptDatasource.TX_MODE_CANONICAL, vt, ref, alt, start, end)
        self.assertTrue(len(recs) != 0, "Issue with test...No transcripts found for: " + str([chr, start, end]))

        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, vt, ref, alt, start, end)
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))

    @data_provider(variants_indels)
    def test_variant_classification_indels(self, chr, start, end, gt_vc, vt, ref, alt):
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)

    variants_snps_missense = lambda: (
        ("22", "22127164", "22127164", "Missense_Mutation", "SNP", "C", "T"),
        ("22", "22143049", "22143049", "Missense_Mutation", "SNP", "C", "T"),
        ("22", "22162014", "22162014", "Missense_Mutation", "SNP", "C", "T")
    )
    @data_provider(variants_snps_missense)
    def test_variant_classification_missense(self, chr, start, end, gt_vc, vt, ref, alt):
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)

    def test_region_queries(self):
        b = region2bin(22221612, 22221919)
        bins = region2bins(22221645, 22221645)
        self.assertTrue(b in bins)

    frameshift_indels = lambda : (
        ("INS", 10, 11,  "A", True),
        ("INS", 10, 12,  "ATC", False),
        ("DEL", 10, 10,  "-", True),
        ("DEL", 10, 12,  "-", False)
    )

    @data_provider(frameshift_indels)
    def test_is_framshift_indel(self, vt, s, e, alt, gt):
        vcer = VariantClassifier()
        self.assertTrue(vcer.is_framshift_indel(vt, s, e, alt) == gt)

    exon_tests = lambda: (
        ([(22221611, 22221919, '1'), (22161952, 22162135, '2'), (22160138, 22160328, '3'), (22153300, 22153417, '4'), (22142982, 22143097, '5'), (22142545, 22142677, '6'), (22127161, 22127271, '7'), (22123483, 22123609, '8'), (22108788, 22118529, '9')],
            22108800, 22108800, 8, 12, 9729),
        ([(5, 10, 1), (15, 20, 2), (25, 40, 3)], 6, 8, 0, 1, 2),
        ([(25, 40, 1), (15, 20, 2), (5, 10, 3)], 14, 16, 1, 1, 1)
    )

    @data_provider(exon_tests)
    def test_determine_closest_distances_from_exon(self, exons, start_genomic_space, end_genomic_space, gt_exon_i, gt_exon_ld, gt_exon_rd):
        vcer = VariantClassifier()
        guess_exon_i, guess_exon_ld, guess_exon_rd = vcer._determine_closest_distances_from_exon(exons, start_genomic_space, end_genomic_space)
        self.assertTrue(guess_exon_i == gt_exon_i, "guess, gt exon did not match: " + str([guess_exon_i, gt_exon_i]))
        self.assertTrue(gt_exon_ld == guess_exon_ld, "guess, gt ld did not match: " + str([guess_exon_ld, gt_exon_ld]))
        self.assertTrue(gt_exon_rd == guess_exon_rd, "guess, gt rd did not match: " + str([guess_exon_rd, gt_exon_rd]))

    @data_provider(muc16testdata)
    def test_muc16_snps(self, chr, start, end, gt_vc, vt, ref, alt):
        """ Test all of the MUC16 SNPs."""
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, gene="MUC16")

    def test_snp_vc_on_one_transcript_5UTR(self):
        """Take the test transcript (ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand) and test the entire 5'UTR"""
        tx = self.retrieve_test_transcript()
        vcer = VariantClassifier()

        chr = 22
        for i in range(0, 200):
            start = end = (22221919 - i)
            ref = tx.get_seq()[i]
            alt = 'G'
            vc = vcer.variant_classify(tx, "SNP", ref, alt, start, end)
            if i < 189:
                self.assertTrue(vc == "5'UTR", "Should be 5'UTR, but saw " + vc + ".  For " + str([ref, alt, start, end]))
            if i < 200 and i >= 189:
                self.assertTrue(vc != "5'UTR", "Should not be 5'UTR, but saw " + vc + ".  For " + str([ref, alt, start, end]))

    def retrieve_test_transcript(self):
        ensembl_ds = self._create_ensembl_ds_from_testdata("MAPK1")
        tx = ensembl_ds.transcript_db['ENST00000215832.6']
        self.assertTrue(tx is not None,
                        "Unit test appears to be misconfigured or a bug exists in the ensembl datasource code.")
        return tx

    variants_snps_splice_sites = lambda: (
        ("22", "22162138", "22162138", "Intron", "SNP", "T", "G"),
        ("22", "22162137", "22162137", "Splice_Site", "SNP", "T", "T"),
        ("22", "22162136", "22162136", "Splice_Site", "SNP", "C", "T"),
        ("22", "22162135", "22162135", "Splice_Site", "SNP", "G", "T"),
        ("22", "22162134", "22162134", "Splice_Site", "SNP", "A", "T"),
        ("22", "22162133", "22162133", "Missense_Mutation", "SNP", "G", "T"),
        ("22", "22162132", "22162132", "Missense_Mutation", "SNP", "A", "T")
    )
    @data_provider(variants_snps_splice_sites)
    def test_snp_vc_on_one_transcript_Splice_Site(self, chr, start, end, gt_vc, vt, ref,alt):
        """Test a lot of positions on one MAPK1 transcript: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self.retrieve_test_transcript()

        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, vt, ref, alt, start, end)
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))


if __name__ == '__main__':
    unittest.main()
