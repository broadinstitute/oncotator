import shutil
from oncotator.DatasourceCreator import DatasourceCreator
from oncotator.MutationData import MutationData
from oncotator.datasources import EnsemblTranscriptDatasource
from oncotator.index.gaf import region2bin, region2bins
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from test.TestUtils import TestUtils

__author__ = 'lichtens'

import unittest
from unittest_data_provider import data_provider
from oncotator.utils.VariantClassifier import VariantClassifier

TestUtils.setupLogging(__file__, __name__)
class VariantClassifierTest(unittest.TestCase):

    variants_indels = lambda: (
        ("22", "22221645", "22221647", "In_Frame_Del", "DEL", "GAG", "-"),
        ("22", "22221645", "22221645", "Frame_Shift_Del", "DEL", "G", "-"),
        ("22",	"22221645", "22221645", "Frame_Shift_Ins", "INS", "-", "A")
    )
    # TODO: Get recently downloaded test data and use that.

    @data_provider(variants_indels)
    def test_variant_classification_indels(self, chr, start, end, gt_vc, vt, ref, alt):
        gencode_input_gtf = "testdata/gencode/MAPK1.gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/MAPK1.gencode.v18.pc_transcripts.fa"
        base_output_filename = "out/test_variant_classification"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(gencode_input_gtf, gencode_input_fasta, base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, version="TEST")

        vcer = VariantClassifier()
        recs = ensembl_ds.get_overlapping_transcripts(chr, start, end)
        vc = vcer.variant_classify(recs[0], vt, ref, alt, start, end)
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc)

    def test_region_queries(self):
        b = region2bin(22221612, 22221919)
        bins = region2bins(22221645, 22221645)
        self.assertTrue(b in bins)

    def test_get_protein_sequence(self):

        gencode_input_gtf = "testdata/gencode/MAPK1.gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/MAPK1.gencode.v18.pc_transcripts.fa"
        base_output_filename = "out/test_get_protein_sequence"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(gencode_input_gtf, gencode_input_fasta, base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, version="TEST")

        # Never do this in real code outside of the datasource
        tx = ensembl_ds.transcript_db["ENST00000215832.6"]
        vcer = VariantClassifier()
        gt_seq = "ATGGCGGCGGCGGC"
        gt_prot = "MAAA"
        protein_seq = vcer.get_protein_sequence(tx, 22221730 - len(gt_seq), 22221730)
        self.assertIsNotNone(protein_seq)
        self.assertEqual(gt_prot, protein_seq)

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


if __name__ == '__main__':
    unittest.main()
