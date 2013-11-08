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

    variants = lambda: (
        (("22", "22221645", "22221645", "G", "-"), "Frame_Shift_Del"),
        (("22",	"22221645", "22221645", "-", "A"), "Frame_Shift_Ins")
    )

    def _convert_variant_tuple_to_mut(self, variant_tuple):
        m = MutationData()
        m.chr = variant_tuple[0]
        m.start = variant_tuple[1]
        m.end = variant_tuple[2]
        m.ref_allele = variant_tuple[3]
        m.alt_allele = variant_tuple[4]
        return m

    @data_provider(variants)
    def test_variant_classification(self, variant_tuple, gt_vc):
        gencode_input_gtf = "testdata/gencode/MAPK1.gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/MAPK1.gencode.v18.pc_transcripts.fa"
        base_output_filename = "out/test_variant_classification"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(gencode_input_gtf, gencode_input_fasta, base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, version="TEST")

        m = self._convert_variant_tuple_to_mut(variant_tuple)
        m2 = ensembl_ds.annotate_mutation(m)
        vcer = VariantClassifier()
        recs = ensembl_ds.get_overlapping_transcripts(m.chr, m2.start, m2.end)
        vcer.variant_classify(recs[0], m2['variant_type'], m2.ref_allele, m2.alt_allele, m2.start, m2.end)
        self.assertTrue(gt_vc == "TODO: Replace with annotated version")

    def test_region_queries(self):
        b = region2bin(22221612, 22221919)
        bins = region2bins(22221645, 22221645)
        self.assertTrue(b in bins)

if __name__ == '__main__':
    unittest.main()
