import shutil
from oncotator.DatasourceCreator import DatasourceCreator
import unittest
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources import EnsemblTranscriptDatasource
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from unittest_data_provider import data_provider
from test.TestUtils import TestUtils

TestUtils.setupLogging(__file__, __name__)
class TranscriptProviderUtilsTest(unittest.TestCase):

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

    @data_provider(locs)
    def test_convert_genomic_space_to_exon_space(self, loc, gt_d):
        """Test genomic --> exon transform on real data. """
        gencode_input_gtf = "testdata/gencode/MAPK1.gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/MAPK1.gencode.v18.pc_transcripts.fa"
        base_output_filename = "out/test_variant_classification"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(gencode_input_gtf, gencode_input_fasta, base_output_filename)
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

    @data_provider(simple_locs)
    def test_transform_to_feature_space(self, exons, s, gt, strand):
        """Run some basic tests transforming genomic coordinates to exon coordinates, taking strand into account. """

        guess = TranscriptProviderUtils._transform_to_feature_space(exons, s, strand)
        self.assertTrue(guess == gt, "Did not transform genomic to exon space properly: " + str(exons) +  "   pos: " + str(s) + "  strand: " + strand + "  guess/gt: " + str(guess) + "/" + str(gt))


if __name__ == '__main__':
    unittest.main()
