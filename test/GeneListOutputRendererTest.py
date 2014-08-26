import unittest
import os
from oncotator.Annotator import Annotator
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.utils.RunSpecificationFactory import RunSpecificationFactory
from oncotator.utils.RunSpecification import RunSpecification
from test.TestUtils import TestUtils


class GeneListOutputRendererTest(unittest.TestCase):
    _multiprocess_can_split_ = True
    def setUp(self):
        self.config = TestUtils.createUnitTestConfig()
        pass

    @TestUtils.requiresDefaultDB()
    def test_basic_rendering(self):
        """Test that we can render a basic seg file as a gene list"""
        inputFilename = "testdata/seg/Patient0.seg.txt"
        output_filename = "out/test_basic_rendering.gene_list.tsv"
        db_dir = self.config.get('DEFAULT',"dbDir")
        if os.path.exists(output_filename):
            os.remove(output_filename)

        annotator = Annotator()
        run_spec = RunSpecificationFactory.create_run_spec("SEG_FILE", "GENE_LIST", inputFilename, output_filename,
                                                           datasourceDir=db_dir, annotating_type=RunSpecification.ANNOTATE_SEGMENTS)
        annotator.initialize(run_spec)
        annotator.annotate()

        # Now check the output
        output_reader = GenericTsvReader(output_filename)

        headers = output_reader.getFieldNames()

        for line_dict in output_reader:
            self.assertTrue(line_dict['segment_start'] is not None)
            self.assertTrue(line_dict['segment_start'].strip() != "")
            self.assertTrue(line_dict['segment_end'] is not None)
            self.assertTrue(line_dict['segment_end'].strip() != "")
            self.assertTrue("gene" in line_dict.keys())
            self.assertTrue(len(line_dict["gene"]) > 0)
            self.assertTrue(float(line_dict["segment_num_probes"]))
            self.assertTrue(line_dict['sample'] == "Patient0")

    @TestUtils.requiresDefaultDB()
    def test_rendering_with_exons(self):
        """Test that we can render a seg file that includes exons at end points"""
        inputFilename = "testdata/seg/Middle_of_exon.seg.txt"
        output_filename = "out/test_exon_seg2.gene_list.tsv"
        db_dir = self.config.get('DEFAULT',"dbDir")
        if os.path.exists(output_filename):
            os.remove(output_filename)

        annotator = Annotator()
        run_spec = RunSpecificationFactory.create_run_spec("SEG_FILE", "GENE_LIST", inputFilename, output_filename,
                                                           datasourceDir=db_dir, annotating_type=RunSpecification.ANNOTATE_SEGMENTS)
        annotator.initialize(run_spec)
        annotator.annotate()

        # Now check the output
        output_reader = GenericTsvReader(output_filename)

        headers = output_reader.getFieldNames()

        for line_dict in output_reader:
            self.assertTrue(line_dict['segment_start'] is not None)
            self.assertTrue(line_dict['segment_start'].strip() != "")
            if line_dict['segment_end_gene'] == "MAPK1":
                self.assertTrue(line_dict['segment_end_exon'].strip() == "8+", "Should have been 8+, but saw: %s" % line_dict['segment_end_exon'].strip())

if __name__ == '__main__':
    unittest.main()
