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

        required_cols = ["Sample", "Num_Probes", "Segment_Mean"]
        headers = output_reader.getFieldNames()

        for line_dict in output_reader:
            self.assertTrue(line_dict['start'] is not None)
            self.assertTrue(line_dict['start'].strip() != "")
            self.assertTrue(line_dict['end'] is not None)
            self.assertTrue(line_dict['end'].strip() != "")
            self.assertTrue("genes" in line_dict.keys())
            self.assertTrue(len(line_dict["genes"].split(",")) > 0)



if __name__ == '__main__':
    unittest.main()
