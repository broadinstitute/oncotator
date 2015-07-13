
import unittest
from oncotator.input.InputMutationCreator import InputMutationCreator
from oncotator.output.OutputRenderer import OutputRenderer
from test.TestUtils import TestUtils
from oncotator.utils.RunSpecificationFactory import RunSpecificationFactory


TestUtils.setupLogging(__file__, __name__)

class RunSpecificationFactoryTest(unittest.TestCase):

    def test_run_spec_creation_no_datasources(self):
        """Test that we can create a run spec with no datasources"""
        run_spec = RunSpecificationFactory.create_run_spec_given_datasources(input_format="VCF",
                                                                             input_filename="testdata/m2_support/phasingExample.vcf",
                                                                    output_format="TCGAMAF",
                                                                    output_filename="out/foo.maf.annotated",
                                                                    datasource_list=[])
        self.assertTrue(isinstance(run_spec.inputCreator, InputMutationCreator))
        self.assertTrue(isinstance(run_spec.outputRenderer, OutputRenderer))
        self.assertTrue(run_spec.is_allow_annotation_overwriting==False)

    def test_tcgamaf_invalid_input_file(self):
        """Test a case where TCGAMAF specified as input and we get an error (as we should) for a missing file"""
        is_exception_seen = False
        try:
            run_spec = RunSpecificationFactory.create_run_spec_given_datasources(input_format="TCGAMAF",
                                                                             input_filename="testdata/Idonotexist",
                                                                    output_format="TCGAMAF",
                                                                    output_filename="out/foo.maf.annotated",
                                                                    datasource_list=[])
        except IOError as ie:
            is_exception_seen = True

        self.assertTrue(is_exception_seen)
