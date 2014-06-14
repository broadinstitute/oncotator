# LICENSE_GOES_HERE
from TestUtils import TestUtils
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MutationData import MutationData
from oncotator.utils.OncotatorCLIUtils import OncotatorCLIUtils, RunSpecification

"""
Created on Nov 7, 2012

@author: lichtens
"""
import unittest
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.Annotator import Annotator
import logging
import os
from oncotator.utils.GenericTsvReader import GenericTsvReader

TestUtils.setupLogging(__file__, __name__)


class AnnotatorTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def testBlankAnnotatorInit(self):
        """ Test an extremely simple scenario, where no additional annotations are needed.  I.e. no data sources """
        self.logger.info("Starting Blank Annotator Init Test...")

        inputCreator = MafliteInputMutationCreator('testdata/maflite/tiny_maflite.maf.txt')
        outputRenderer = SimpleOutputRenderer("out/testBlankAnnotatorTestFile.tsv")

        # Assumed myIC and myOC have been initialized as the proper Input and Output Creators, respectively.
        # 1) Initialize the Annotator
        annotator = Annotator()
        annotator.setInputCreator(inputCreator)
        annotator.setOutputRenderer(outputRenderer)
        testOutputFilename = annotator.annotate()

        # Test that file exists and that it has correct # of mutations (+1 for header +1 for annotator comment line).
        numSamples = 1
        numExtraLines = 3  # one for header, two for comment lines
        numDoubleLines = 0  # Number of lines with two alt alleles
        numVariants = 9
        gt = numSamples * numVariants + numDoubleLines * numSamples + numExtraLines
        fp = file(testOutputFilename, 'r')
        ctr = 0
        for line in fp:
            ctr += 1
        fp.close()
        self.assertEqual(ctr, gt,
                         "Number of lines read was not correct: " + str(ctr) + " -- should have been: " + str(gt))

    def testVersionHeader(self):
        """ This method simply tests that the version string returned by the annotator does not cause an exception.
            Minimal checking that the returned sting is actually correct.
            Does not attempt to initialize input or output.  Only a gaf datasource.
         """
        annotator = Annotator()
        annotator.addDatasource(TestUtils.createTranscriptProviderDatasource(self.config))
        tmp = annotator.createHeaderString()
        self.assertTrue(tmp.find("Gaf ") != -1 or tmp.find("GENCODE") != -1, "Could not find Gaf or GENCODE version in header string.")
        self.assertTrue(tmp.find("Oncotator") != -1, "Could not find the word Oncotator in header string.")

    def testManualAnnotations(self):
        """ Test that the manual annotation facility in the Annotator is working properly. """
        annotator = Annotator()
        overrides = {'source': 'Capture', 'status': 'Somatic', 'phase': 'Phase_I', 'sequencer': 'Illumina GAIIx'}
        annotator.setManualAnnotations(overrides)
        inputCreator = MafliteInputMutationCreator('testdata/maflite/Patient0.snp.maf.txt')
        outputRenderer = SimpleOutputRenderer("out/testManualAnnotationsFile.tsv")
        annotator.setInputCreator(inputCreator)
        annotator.setOutputRenderer(outputRenderer)

        testOutputFilename = annotator.annotate()

        keysOfInterest = overrides.keys()

        statinfo = os.stat(testOutputFilename)
        self.assertTrue(statinfo.st_size > 0, "Generated TSV file (" + testOutputFilename + ") is empty.")

        tsvReader = GenericTsvReader(testOutputFilename)

        ctr = 1
        for lineDict in tsvReader:
            for k in keysOfInterest:
                self.assertTrue(lineDict[k] != "__UNKNOWN__",
                                "__UNKNOWN__ value seen on line " + str(ctr) + ", when it should be populated: " + k)
                self.assertTrue(lineDict[k] != "",
                                "Blank value seen on line " + str(ctr) + ", when it should be populated: " + k)
                self.assertTrue(lineDict[k] == overrides[k],
                                "Value for " + k + " on line " + str(ctr) + " did not match override: " + str(
                                    lineDict[k]) + " <> " + str(overrides[k]))
            ctr += 1

    def testDefaultAnnotations(self):
        """Test that the default annotation values populate properly. """
        annotator = Annotator()
        default_annotations = {"test2": "foo2", "test3": "Should not be seen"}
        overrides = {'test3': 'foo3'}

        m1 = MutationData()
        m1.createAnnotation("test1", "foo1")
        m1.createAnnotation("test2", "")

        m2 = MutationData()
        m2.createAnnotation("test1", "")


        m3 = MutationData()
        m3.createAnnotation("test1", "")
        m3.createAnnotation("test2", "foo2-original")

        muts = [m1, m2, m3]

        muts2 = annotator._applyManualAnnotations(muts, overrides)
        muts_final_gen = annotator._applyDefaultAnnotations(muts2, default_annotations)

        muts_final = []
        for m in muts_final_gen:
            self.assertTrue(m['test3'] == "foo3", "Override did not work")
            muts_final.append(m)

        self.assertTrue(muts_final[0]['test1'] == "foo1")
        self.assertTrue(muts_final[0]['test2'] == "foo2")
        self.assertTrue(muts_final[0]['test3'] == "foo3")

        self.assertTrue(muts_final[1]['test1'] == "")
        self.assertTrue(muts_final[1]['test2'] == "foo2")
        self.assertTrue(muts_final[1]['test3'] == "foo3")

        self.assertTrue(muts_final[2]['test1'] == "")
        self.assertTrue(muts_final[2]['test2'] == "foo2-original")
        self.assertTrue(muts_final[2]['test3'] == "foo3")

    def testAnnotateListOfMutations(self):
        """Test that we can initialize an Annotator, without an input or output and then feed mutations,
        one at a time... using a runspec"""

        # Locate the datasource directory and create a runspec
        dbDir = self.config.get("DEFAULT", "dbDir")
        ds = DatasourceFactory.createDatasources(dbDir)
        runSpec = RunSpecification()
        runSpec.initialize(None, None, datasources=ds)

        # Initialize the annotator with the runspec
        annotator = Annotator()
        annotator.initialize(runSpec)

        m = MutationData()
        m.chr = "1"
        m.start = "12941796"
        m.end = "12941796"
        m.alt_allele = "G"
        m.ref_allele = "T"

        muts = [m]

        muts = annotator.annotate_mutations(muts)
        m2 = muts.next()
        self.assertTrue(m2.get("gene", None) is not None)

    def testSkippingAltsForSingleMut(self):
        """Test a simple case where a single mutation with alt_allele_seen of False is not produced."""

        runSpec = RunSpecification()
        runSpec.initialize(None, None, datasources=[], is_skip_no_alts=True)

        # Initialize the annotator with the runspec
        annotator = Annotator()
        annotator.initialize(runSpec)

        m = MutationData()
        m.chr = "1"
        m.start = "12941796"
        m.end = "12941796"
        m.alt_allele = "G"
        m.ref_allele = "T"
        m.createAnnotation("alt_allele_seen", "False")

        muts = [m]

        muts = annotator.annotate_mutations(muts)
        self.assertRaises(StopIteration, muts.next)

    def _simple_annotate(self, is_skip_no_alts):
        runSpec = RunSpecification()
        runSpec.initialize(None, None, datasources=[], is_skip_no_alts=is_skip_no_alts)
        # Initialize the annotator with the runspec
        annotator = Annotator()
        annotator.initialize(runSpec)
        m = MutationData()
        m.chr = "1"
        m.start = "12941796"
        m.end = "12941796"
        m.alt_allele = "G"
        m.ref_allele = "T"
        m.createAnnotation("alt_allele_seen", "False")
        m2 = MutationData()
        m2.chr = "1"
        m2.start = "12941796"
        m2.end = "12941796"
        m2.alt_allele = "G"
        m2.ref_allele = "T"
        muts = [m, m2]
        muts = annotator.annotate_mutations(muts)
        ctr = 0
        for m in muts:
            ctr += 1
        return ctr

    def testSkippingAlts(self):
        """Test a simple case where mutations with alt_allele_seen of False are not produced."""

        ctr = self._simple_annotate(True)
        self.assertTrue(ctr == 1)

    def testSkippingAltsFalse(self):
        """Test a simple case that is_skip_alts of False does not affect anything."""

        ctr = self._simple_annotate(False)
        self.assertTrue(ctr == 2)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testBasicAnnotatorInit']
    unittest.main()