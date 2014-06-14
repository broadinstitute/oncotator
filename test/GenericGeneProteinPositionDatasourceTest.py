# LICENSE_GOES_HERE


from oncotator.Annotator import Annotator
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MissingAnnotationException import MissingAnnotationException
from oncotator.MutationData import MutationData
from oncotator.datasources.GenericGeneProteinPositionDatasource import GenericGeneProteinPositionDatasource
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.utils.GenericTsvReader import GenericTsvReader
from TestUtils import TestUtils

__author__ = 'lichtens'
import logging
import unittest
import os

TestUtils.setupLogging(__file__, __name__)
class GenericGeneProteinPositionDatasourceTest(unittest.TestCase):
    _multiprocess_can_split_ = True
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()

    def tearDown(self):
        pass

    def testBasicAnnotation(self):
        ''' Test an extremely simple case.
        '''
        datasource = GenericGeneProteinPositionDatasource("testdata/simple_uniprot_natvar/simple_uniprot_natvar.tsv", title="UniProt_NatVar", version="2011_09")

        m = MutationData()
        m.createAnnotation("gene", "TP53")
        m.createAnnotation("protein_change", "p.S376C")
        m.createAnnotation("other_transcripts", "TP53_uc002gig.1_Intron|TP53_uc002gih.2_Intron|TP53_uc010cne.1_RNA|TP53_uc010cnf.1_3'UTR|TP53_uc010cng.1_3'UTR|TP53_uc002gii.1_Missense_Mutation_p.S244C|TP53_uc010cnh.1_3'UTR|TP53_uc010cni.1_3'UTR|TP53_uc002gij.2_Missense_Mutation_p.S376C")

        m2 = datasource.annotate_mutation(m)
        annotationName= "UniProt_NatVar_natural_variations"
        self.assertTrue(sorted(m[annotationName].split("|")) == sorted("S -> T (in a sporadic cancer; somatic mutation).|S -> A (in a sporadic cancer; somatic mutation).".split("|")), "Incorrect annotation value seen: " + m[annotationName])

    def testCreationAndAnnotation(self):
        """ Test the datasource creation and then do a simple annotation
        """
        outputFilename = 'out/genericGeneProteinPositionTest.out.tsv'

        gafDS = TestUtils.createTranscriptProviderDatasource(self.config)
        gppDS = DatasourceFactory.createDatasource("testdata/simple_uniprot_natvar/simple_uniprot_natvar.config", "testdata/simple_uniprot_natvar/")

        annotator = Annotator()
        annotator.setInputCreator(MafliteInputMutationCreator('testdata/maflite/tiny_maflite_natvar.maf.tsv'))
        annotator.setOutputRenderer(SimpleOutputRenderer(outputFilename))
        annotator.addDatasource(gafDS)
        annotator.addDatasource(gppDS)
        testFilename = annotator.annotate()

        # Make sure that some values were populated
        self.assertTrue(os.path.exists(testFilename))
        tsvReader = GenericTsvReader(testFilename)

        ctr = 0
        for lineDict in tsvReader:
            colName = "UniProt_NatVar_natural_variations"
            self.assertTrue(sorted(lineDict[colName].split("|")) == sorted("R -> RR (in EDMD2).|R -> Q (in EDMD2).".split("|")), "Annotation value did not match: " + lineDict[colName])
            ctr += 1

        self.assertTrue(ctr == 1, "Number of mutations incorrect (1): " + str(ctr) )

    def testMissingAnnotations(self):
        ''' Tests that if the required annotations ("gene", "protein_change", and "other_transcripts") are missing, an excpetion is thrown.
        '''
        datasource = GenericGeneProteinPositionDatasource("testdata/simple_uniprot_natvar/simple_uniprot_natvar.tsv", title="SmallNatVar", version="test")

        m = MutationData()
        m.createAnnotation("gene", "TP53")
        #m.createAnnotation("protein_change", "p.S376C")

        self.assertRaisesRegexp(MissingAnnotationException, "protein_change", datasource.annotate_mutation, m)

