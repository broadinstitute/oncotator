# LICENSE_GOES_HERE
from TestUtils import TestUtils


'''
Created on Jan 22, 2013

@author: lichtens
'''
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.Annotator import Annotator
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
import unittest
import logging
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.MutationData import MutationData

TestUtils.setupLogging(__file__, __name__)


class GenericTranscriptDatasourceTest(unittest.TestCase):

    _multiprocess_can_split_ = True
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()


    def tearDown(self):
        pass


    def testBasicAnnotation(self):
        ''' Test annotation from a generic TSV based on a transcript annotation.  Only confirms the proper headers of the output. '''
        # We need a gaf data source to annotate gene

        gafDatasource = TestUtils.createTranscriptProviderDatasource(config=self.config)
        transcriptDS = DatasourceFactory.createDatasource("testdata/small_transcript_tsv_ds/small_transcript_tsv_ds.config", "testdata/small_transcript_tsv_ds/")
        outputFilename = 'out/genericTranscriptTest.out.tsv'
        
        annotator = Annotator()
        annotator.setInputCreator(MafliteInputMutationCreator('testdata/maflite/Patient0.snp.maf.txt'))
        annotator.setOutputRenderer(SimpleOutputRenderer(outputFilename))
        annotator.addDatasource(gafDatasource)
        annotator.addDatasource(transcriptDS)
        outputFilename = annotator.annotate()

        tsvReader = GenericTsvReader(outputFilename) 
        headers = tsvReader.getFieldNames()        
        self.assertTrue("refseq_test_mRNA_Id" in headers, "refseq_test_mRNA_Id not found in headers: " + str(headers))
        self.assertTrue("refseq_test_prot_Id" in headers, "refseq_test_prot_Id not found in headers: " + str(headers))
        

    def testSimpleAnnotation(self):
        ''' Create a dummy mutation and make sure it gets annotated properly '''
        m = MutationData()
        m.createAnnotation('transcript_id', 'uc001hms.3')
        transcriptDS = DatasourceFactory.createDatasource("testdata/small_transcript_tsv_ds/small_transcript_tsv_ds.config", "testdata/small_transcript_tsv_ds/")
        m = transcriptDS.annotate_mutation(m)
        self.assertTrue(m['refseq_test_mRNA_Id'] == 'NM_022746', "Transcript-based annotation did not populate properly: " + m['refseq_test_mRNA_Id'])
        self.assertTrue(m['refseq_test_prot_Id'] == 'NP_073583', "Transcript-based annotation did not populate properly: " + m['refseq_test_prot_Id'])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()