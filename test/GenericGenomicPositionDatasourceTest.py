# LICENSE_GOES_HERE
from TestUtils import TestUtils


'''
Created on Jan 16, 2013

@author: lichtens
'''
import unittest

import os

from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.Annotator import Annotator
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator 
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
import logging
from oncotator.utils.GenericTsvReader import GenericTsvReader

TestUtils.setupLogging(__file__, __name__)
class GenericGenomicPositionDatasourceTest(unittest.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()


    def tearDown(self):
        pass


    def testBasicAnnotation(self):
        ''' Annotate from a basic tsv of Genomic positions.  This tests both single- and multiple-nucleotide variants.  The tsv is already installed (i.e. proper config file created).
        '''
        outputFilename = 'out/genericGenomePositionTest.out.tsv'
        
        gpDS = DatasourceFactory.createDatasource("testdata/small_genome_position_tsv_ds/oreganno_trim.config", "testdata/small_genome_position_tsv_ds/")
        
        annotator = Annotator()
        annotator.setInputCreator(MafliteInputMutationCreator('testdata/maflite/tiny_maflite.maf.txt'))
        annotator.setOutputRenderer(SimpleOutputRenderer(outputFilename))
        annotator.addDatasource(gpDS)
        testFilename = annotator.annotate()
        
        # Make sure that some values were populated
        self.assertTrue(os.path.exists(testFilename))
        tsvReader = GenericTsvReader(testFilename) 
        
        ctr = 1
        # Two overlap, one does not.  Repeat...
        for lineDict in tsvReader:
            if (ctr % 3 == 0):
                self.assertTrue(lineDict["ORegAnno_hg19.oreganno.id"] == '', "Line " + str(ctr) + " should have had blank value, but did not: " + lineDict["ORegAnno_hg19.oreganno.id"])
            else:
                self.assertFalse(lineDict["ORegAnno_hg19.oreganno.id"] == '', "Line " + str(ctr) + " should not have had blank value, but did.")
                self.assertTrue(lineDict["ORegAnno_hg19.oreganno.id"] == 'OREG0013034', "Line " + str(ctr) + " did not have correct value: " + lineDict["ORegAnno_hg19.oreganno.id"])
            ctr = ctr + 1
    
    def testDoubleAnnotationError(self):
        ''' Given a maf file that used to cause a duplicate annotation exception, do not throw that (or any) exception. '''
        outputFilename = 'out/genericGenomePositionDoubleAnnotationTest.out.tsv'
        
        gpDS = DatasourceFactory.createDatasource("testdata/small_genome_position_tsv_ds/oreganno_trim.config", "testdata/small_genome_position_tsv_ds/")
        
        
        annotator = Annotator()
        annotator.setInputCreator(MafliteInputMutationCreator('testdata/maflite/testDoubleAnnotate.maf.tsv'))
        annotator.setOutputRenderer(SimpleOutputRenderer(outputFilename))
        annotator.addDatasource(gpDS)
        testFilename = annotator.annotate()
        
        
        # Make sure that some values were populated
        self.assertTrue(os.path.exists(testFilename))
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()