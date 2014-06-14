# LICENSE_GOES_HERE
from TestUtils import TestUtils


'''
Created on Jan 7, 2013

@author: lichtens
'''
from math import fabs
import unittest
import logging
from oncotator.datasources.ReferenceDatasource import ReferenceDatasource
from oncotator.MutationData import MutationData

TestUtils.setupLogging(__file__, __name__)
class ReferenceDatasourceTest(unittest.TestCase):


    def setUp(self):
        self.logger = logging.getLogger(__name__)


    def tearDown(self):
        pass


    def testSimpleGLAnnotate(self):
        ''' Test a simple annotation case.  Make sure that the ref_context and gc_content annotations are correct. '''
        ds = ReferenceDatasource('testdata/reference_ds', windowSizeGCContent=5)
        m = MutationData()
        m.chr = "GL000211.1"
        m.start = "11"
        m.end = "11"
        
        groundTruth = "gaattctttttcaagtaagtc"
        
        guess = ds.annotate_mutation(m)
        
        self.assertTrue(guess['ref_context'] == groundTruth, "ref_context was not populated properly: " + str(m['ref_context']))

        # gc_content is rounded to 3 decimal places
        self.assertTrue(fabs(float(guess['gc_content']) - (float(3)/float(11))) < .001, "gc_content was not populated properly: " + str(m['gc_content']))

    def testExtentOutOfRangeError(self):
        ''' If a window is specified that extends beyond the beginning or end of a file, truncate the ref_context.  
        Use what is left for gc_content as well.'''
        ds = ReferenceDatasource('testdata/reference_ds', windowSizeRef=6, windowSizeGCContent=5)
        m = MutationData()
        m.chr = "22"
        m.start = "4"
        m.end = "4"
        
        # "CCCAAGCTAAACCCAGGCCAC"
        groundTruth = "CCCAAGCTAA"
        
        guess = ds.annotate_mutation(m)
        
        self.assertTrue(guess['ref_context'] == groundTruth, "ref_context was not populated properly: " + str(guess['ref_context']))

        # gc_content is rounded to 3 decimal places
        self.assertTrue(fabs(float(guess['gc_content']) - (float(5)/float(9))) < .001, "gc_content was not populated properly: " + str(guess['gc_content']))

    def testSimpleAnnotate(self):
        ''' Perform a simple test of one of the aligned chromosomes (chr22) and make sure that we get a reasonable answer.
        '''
        ds = ReferenceDatasource('testdata/reference_ds', windowSizeGCContent=5)
        m = MutationData()
        m.chr = "22"
        m.start = "11"
        m.end = "11"
        
        groundTruth = "CCCAAGCTAAACCCAGGCCAC"

        guess = ds.annotate_mutation(m)
        
        self.assertTrue(guess['ref_context'] == groundTruth, "ref_context was not populated properly: " + str(guess['ref_context']))

        # gc_content is rounded to 3 decimal places
        self.assertTrue(fabs(float(guess['gc_content'])- (float(6)/float(11))) < .001, "gc_content was not populated properly: " + str(guess['gc_content']))
    
    def testFilenameDetermination(self):
        ''' Test that proper conversions are being done for chromosome to flat filename '''
        ds = ReferenceDatasource('testdata/reference_ds')
        self.assertTrue(ds.convertMutationChrToFilename("GL000211.1") == 'chrUn_gl000211.txt', "Did not find GL file: " + str(ds.convertMutationChrToFilename("GL000211.1"))) 
        self.assertTrue(ds.convertMutationChrToFilename("X") == 'chrX.txt', "Did not find chrX file: " + str(ds.convertMutationChrToFilename("X")))
        self.assertTrue(ds.convertMutationChrToFilename("GL000209.1") == 'chr19_gl000209_random.txt', "Did not find GL chr19 file: " + str(ds.convertMutationChrToFilename("GL000209.1")))
    
    def testEmptyAnswer(self):
        ''' The Reference Datasource should return a blank result if the chromosome is not found.
        Note: A log entry should also be written, but this is not tested. '''
        self.logger.info("Please ignore the next logging warning: testdata/reference_ds/chrTHIS_DOES_NOT_EXIST.txt not found.  Please add it.")
        ds = ReferenceDatasource('testdata/reference_ds')
        m = MutationData()
        m.chr = "THIS_DOES_NOT_EXIST"
        m.start = "11"
        m.end = "11"
        
        groundTruth = ""
        # remember that the annotate_mutation returns a generator, so we use an iterator
        guess = ds.annotate_mutation(m)
        self.assertTrue(guess['ref_context'] == groundTruth, "ref_context was not populated properly -- should be blank: " + str(guess['ref_context']))
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()