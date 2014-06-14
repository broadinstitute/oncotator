# LICENSE_GOES_HERE
from oncotator.Metadata import Metadata
from TestUtils import TestUtils


'''
Created on Feb 14, 2013

@author: lichtens
'''
import unittest
from oncotator.MutationData import MutationData
from oncotator.output.SimpleBedOutputRenderer import SimpleBedOutputRenderer

TestUtils.setupLogging(__file__, __name__)
class SimpleBedOutputRendererTest(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testSimpleRendering(self):
        m = MutationData()
        m.chr = '1'
        m.start = 1000000
        m.end = 1000000
        outputFilename = "out/simpleBEDTest.bed"
        outputRenderer = SimpleBedOutputRenderer(outputFilename)

        outputRenderer.renderMutations([m], Metadata())
        
        fp = file(outputFilename,'r')
        mOut = fp.readline().strip().split(' ')
        self.assertTrue(mOut[0] == "chr1")
        self.assertTrue(mOut[1] == "999999")
        self.assertTrue(mOut[2] == "1000000")
        fp.close()
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSimpleRendering']
    unittest.main()