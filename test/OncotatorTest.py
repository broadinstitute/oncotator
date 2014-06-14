# LICENSE_GOES_HERE
from TestUtils import TestUtils


'''
Created on Oct 23, 2012

@author: gavlee
'''
import unittest

TestUtils.setupLogging(__file__, __name__)
class OncotatorTest(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testImports(self):
        ''' A simple test of import statements to make sure environment is set up correctly. '''
        missingImports = []
        try:
            import vcf
        except:
            missingImports.append("vcf")
            
        try:
            import pysam
        except:
            missingImports.append("pysam")
        try:
            import pandas
        except:
            missingImports.append("pandas")
        try:
            import Bio
        except:
            missingImports.append("biopython")
        try:
            import shelve
        except:
            missingImports.append("shelve")
        try:
            import numpy
        except:
            missingImports.append("numpy")
        
        self.assertTrue(len(missingImports) == 0, "The following modules could not be imported: " + str(missingImports) + ".  Try 'pip install <packagename>'")
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testImports']
    unittest.main()