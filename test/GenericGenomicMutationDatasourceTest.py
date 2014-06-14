# LICENSE_GOES_HERE
from TestUtils import TestUtils


'''
Created on Jan 16, 2013

@author: aramos
'''
import unittest

import os

from oncotator.datasources.GenericGenomicMutationDatasource import GenericGenomicMutationDatasource
from oncotator.MutationData import MutationData
import logging

TestUtils.setupLogging(__file__, __name__)
class GenericGenomicMutationDatasourceTest(unittest.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()


    def tearDown(self):
        pass

    @unittest.skip("Currently fails and code has not been debugged.  Also, this datasource is not used yet.")
    def testBasicAnnotation(self):
        ds = GenericGenomicMutationDatasource('testdata/small_cosmic_2/cosmic_v65_chr18.tsv')
    
        m = MutationData()
        m.chr = '18'
        m.start = '48604683'
        m.end = '48604683'
        m.ref_allele = 'G'
        m.alt_allele = 'A'
        m.createAnnotation('strand', '+')
    
        guess = ds.annotate_mutation(m)
        self.assertTrue(guess['_cosmic_muts_disease_counts'], 'Unable to annotate mutation correctly')
            
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()