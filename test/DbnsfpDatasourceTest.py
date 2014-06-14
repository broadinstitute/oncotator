# LICENSE_GOES_HERE

import unittest
import logging

from oncotator.MutationData import MutationData
from oncotator.datasources.dbNSFP import dbNSFP
from TestUtils import TestUtils


TestUtils.setupLogging(__file__, __name__)
class DbnsfpDatasourceTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()

    def testBasicAnnotate(self):
        ds = dbNSFP('testdata/dbnsfp_ds')

        m = MutationData()
        m.chr = 'Y'
        m.start = '2655175'
        m.end = '2655175'
        m.ref_allele = 'A'
        m.alt_allele = 'T'
        m.createAnnotation('variant_classification', 'Missense_Mutation')
        m.createAnnotation('protein_change', 'p.V157E')

        guess = ds.annotate_mutation(m)
        mutation_assessor_prediction = guess['MutationAssessor_pred']

        self.assertTrue(mutation_assessor_prediction == 'medium', 'Unable to retrieve correct Mutation Assessor prediction')

if __name__ == '__main__':
    unittest.main()
