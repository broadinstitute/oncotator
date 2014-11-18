import unittest
import logging

from oncotator.MutationData import MutationData
import oncotator.datasources.BigWigDatasource
from TestUtils import TestUtils

TestUtils.setupLogging(__file__, __name__)

class BigWigDatasourceTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.bigwig_datasource = oncotator.datasources.BigWigDatasource.BigWigDatasource('testdata/bigwig/small.bigWig', title='TestBigWig')

    def tearDown(self):
        self.bigwig_datasource.close()

    def test_basic_fetch(self):
        m = MutationData()
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 78978)
        m.createAnnotation('end', 78978)
        
        self.bigwig_datasource.annotate_mutation(m)
        self.assertEqual(m.get('TestBigWig_score'), 0.5)

    def test_range_fetch(self):
        m = MutationData()
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 78978)
        m.createAnnotation('end', 79000)
        
        self.bigwig_datasource.annotate_mutation(m)
        self.assertEqual(m.get('TestBigWig_score'), 0.75)

    def test_no_data_fetch(self):
        """Test for value not found in bigwig.  In this case, our test bigwig only has data for 
        chr1 so None is expected return value.
        """
        m = MutationData()
        m.createAnnotation('chr', '13')
        m.createAnnotation('start', 78978)
        m.createAnnotation('end', 79000)
        
        self.bigwig_datasource.annotate_mutation(m)
        self.assertEqual(m.get('TestBigWig_score'), None)

if __name__ == '__main__':
    unittest.main()