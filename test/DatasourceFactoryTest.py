# LICENSE_GOES_HERE
from TestUtils import TestUtils


"""
Created on Jan 4, 2013

@author: lichtens
"""
import os
import unittest
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MutationData import MutationData
import logging
from oncotator.datasources.Datasource import Datasource

TestUtils.setupLogging(__file__, __name__)
class DatasourceFactoryTest(unittest.TestCase):

    _multiprocess_can_split_ = True

    # HACK: Allow config to be viewed by unittest decorators.
    globalConfig = TestUtils.createUnitTestConfig()

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()

    def tearDown(self):
        pass


    def testBasicCosmicInit(self):
        """ Very simple test that will create a datasource from a sample datasource directory.  
        The directory conforms to the standard datasource structure, including placement of the config file.
        """
        ds = DatasourceFactory.createDatasource('testdata/small_cosmic/small_cosmic.config', "testdata/small_cosmic")
        
        m = MutationData()
        m.chr = 19
        m.start = 58858921
        m.end = 58858921
        
        m = ds.annotate_mutation(m)
        
        self.assertTrue(m['COSMIC_overlapping_mutation_AAs'] == 'p.P426P(1)', "Did not properly annotate mutation: " + m['COSMIC_overlapping_mutation_AAs'])

    def testBasicRefInit(self):
        """ Very simple test that will create a reference datasource from a sample datasource directory.  
        The directory conforms to the standard datasource structure, including placement of the config file.
        """
        ds = DatasourceFactory.createDatasource('testdata/reference_ds/reference_ds.config', "testdata/reference_ds")
        
        m = MutationData()
        m.chr = "22"
        m.start = "11"
        m.end = "11"
        
        groundTruth = "CCCAAGCTAAACCCAGGCCAC"
        
        # remember that the annotate_mutation returns a generator, so we use an iterator
        m = ds.annotate_mutation(m)
        
        self.assertTrue(m['ref_context'] == groundTruth, "ref_context was not populated properly: " + str(m['ref_context']))

    def testBasicGeneTSVInit(self):
        """ Make sure that we can initialize a simple tsv data source """

        geneDS = DatasourceFactory.createDatasource("testdata/small_tsv_ds/small_tsv_ds.config", "testdata/small_tsv_ds/")
        self.assertTrue(geneDS <> None, "gene indexed datasource was None.")
        
        m = MutationData()
        m.createAnnotation('gene',"ABL1")
        m = geneDS.annotate_mutation(m)
        self.assertTrue(m['CGC_Abridged_Name'] == "v-abl Abelson murine leukemia viral oncogene homolog 1","Test gene TSV datasource did not annotate properly.")
    
    def testBasicDatasourceSorting(self):
        """Test that the GAF datasource is sorted before a gene-based datasource"""

        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
        geneDS = DatasourceFactory.createDatasource("testdata/small_tsv_ds/small_tsv_ds.config", "testdata/small_tsv_ds/")
        
        incorrectSortList = [geneDS, gafDatasource]
        guessSortList =  DatasourceFactory.sortDatasources(incorrectSortList)
        self.assertTrue(guessSortList[1] == geneDS, "Sorting is incorrect.")
        self.assertTrue(len(guessSortList) == 2, "Sorting altered number of datasources (gt: 2): " + str(len(guessSortList)))

    @unittest.skipIf(not os.path.exists(globalConfig.get("DEFAULT", "dbDir")), "Default Datasource corpus is needed to run this test")
    def testInitializingDatasources(self):
        """ Test initializing a database dir, both single and multicore.  This test is RAM intensive and requires default data corpus."""
        
        multiDS = DatasourceFactory.createDatasources(self.config.get("DEFAULT", "dbDir"), "hg19", isMulticore=True)
        self.assertTrue(multiDS is not None, "Datasource list was None")
        self.assertTrue(len(multiDS) != 0, "Datasource list was empty")
        for i in range(0,len(multiDS)):
            self.assertTrue(multiDS[i] is not None, "multi core datasource was None:  " + str(i))
            self.assertTrue(isinstance(multiDS[i],Datasource))

        # This test can be memory intensive, so get rid of the multiDS, but record how many datasources were created.
        numMultiDS = len(multiDS)
        del multiDS

        singleCoreDS = DatasourceFactory.createDatasources(self.config.get("DEFAULT", "dbDir"), "hg19", isMulticore=False)
        self.assertTrue(singleCoreDS is not None, "Datasource list was None")
        self.assertTrue(len(singleCoreDS) != 0, "Datasource list was empty")
        for i in range(0,len(singleCoreDS)):
            self.assertTrue(singleCoreDS[i] is not None, "single core datasource was None:  " + str(i))
            self.assertTrue(isinstance(singleCoreDS[i],Datasource))
            
        self.assertTrue(numMultiDS == len(singleCoreDS), "Length of single core datasource list was not the same as multicore")
        del singleCoreDS

    def testMulticoreExceptionCatching(self):
        """ Test that a datasource throws an exception during initialization, the DatasourceCreator does not freeze. """
        datasourceTuples = [("testdata/mock_exception_throwing_ds/mock_exception_throwing_ds.config",
                             "testdata/mock_exception_throwing_ds/")]


        #DatasourceCreator._createDatasourcesMulticore(4, datasourceTuples)
        self.assertRaisesRegexp(NotImplementedError,"This class throws exception",
                                DatasourceFactory._createDatasourcesMulticore,4, datasourceTuples)

    def testMulticoreNoDatasources(self):
        """ If using multicore, does not hang when no datasources are in the db dir"""
        multiDS = DatasourceFactory.createDatasources('testdata/maflite/', "hg19", True)
        self.assertTrue(len(multiDS) == 0, "Length of multiDS when there were no datasources was not zero.")

    def test_hashcode_generation(self):
        """Test that we can read a hashcode for a datasource, if available."""
        geneDS = DatasourceFactory.createDatasource("testdata/thaga_janakari_gene_ds/hg19/tj_data.config", "testdata/thaga_janakari_gene_ds/hg19/")
        self.assertTrue(geneDS is not None, "gene indexed datasource was None.")

        self.assertTrue(geneDS.get_hashcode() == "7120edfdc7b29e45191c81c99894afd5")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testBasicInit']
    unittest.main()