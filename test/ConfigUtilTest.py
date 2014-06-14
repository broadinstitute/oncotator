# LICENSE_GOES_HERE


'''
Created on Jan 16, 2013

@author: lichtens
'''
import unittest
from oncotator.utils.ConfigUtils import ConfigUtils
from ConfigParser import SafeConfigParser
from TestUtils import TestUtils
from TestUtils import TestUtils
TestUtils.setupLogging(__file__, __name__)
class ConfigUtilTest(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass

    def testLocateWithAdditionalDir(self):
        ''' Call the locate command on a config file and make sure the location is proper. '''
        fp = ConfigUtils._locateConfigFile("dummy.config", isRelaxedLogging=True, additional_config_dir="testdata/dummy_configs/")
        config = SafeConfigParser()
        config.readfp(fp)
        self.assertTrue(config.get("general", "dummy2") == "world")
        self.assertTrue(config.get("general", "dummy1") == "Hello")

    def testLocate(self):
        ''' Call the locate command on a config file and make sure the location is proper. '''
        fp = ConfigUtils._locateConfigFile("testdata/dummy_configs/dummy.config", isRelaxedLogging=True)
        config = SafeConfigParser()
        config.readfp(fp)
        self.assertTrue(config.get("general", "dummy1") == "Hello")
        self.assertTrue(config.get("general", "dummy2") == "world")
    
    def testLocalConfig(self):
        ''' Get a key from a local config and a config file that are the same basic name, but different values. '''
        config = ConfigUtils.createConfigParser("testdata/dummy_configs/dummy.config")
        self.assertTrue(config.get("general", "dummy1") == "Super")
        self.assertTrue(config.get("general", "dummy2") == "world")
     


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testLocate']
    unittest.main()