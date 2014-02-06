from oncotator.index.TabixIndexedVcfDatasourceCreator import TabixIndexedVcfDatasourceCreator
from test.TestUtils import TestUtils
from oncotator.utils.ConfigUtils import ConfigUtils
import os
import string
import vcf

__author__ = 'lichtens'

import unittest

TestUtils.setupLogging(__file__, __name__)


class IndexedVcfDatasourceBuilderTest(unittest.TestCase):

    def testCreateDatabase(self):
        dsFile = os.path.join(*["testdata", "vcf", "example.vcf"])
        destDir = "out"
        datasourceBuilder = TabixIndexedVcfDatasourceCreator()
        datasourceFilename = datasourceBuilder._createDatabase(destDir=destDir, ds_file=dsFile)
        tabixIndexedFilename = string.join([destDir, os.sep, datasourceFilename], "")

        self.assertTrue(os.path.exists(tabixIndexedFilename), "No index file was generated.")

        vcfReader = vcf.Reader(filename=tabixIndexedFilename, compressed=True, strict_whitespace=True)
        vcfRecords = vcfReader.fetch(chrom=20, start=1230237, end=1230237)
        for vcfRecord in vcfRecords:
            self.assertEqual(vcfRecord.INFO["NS"], 3, "Expected %s but got %s." % (3, vcfRecord.INFO["NS"]))
            self.assertEqual(vcfRecord.INFO["DP"], 13, "Expected %s but got %s." % (13, vcfRecord.INFO["DP"]))

    def testCreateConfigFile(self):
        configFilename = os.path.join("out", "esp.config")
        datasourceFilename = "ESP6500SI-V2.vcf.gz"
        dataSourceType = "indexed_vcf"
        datasourceMatchMode = "avg"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"

        datasourceBuilder = TabixIndexedVcfDatasourceCreator()
        datasourceBuilder._createConfigFile(configFilename=configFilename, baseDSFile=datasourceFilename,
                                            ds_type=dataSourceType, ds_name=dataSourceName,
                                            ds_version=dataSourceVersion, ds_match_mode=datasourceMatchMode)
        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "match_mode"),
                        "match_mode option is missing in general section.")

        self.assertEqual(configParser.get("general", "type"), dataSourceType,
                         "Expected data source type is %s but was %s."
                         % (dataSourceType, configParser.get("general", "type")))
        self.assertEqual(configParser.get("general", "src_file"), datasourceFilename,
                         "Expected data source src_file is %s but was %s."
                         % (datasourceFilename, configParser.get("general", "src_file")))
        self.assertEqual(configParser.get("general", "title"), dataSourceName,
                         "Expected data source title is %s but was %s."
                         % (dataSourceName, configParser.get("general", "title")))
        self.assertEqual(configParser.get("general", "version"), dataSourceVersion,
                         "Expected data source version is %s but was %s."
                         % (dataSourceVersion, configParser.get("general", "version")))
        self.assertEqual(configParser.get("general", "match_mode"), datasourceMatchMode,
                         "Expected data source match mode is %s but was %s."
                         % (datasourceMatchMode, configParser.get("general", "match_mode")))

    def testCreateDatasource(self):
        dsFile = os.path.join(*["testdata", "vcf", "example.vcf"])
        destDir = "out"
        configFilename = os.path.join("out", "esp.config")
        datasourceFilename = "example.tabix_indexed.vcf.gz"
        dataSourceType = "indexed_vcf"
        datasourceMatchMode = "avg"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"

        datasourceBuilder = TabixIndexedVcfDatasourceCreator()
        datasourceBuilder.createDatasource(destDir, dsFile, configFilename, dataSourceType, dataSourceName,
                                           dataSourceVersion, datasourceMatchMode)

        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "match_mode"),
                        "match_mode option is missing in general section.")

        self.assertEqual(configParser.get("general", "type"), dataSourceType,
                         "Expected data source type is %s but was %s."
                         % (dataSourceType, configParser.get("general", "type")))
        self.assertEqual(configParser.get("general", "src_file"), datasourceFilename,
                         "Expected data source src_file is %s but was %s."
                         % (datasourceFilename, configParser.get("general", "src_file")))
        self.assertEqual(configParser.get("general", "title"), dataSourceName,
                         "Expected data source title is %s but was %s."
                         % (dataSourceName, configParser.get("general", "title")))
        self.assertEqual(configParser.get("general", "version"), dataSourceVersion,
                         "Expected data source version is %s but was %s."
                         % (dataSourceVersion, configParser.get("general", "version")))
        self.assertEqual(configParser.get("general", "match_mode"), datasourceMatchMode,
                         "Expected data source match mode is %s but was %s."
                         % (datasourceMatchMode, configParser.get("general", "match_mode")))

    def testCreateDatasourceFromGZFile(self):
        dsFile = os.path.join("testdata", "example.vcf.gz")
        destDir = "out"
        configFilename = os.path.join("out", "esp.config")
        datasourceFilename = "example.vcf.gz"
        dataSourceType = "indexed_vcf"
        datasourceMatchMode = "avg"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"

        datasourceBuilder = TabixIndexedVcfDatasourceCreator()
        datasourceBuilder.createDatasource(destDir, dsFile, configFilename, dataSourceType, dataSourceName,
                                           dataSourceVersion, datasourceMatchMode)

        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "match_mode"),
                        "match_mode option is missing in general section.")

        self.assertEqual(configParser.get("general", "type"), dataSourceType,
                         "Expected data source type is %s but was %s."
                         % (dataSourceType, configParser.get("general", "type")))
        self.assertEqual(configParser.get("general", "src_file"), datasourceFilename,
                         "Expected data source src_file is %s but was %s."
                         % (datasourceFilename, configParser.get("general", "src_file")))
        self.assertEqual(configParser.get("general", "title"), dataSourceName,
                         "Expected data source title is %s but was %s."
                         % (dataSourceName, configParser.get("general", "title")))
        self.assertEqual(configParser.get("general", "version"), dataSourceVersion,
                         "Expected data source version is %s but was %s."
                         % (dataSourceVersion, configParser.get("general", "version")))
        self.assertEqual(configParser.get("general", "match_mode"), datasourceMatchMode,
                         "Expected data source match mode is %s but was %s."
                         % (datasourceMatchMode, configParser.get("general", "match_mode")))
