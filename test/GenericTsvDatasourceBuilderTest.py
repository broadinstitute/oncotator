import os
from oncotator.index.GenericTsvDatasourceCreator import GenericTsvDatasourceCreator
from test.TestUtils import TestUtils
from oncotator.utils.install.DatasourceInstallUtils import DatasourceInstallUtils
from oncotator.utils.ConfigUtils import ConfigUtils
import pysam
import string

__author__ = 'lichtens'

import unittest

TestUtils.setupLogging(__file__, __name__)


class GenericTsvDatasourceBuilderTest(unittest.TestCase):

    def testCreateGPTsvDatasource(self):
        dsFile = "testdata/small_genome_position_tsv_ds/oreganno_trim.hg19.txt"
        destDir = "out"
        datasourceBuilder = GenericTsvDatasourceCreator()
        datasourceFilename = datasourceBuilder.createDatasource(destDir=destDir, ds_file=dsFile)
        datasourceFilename = string.join([destDir, os.sep, datasourceFilename], "")

        self.assertTrue(os.path.exists(datasourceFilename), "No data source file was generated.")

    def testCreateGPTsvConfigFile(self):
        configFilename = "out/ccle_by_gp.config"
        datasourceFilename = "ccle_results_by_pos.hg19.import.txt"
        dataSourceType = "gp_tsv"
        dataSourceName = "CCLE_By_GP"
        dataSourceVersion = "09292010"
        genomicPositionColumnNames = "chr,start,end"

        datasourceBuilder = GenericTsvDatasourceCreator()
        datasourceBuilder.createConfigFile(configFilename=configFilename, baseDSFile=datasourceFilename,
                                           ds_name=dataSourceName, ds_type=dataSourceType, ds_version=dataSourceVersion,
                                           indexCols=DatasourceInstallUtils.getIndexCols("gp_tsv",
                                                                                         genomicPositionColumnNames))

        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "genomic_position_cols"),
                        "genomic_position_cols option is missing in general section.")

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
        self.assertEqual(configParser.get("general", "genomic_position_cols"), genomicPositionColumnNames,
                         "Expected data source genomic_position_cols is %s but was %s."
                         % (genomicPositionColumnNames, configParser.get("general", "genomic_position_cols")))

    def getGeneTsvConfigFile(self):
        configFilename = "out/simple_uniprot.config"
        datasourceFilename = "simple_uniprot.out.2011_09.tsv"
        dataSourceType = "gene_tsv"
        dataSourceName = "UniProt"
        dataSourceVersion = "2011_09"
        geneColumnName = "gene"

        datasourceBuilder = GenericTsvDatasourceCreator()
        datasourceBuilder.createConfigFile(configFilename=configFilename, baseDSFile=datasourceFilename,
                                           ds_name=dataSourceName, ds_type=dataSourceType, ds_version=dataSourceVersion,
                                           indexCols=DatasourceInstallUtils.getIndexCols("gene_tsv", geneColumnName))

        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "gene_col"),
                        "gene_col option is missing in general section.")

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
        self.assertEqual(configParser.get("general", "gene_col"), geneColumnName,
                         "Expected data source gene_col is %s but was %s."
                         % (geneColumnName, configParser.get("general", "gene_col")))
