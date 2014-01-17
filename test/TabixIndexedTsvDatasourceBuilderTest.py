import os
from oncotator.index.TabixIndexedTsvDatasourceCreator import TabixIndexedTsvDatasourceCreator
from test.TestUtils import TestUtils
from oncotator.utils.ConfigUtils import ConfigUtils
import pysam
import string
from oncotator.utils.install.DatasourceInstallUtils import DatasourceInstallUtils
from oncotator.index.InputMismatchException import InputMismatchException

__author__ = 'lichtens'

import unittest

TestUtils.setupLogging(__file__, __name__)


class IndexedTsvDatasourceBuilderTest(unittest.TestCase):

    def testCreateDatabase(self):
        dsFile = "testdata/ESP6500SI-V2.chr1.snps_indels.head.25.txt"
        destDir = "out"
        indexColumnNames = "CHROM,POS,POS"
        columnNames = "CHROM,POS,REF,ALT,DBSNP,EA_AC,AA_AC,TAC,MAF,GTS,EA_GTC,AA_GTC,GTC,DP,FG,GM,AA,AAC,PP,CDP,PH,CP,CG,GL,GS,CA,EXOME_CHIP,GWAS_PUBMED"
        annotationColumnNames = "CHROM,POS,REF,ALT,DBSNP"

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        datasourceFilename = datasourceBuilder._createDatabase(destDir=destDir, ds_file=dsFile,
                                                               index_column_names=indexColumnNames,
                                                               column_names=columnNames,
                                                               annotation_column_names=annotationColumnNames)
        tabixIndexedFilename = string.join([destDir, os.sep, datasourceFilename], "")

        self.assertTrue(os.path.exists(tabixIndexedFilename), "No index file was generated.")

        chrom = "1"
        start = "69594"
        end = "69594"
        tsvRecords = None
        tsvReader = pysam.Tabixfile(filename=tabixIndexedFilename)  # initialize the tsv reader
        try:
            tsvRecords = tsvReader.fetch(chrom, int(start)-1, int(end), parser=pysam.asTuple())
        except ValueError:
            pass

        tsvRecord = None
        for tsvRecord in tsvRecords:
            self.assertEqual(tsvRecord[5], "2,6190", "Value in column sixth does not match the expected value.")

        self.assertIsNotNone(tsvRecord, "No record for %s:%s-%s was found." % (chrom, start, end))

    def testCreateConfigFile(self):
        """


        """
        configFilename = "out/esp_coverage.config"
        datasourceFilename = "ESP6500SI-V2.coverage.txt.gz"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        indexColumnNames = "Chromosome,Position,Position"
        columnNames = "Chromosome,Position,TotalSamplesCovered,AvgSampleReadDepth,TotalEAsamplesCovered,AvgEAsampleReadDepth,TotalAAsamplesCovered,AvgAAsampleReadDepth"
        annotationColumnNames = "TotalSamplesCovered,AvgSampleReadDepth,AvgEAsampleReadDepth,TotalAAsamplesCovered,AvgAAsampleReadDepth"

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        datasourceBuilder._createConfigFile(configFilename=configFilename, baseDSFile=datasourceFilename,
                                            ds_type=dataSourceType, ds_name=dataSourceName,
                                            ds_version=dataSourceVersion, column_names=columnNames,
                                            annotation_column_names=annotationColumnNames,
                                            ds_match_mode=dataSourceMatchMode,
                                            indexCols=DatasourceInstallUtils.getIndexCols(dataSourceType,
                                                                                          indexColumnNames))
        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_section("data_types"), "data_types section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "column_names"),
                        "column_names option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "annotation_column_names"),
                        "annotation_column_names option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "match_mode"),
                        "match_mode option is missing in general section")
        self.assertTrue(configParser.has_option("general", "index_column_names"),
                        "index_column_names option is missing in general section.")

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
        self.assertEqual(configParser.get("general", "column_names"), columnNames,
                         "Expected data source column names is %s but was %s."
                         % (columnNames, configParser.get("general", "column_names")))
        self.assertEqual(configParser.get("general", "annotation_column_names"), annotationColumnNames,
                         "Expected data source annotation column names is %s but was %s."
                         % (annotationColumnNames, configParser.get("general", "annotation_column_names")))
        self.assertEqual(configParser.get("general", "match_mode"), dataSourceMatchMode,
                         "Expected data source match mode is %s but was %s."
                         % (dataSourceMatchMode, configParser.get("general", "match_mode")))
        self.assertEqual(configParser.get("general", "index_column_names"), indexColumnNames,
                         "Expected data source index column names is %s but was %s."
                         % (indexColumnNames, configParser.get("general", "index_column_names")))

    def testCreateDatasource(self):
        dsFile = "testdata/ESP6500SI-V2.chr1.snps_indels.head.25.txt"
        destDir = "out"
        datasourceFilename = "ESP6500SI-V2.chr1.snps_indels.head.25.tabix_indexed.txt.gz"
        indexColumnNames = "CHROM,POS,POS"
        columnNames = "CHROM,POS,REF,ALT,DBSNP,EA_AC,AA_AC,TAC,MAF,GTS,EA_GTC,AA_GTC,GTC,DP,FG,GM,AA,AAC,PP,CDP,PH,CP,CG,GL,GS,CA,EXOME_CHIP,GWAS_PUBMED"
        configFilename = "out/esp_coverage.config"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        annotationColumnNames = "DBSNP,EA_GTC,DP"

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        datasourceBuilder.createDatasource(destDir, dsFile, indexColumnNames, columnNames, configFilename,
                                           dataSourceType, dataSourceName, dataSourceVersion, dataSourceMatchMode,
                                           annotationColumnNames,
                                           DatasourceInstallUtils.getIndexCols(dataSourceType, indexColumnNames))

        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_section("data_types"), "data_types section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "column_names"),
                        "column_names option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "annotation_column_names"),
                        "annotation_column_names option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "match_mode"),
                        "match_mode option is missing in general section")
        self.assertTrue(configParser.has_option("general", "index_column_names"),
                        "index_column_names option is missing in general section.")

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
        self.assertEqual(configParser.get("general", "column_names"), columnNames,
                         "Expected data source column names is %s but was %s."
                         % (columnNames, configParser.get("general", "column_names")))
        self.assertEqual(configParser.get("general", "annotation_column_names"), annotationColumnNames,
                         "Expected data source annotation column names is %s but was %s."
                         % (annotationColumnNames, configParser.get("general", "annotation_column_names")))
        self.assertEqual(configParser.get("general", "match_mode"), dataSourceMatchMode,
                         "Expected data source match mode is %s but was %s."
                         % (dataSourceMatchMode, configParser.get("general", "match_mode")))
        self.assertEqual(configParser.get("general", "index_column_names"), indexColumnNames,
                         "Expected data source index column names is %s but was %s."
                         % (indexColumnNames, configParser.get("general", "index_column_names")))

        self.assertEqual(configParser.get("data_types", "EA_GTC"), "String",
                         "Expected EA_GTC data type is %s but was %s."
                         % ("String", configParser.get("data_types", "EA_GTC")))
        self.assertEqual(configParser.get("data_types", "DP"), "Integer",
                         "Expected DP data type is %s but was %s."
                         % ("Integer", configParser.get("data_types", "DP")))

    def testCreateDatasourceWithMissingValues(self):
        dsFile = "testdata/ESP6500SI-V2.chr1.snps_indels.head.25.missing.txt"
        destDir = "out"
        datasourceFilename = "ESP6500SI-V2.chr1.snps_indels.head.25.missing.tabix_indexed.txt.gz"
        indexColumnNames = "CHROM,POS,POS"
        columnNames = "CHROM,POS,REF,ALT,EA_GTC,DP"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        annotationColumnNames = "EA_GTC,DP"
        configFilename = "out/esp_coverage.missing.config"

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        datasourceBuilder.createDatasource(destDir, dsFile, indexColumnNames, columnNames, configFilename,
                                           dataSourceType, dataSourceName, dataSourceVersion, dataSourceMatchMode,
                                           annotationColumnNames,
                                           DatasourceInstallUtils.getIndexCols(dataSourceType, indexColumnNames))

        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertEqual(configParser.get("general", "src_file"), datasourceFilename,
                         "Expected data source src_file is %s but was %s."
                         % (datasourceFilename, configParser.get("general", "src_file")))

        self.assertEqual(configParser.get("data_types", "EA_GTC"), "Float",
                         "Expected EA_GTC data type is %s but was %s."
                         % ("Float", configParser.get("data_types", "EA_GTC")))
        self.assertEqual(configParser.get("data_types", "DP"), "Integer",
                         "Expected DP data type is %s but was %s."
                         % ("Integer", configParser.get("data_types", "DP")))

    def testCreateDatasourceWithMissingColumns(self):
        dsFile = "testdata/ESP6500SI-V2.chr1.snps_indels.head.25.txt"
        destDir = "out"
        indexColumnNames = "CHROM,POS,POS"
        columnNames = "CHROM,POS,REF,ALT,EA_GTC,DP,ESP_DBSNP"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        annotationColumnNames = "EA_GTC,DP"
        configFilename = "out/esp_coverage.missing.config"

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        try:
            datasourceBuilder.createDatasource(destDir, dsFile, indexColumnNames, columnNames, configFilename,
                                               dataSourceType, dataSourceName, dataSourceVersion, dataSourceMatchMode,
                                               annotationColumnNames,
                                               DatasourceInstallUtils.getIndexCols(dataSourceType, indexColumnNames))
        except InputMismatchException:
            pass

    def testCreateDatasourceWithMissingAnnotationColumns(self):
        dsFile = "testdata/ESP6500SI-V2.chr1.snps_indels.head.25.txt"
        destDir = "out"
        indexColumnNames = "CHROM,POS,POS"
        columnNames = "CHROM,POS,REF,ALT,DBSNP,EA_AC,AA_AC,TAC,MAF,GTS,EA_GTC,AA_GTC,GTC,DP,FG,GM,AA,AAC,PP,CDP,PH,CP,CG,GL,GS,CA,EXOME_CHIP,GWAS_PUBMED"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        annotationColumnNames = "EA_GTC,DP,ESP_DBSNP"
        configFilename = "out/esp_coverage.missing.config"

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        try:
            datasourceBuilder.createDatasource(destDir, dsFile, indexColumnNames, columnNames, configFilename,
                                               dataSourceType, dataSourceName, dataSourceVersion, dataSourceMatchMode,
                                               annotationColumnNames,
                                               DatasourceInstallUtils.getIndexCols(dataSourceType, indexColumnNames))
        except InputMismatchException as e:
            print e
            pass