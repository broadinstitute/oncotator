"""
By downloading the PROGRAM you agree to the following terms of use:

BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY

This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").

WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and

WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.

NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:

1. DEFINITIONS
1.1 "PROGRAM" shall mean the object code and source code known as Oncotator 1.0 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.  BROAD acknowledges that the PROGRAM employs one or more public domain code(s) that are freely available for public use.

2. LICENSE
2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.  LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.

3. OWNERSHIP OF INTELLECTUAL PROPERTY
LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.

Copyright 2014 Broad Institute, Inc.
Notice of attribution:  The Oncotator 1.0 program was made available through the generosity of the Broad Institute, Inc.

LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.

4. INDEMNIFICATION
LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorney fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.

5. NO REPRESENTATIONS OR WARRANTIES
THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.

6. ASSIGNMENT
This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.

7. MISCELLANEOUS
7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
"""
import os
from oncotator.index.TabixIndexedTsvDatasourceCreator import TabixIndexedTsvDatasourceCreator
from test.TestUtils import TestUtils
from oncotator.utils.ConfigUtils import ConfigUtils
import pysam
import string
from oncotator.utils.install.DatasourceInstallUtils import DatasourceInstallUtils
from oncotator.index.InputMismatchException import InputMismatchException
import unittest

TestUtils.setupLogging(__file__, __name__)


class TabixIndexedTsvDatasourceBuilderTest(unittest.TestCase):

    def testCreateDatabase(self):
        """

        """
        dsFile = os.path.join("testdata", "ESP6500SI-V2.chr1.snps_indels.head.25.txt")
        destDir = "out"
        indexColumnNames = ["CHROM", "POS", "POS"]
        columnNames = ["CHROM", "POS", "REF", "ALT", "DBSNP", "EA_AC", "AA_AC", "TAC", "MAF", "GTS", "EA_GTC", "AA_GTC",
                       "GTC", "DP", "FG", "GM", "AA", "AAC", "PP", "CDP", "PH", "CP", "CG", "GL", "GS", "CA",
                       "EXOME_CHIP", "GWAS_PUBMED"]
        annotationColumnNames = ["CHROM", "POS", "REF", "ALT", "DBSNP"]
        matchMode = "exact"

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        datasourceFilename = datasourceBuilder._createDatabase(destDir=destDir, ds_file=dsFile,
                                                               index_column_names=indexColumnNames,
                                                               column_names=columnNames, ds_match_mode=matchMode,
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
        configFilename = os.path.join("out", "esp_coverage.config")
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
        """

        """
        dsFile = os.path.join("testdata", "ESP6500SI-V2.chr1.snps_indels.head.25.txt")
        destDir = "out"
        datasourceFilename = "ESP6500SI-V2.chr1.snps_indels.head.25.tabix_indexed.txt.gz"
        indexColumnNames = "CHROM,POS,POS,REF,ALT"
        columnNames = "CHROM,POS,REF,ALT,DBSNP,EA_AC,AA_AC,TAC,MAF,GTS,EA_GTC,AA_GTC,GTC,DP,FG,GM,AA,AAC,PP,CDP,PH,CP,CG,GL,GS,CA,EXOME_CHIP,GWAS_PUBMED"
        configFilename = "out/esp_coverage.config"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        annotationColumnNames = "DBSNP,EA_GTC,DP"

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        datasourceBuilder.createDatasource(destDir, dsFile, indexColumnNames, configFilename, dataSourceType,
                                           dataSourceName, dataSourceVersion, dataSourceMatchMode,
                                           annotationColumnNames, DatasourceInstallUtils.getIndexCols(dataSourceType,
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

        self.assertEqual(configParser.get("data_types", "EA_GTC"), "String",
                         "Expected EA_GTC data type is %s but was %s."
                         % ("String", configParser.get("data_types", "EA_GTC")))
        self.assertEqual(configParser.get("data_types", "DP"), "Integer",
                         "Expected DP data type is %s but was %s."
                         % ("Integer", configParser.get("data_types", "DP")))

    def testCreateDatasourceWithMissingValues(self):
        """

        """
        dsFile = os.path.join("testdata", "ESP6500SI-V2.chr1.snps_indels.head.25.missing.txt")
        destDir = "out"
        datasourceFilename = "ESP6500SI-V2.chr1.snps_indels.head.25.missing.tabix_indexed.txt.gz"
        indexColumnNames = "CHROM,POS,POS"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        annotationColumnNames = "EA_GTC,DP"
        configFilename = os.path.join("out", "esp_coverage.missing.config")

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        datasourceBuilder.createDatasource(destDir, dsFile, indexColumnNames, configFilename, dataSourceType, dataSourceName,
                                           dataSourceVersion, dataSourceMatchMode, annotationColumnNames,
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
        """

        """
        dsFile = os.path.join("testdata", "ESP6500SI-V2.chr1.snps_indels.head.25.txt")
        destDir = "out"
        indexColumnNames = "CHROM,POS,POS"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        annotationColumnNames = "EA_GTC,DP"
        configFilename = os.path.join("out", "esp_coverage.missing.config")

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        try:
            datasourceBuilder.createDatasource(destDir, dsFile, indexColumnNames, configFilename,
                                               dataSourceType, dataSourceName, dataSourceVersion, dataSourceMatchMode,
                                               annotationColumnNames,
                                               DatasourceInstallUtils.getIndexCols(dataSourceType, indexColumnNames))
        except InputMismatchException:
            pass

    def testCreateDatasourceWithMissingAnnotationColumns(self):
        """

        """
        dsFile = os.path.join("testdata", "ESP6500SI-V2.chr1.snps_indels.head.25.txt")
        destDir = "out"
        indexColumnNames = "CHROM,POS,POS"
        dataSourceType = "indexed_tsv"
        dataSourceName = "ESP"
        dataSourceVersion = "6500SI-V2"
        dataSourceMatchMode = "overlap"
        annotationColumnNames = "EA_GTC,DP,ESP_DBSNP"
        configFilename = os.path.join("out", "esp_coverage.missing.config")

        datasourceBuilder = TabixIndexedTsvDatasourceCreator()
        try:
            datasourceBuilder.createDatasource(destDir, dsFile, indexColumnNames, configFilename, dataSourceType,
                                               dataSourceName, dataSourceVersion, dataSourceMatchMode,
                                               annotationColumnNames,
                                               DatasourceInstallUtils.getIndexCols(dataSourceType, indexColumnNames))
        except ValueError:
            pass