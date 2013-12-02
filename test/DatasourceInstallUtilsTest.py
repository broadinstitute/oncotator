"""
# By downloading the PROGRAM you agree to the following terms of use:
#
# BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
# FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
#
# This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
# WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
# WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
# NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
#
# 1. DEFINITIONS
# 1.1	"PROGRAM" shall mean copyright in the object code and source code known as Oncotator and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.
#
# 2. LICENSE
# 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
# LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
# The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
# 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
# 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
#
# 3. OWNERSHIP OF INTELLECTUAL PROPERTY
# LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
#
# Copyright 2012 Broad Institute, Inc.
# Notice of attribution:  The Oncotator program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc.
#
# LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
#
# 4. INDEMNIFICATION
# LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
#
# 5. NO REPRESENTATIONS OR WARRANTIES
# THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
# IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
#
# 6. ASSIGNMENT
# This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
#
# 7. MISCELLANEOUS
# 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
# 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
# 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
# 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
# 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
# 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
# 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
#"""
import tempfile
import unittest
import logging
from oncotator.utils.ConfigUtils import ConfigUtils
from oncotator.utils.install.DatasourceInstallUtils import DatasourceInstallUtils
import os
import string
import vcf
from oncotator.DatasourceCreator import DatasourceCreator
from oncotator.Annotation import Annotation
from oncotator.MutationData import MutationData
from oncotator.utils.TagConstants import TagConstants

from TestUtils import TestUtils
from oncotator.utils.MutUtils import MutUtils

TestUtils.setupLogging(__file__, __name__)


class DatasourceInstallUtilsTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()

    def tearDown(self):
        pass

    def testCreateIndexedVcfDatasource(self):
        datasourceFilename = "testdata/vcf/example.vcf"
        datasourceFoldername = "1000Genomes"
        datasourceName = "1000Genomes"
        datasourceType = "indexed_vcf"
        datasourceVersion = "V4.1"
        genomeBuild = "hg19"
        tmpDir = tempfile.mkdtemp()
        destDir = os.path.join(*[tmpDir, datasourceFoldername, genomeBuild])
        os.makedirs(destDir)

        DatasourceInstallUtils.create_datasource(destDir, datasourceFilename, datasourceFoldername, datasourceName,
                                                 datasourceType, datasourceVersion)

        datasourceFilename = "example.tabix_indexed.vcf.gz"
        configFilename = string.join([destDir, "1000Genomes.config"], os.sep)
        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")

        self.assertEqual(configParser.get("general", "type"), datasourceType,
                         "Expected data source type is %s but was %s."
                         % (datasourceType, configParser.get("general", "type")))
        self.assertEqual(configParser.get("general", "src_file"), datasourceFilename,
                         "Expected data source src_file is %s but was %s."
                         % (datasourceFilename, configParser.get("general", "src_file")))
        self.assertEqual(configParser.get("general", "title"), datasourceName,
                         "Expected data source title is %s but was %s."
                         % (datasourceName, configParser.get("general", "title")))
        self.assertEqual(configParser.get("general", "version"), datasourceVersion,
                         "Expected data source version is %s but was %s."
                         % (datasourceVersion, configParser.get("general", "version")))

        self.assertTrue(os.path.exists(string.join([tmpDir, datasourceFoldername, genomeBuild + ".md5"], os.sep)),
                        "No md5 file was generated.")

        # Data source was created correctly
        tabixIndexedFilename = string.join([destDir, "example.tabix_indexed.vcf.gz"], os.sep)
        self.assertTrue(os.path.exists(tabixIndexedFilename), "No index file was generated.")

        vcfReader = vcf.Reader(filename=tabixIndexedFilename, compressed=True, strict_whitespace=True)
        vcfRecords = vcfReader.fetch(chrom=20, start=1230237, end=1230237)
        for vcfRecord in vcfRecords:
            self.assertEqual(vcfRecord.INFO["NS"], 3, "Expected %s but got %s." % (3, vcfRecord.INFO["NS"]))
            self.assertEqual(vcfRecord.INFO["DP"], 13, "Expected %s but got %s." % (13, vcfRecord.INFO["DP"]))

        MutUtils.removeDir(tmpDir)

    def testCreateIndexedTsvDatasource(self):
        datasourceFilename = "testdata/ESP6500SI-V2.chr1.snps_indels.head.25.txt"
        datasourceFoldername = "1000Genomes"
        datasourceName = "1000Genomes"
        datasourceType = "indexed_tsv"
        datasourceVersion = "V4.1"
        genomeBuild = "hg19"
        indexColumnNames = "CHROM,POS,POS"
        columnNames = "CHROM,POS,REF,ALT,DBSNP,EA_AC,AA_AC,TAC,MAF,GTS,EA_GTC,AA_GTC,GTC,DP,FG,GM,AA,AAC,PP,CDP,PH,CP,CG,GL,GS,CA,EXOME_CHIP,GWAS_PUBMED"
        annotationColumnNames = "DBSNP,EA_AC,AA_AC,TAC"

        tmpDir = tempfile.mkdtemp()
        destDir = string.join([tmpDir, datasourceFoldername, genomeBuild], os.sep)
        os.makedirs(destDir)

        DatasourceInstallUtils.create_datasource(destDir, datasourceFilename, datasourceFoldername, datasourceName,
                                                 datasourceType, datasourceVersion, indexColumnNames, columnNames,
                                                 annotationColumnNames)

        datasourceFilename = "ESP6500SI-V2.chr1.snps_indels.head.25.tabix_indexed.txt.gz"
        configFilename = string.join([destDir, "1000Genomes.config"], os.sep)
        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "column_names"),
                        "column_names option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "annotation_column_names"),
                        "annotation_column_names option is missing in general section.")

        self.assertEqual(configParser.get("general", "type"), datasourceType,
                         "Expected data source type is %s but was %s."
                         % (datasourceType, configParser.get("general", "type")))
        self.assertEqual(configParser.get("general", "src_file"), datasourceFilename,
                         "Expected data source src_file is %s but was %s."
                         % (datasourceFilename, configParser.get("general", "src_file")))
        self.assertEqual(configParser.get("general", "title"), datasourceName,
                         "Expected data source title is %s but was %s."
                         % (datasourceName, configParser.get("general", "title")))
        self.assertEqual(configParser.get("general", "version"), datasourceVersion,
                         "Expected data source version is %s but was %s."
                         % (datasourceVersion, configParser.get("general", "version")))
        self.assertEqual(configParser.get("general", "column_names"), columnNames,
                         "Expected data source column names is %s but was %s."
                         % (columnNames, configParser.get("general", "column_names")))
        self.assertEqual(configParser.get("general", "annotation_column_names"), annotationColumnNames,
                         "Expected data source annotation column names is %s but was %s."
                         % (annotationColumnNames, configParser.get("general", "annotation_column_names")))

        self.assertTrue(os.path.exists(string.join([tmpDir, datasourceFoldername, genomeBuild + ".md5"], os.sep)),
                        "No md5 file was generated.")

        datasource = DatasourceCreator.createDatasource(configFilename, destDir)

        m1 = MutationData()
        m1.chr = "1"
        m1.start = "802177"
        m1.end = "802177"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = datasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("1000Genomes_AA_AC")
        cur_annotation = Annotation(value="2,866", datasourceName="1000Genomes", dataType="String",
                                    description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        annotationNames = ["1000Genomes_CHROM", "1000Genomes_POS", "1000Genomes_REF", "1000Genomes_ALT",
                           "1000Genomes_GWAS_PUBMED"]
        for annotationName in annotationNames:
            self.assertTrue(annotationName not in m1_annotated, "m1_annotated was annotated with %s." % annotationName)

        annotationNames = ["1000Genomes_DBSNP", "1000Genomes_EA_AC", "1000Genomes_AA_AC", "1000Genomes_TAC"]
        for annotationName in annotationNames:
            self.assertTrue(annotationName in m1_annotated, "m1_annotated was not annotated with %s value."
                                                            % annotationName)

        MutUtils.removeDir(tmpDir)

    def testCreateGPTsvDatasource(self):
        """


        """
        datasourceFilename = "testdata/small_genome_position_tsv_ds/oreganno_trim.hg19.txt"
        datasourceType = "gp_tsv"
        datasourceName = "ORegAnno"
        datasourceFoldername = "ORegAnno"
        datasourceVersion = "UCSC Track"
        genomeBuild = "hg19"
        genomicPositionColumnNames = "hg19.oreganno.chrom,hg19.oreganno.chromStart,hg19.oreganno.chromEnd"

        tmpDir = tempfile.mkdtemp()
        destDir = string.join([tmpDir, datasourceFoldername, genomeBuild], os.sep)
        os.makedirs(destDir)

        DatasourceInstallUtils.create_datasource(destDir, datasourceFilename, datasourceFoldername, datasourceName,
                                                 datasourceType, datasourceVersion, genomicPositionColumnNames)

        datasourceFilename = "oreganno_trim.hg19.txt"
        configFilename = string.join([destDir, "ORegAnno.config"], os.sep)
        configParser = ConfigUtils.createConfigParser(configFilename)
        self.assertTrue(configParser.has_section("general"), "general section is missing.")
        self.assertTrue(configParser.has_option("general", "type"), "type option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "src_file"),
                        "src_file option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "title"), "title option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "version"), "version option is missing in general section.")
        self.assertTrue(configParser.has_option("general", "genomic_position_cols"),
                        "genomic_position_cols option is missing in general section.")

        self.assertEqual(configParser.get("general", "type"), datasourceType,
                         "Expected data source type is %s but was %s."
                         % (datasourceType, configParser.get("general", "type")))
        self.assertEqual(configParser.get("general", "src_file"), datasourceFilename,
                         "Expected data source src_file is %s but was %s."
                         % (datasourceFilename, configParser.get("general", "src_file")))
        self.assertEqual(configParser.get("general", "title"), datasourceName,
                         "Expected data source title is %s but was %s."
                         % (datasourceName, configParser.get("general", "title")))
        self.assertEqual(configParser.get("general", "version"), datasourceVersion,
                         "Expected data source version is %s but was %s."
                         % (datasourceVersion, configParser.get("general", "version")))
        self.assertEqual(configParser.get("general", "genomic_position_cols"), genomicPositionColumnNames,
                         "Expected data source genomic_position_cols is %s but was %s."
                         % (genomicPositionColumnNames, configParser.get("general", "genomic_position_cols")))

        self.assertTrue(os.path.exists(string.join([tmpDir, datasourceFoldername, genomeBuild + ".md5"], os.sep)),
                        "No md5 file was generated.")

        MutUtils.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()