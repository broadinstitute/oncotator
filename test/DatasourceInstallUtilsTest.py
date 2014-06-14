# LICENSE_GOES_HERE
import tempfile
import unittest
import logging
from oncotator.utils.ConfigUtils import ConfigUtils
from oncotator.utils.install.DatasourceInstallUtils import DatasourceInstallUtils
import os
import vcf
from oncotator.DatasourceFactory import DatasourceFactory
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
        configFilename = os.path.join(*[destDir, "1000Genomes.config"])
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

        self.assertTrue(os.path.exists(os.path.join(*[tmpDir, datasourceFoldername, genomeBuild + ".md5"])),
                        "No md5 file was generated.")

        # Data source was created correctly
        tabixIndexedFilename = os.path.join(*[destDir, "example.tabix_indexed.vcf.gz"])
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
        destDir = os.path.join(*[tmpDir, datasourceFoldername, genomeBuild])
        os.makedirs(destDir)

        DatasourceInstallUtils.create_datasource(destDir=destDir, ds_file=datasourceFilename,
                                                 ds_foldername=datasourceFoldername, ds_name=datasourceName,
                                                 ds_type=datasourceType, ds_version=datasourceVersion,
                                                 index_columns=indexColumnNames,
                                                 ds_annotation_columns=annotationColumnNames)

        datasourceFilename = "ESP6500SI-V2.chr1.snps_indels.head.25.tabix_indexed.txt.gz"
        configFilename = os.path.join(*[destDir, "1000Genomes.config"])
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

        self.assertTrue(os.path.exists(os.path.join(*[tmpDir, datasourceFoldername, genomeBuild + ".md5"])),
                        "No md5 file was generated.")

        datasource = DatasourceFactory.createDatasource(configFilename, destDir)

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
        destDir = os.path.join(*[tmpDir, datasourceFoldername, genomeBuild])
        os.makedirs(destDir)

        DatasourceInstallUtils.create_datasource(destDir, datasourceFilename, datasourceFoldername, datasourceName,
                                                 datasourceType, datasourceVersion, genomicPositionColumnNames)

        datasourceFilename = "oreganno_trim.hg19.txt"
        configFilename = os.path.join(*[destDir, "ORegAnno.config"])
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

        self.assertTrue(os.path.exists(os.path.join(*[tmpDir, datasourceFoldername, genomeBuild + ".md5"])),
                        "No md5 file was generated.")

        MutUtils.removeDir(tmpDir)

if __name__ == '__main__':
    unittest.main()