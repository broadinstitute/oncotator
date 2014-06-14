# LICENSE_GOES_HERE
import unittest
import logging
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.Annotation import Annotation
from oncotator.MutationData import MutationData
from oncotator.utils.TagConstants import TagConstants
from TestUtils import TestUtils
import os

TestUtils.setupLogging(__file__, __name__)


class TabixIndexedTsvDatasourceTest(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def testESPCoverageAnnotationWithSNPExactMatch(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_exact_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_exact_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075334"
        m1.end = "100075334"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="75.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithSNPAvgMatch(self):
        """
        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_avg_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_avg_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075334"
        m1.end = "100075334"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="75.0", datasourceName="ESP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692.0", datasourceName="ESP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithSNPOverlapMatch(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_overlap_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_overlap_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075334"
        m1.end = "100075334"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="75.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingSNPExactMatch(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_exact_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_exact_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075333"
        m1.end = "100075333"

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingSNPAvgMatch(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_avg_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_avg_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075333"
        m1.end = "100075333"

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingSNPOverlapMatch(self):
        """
        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_overlap_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_overlap_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075333"
        m1.end = "100075333"

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithIndelExactMatch(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_exact_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_exact_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075337"
        m1.end = "100075340"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="85.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithIndelAvgMatch(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_avg_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_avg_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075337"
        m1.end = "100075340"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="83.0", datasourceName="ESP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692.0", datasourceName="ESP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X|X|X|X|X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithIndelOverlapMatch(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_overlap_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_overlap_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075337"
        m1.end = "100075340"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="81.0|85.0|84.0|84.0|81.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692|692|692|692|692", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X|X|X|X|X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingIndelExactMatch(self):
        """


        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_exact_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_exact_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075300"
        m1.end = "100075336"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingIndelAvgMatch(self):
        """


        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_avg_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_avg_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075300"
        m1.end = "100075336"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="79.0", datasourceName="ESP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692.0", datasourceName="ESP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X|X|X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingIndelOverlapMatch(self):
        """


        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_overlap_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_overlap_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075300"
        m1.end = "100075336"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="75.0|81.0|81.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692|692|692", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X|X|X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingAnnotationValuesIndelAvgMatch(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_avg_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_avg_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075350"
        m1.end = "100075356"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        cur_annotation = Annotation(value="91.25", datasourceName="ESP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testdbNSFPAnnotationWithExactMatch(self):  # SNPs only
        """

        """
        self.logger.info("Initializing dbNSFP")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "dbNSFP_chr1_6vars_exact_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "dbNSFP_chr1_6vars_exact_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "1"
        m1.start = "35138"
        m1.end = "35138"
        m1.ref_allele = "T"
        m1.alt_allele = "G"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("dbNSFP_codonpos")
        cur_annotation = Annotation(value="3", datasourceName="dbNSFP", dataType="Integer",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_refcodon")
        cur_annotation = Annotation(value="TAA", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_cds_strand")
        cur_annotation = Annotation(value="-", datasourceName="dbNSFP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testdbNSFPAnnotationWithMissingExactMatch(self):  # SNPs only
        """

        """
        self.logger.info("Initializing dbNSFP")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "dbNSFP_chr1_6vars_exact_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "dbNSFP_chr1_6vars_exact_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "1"
        m1.start = "35138"
        m1.end = "35138"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("dbNSFP_codonpos")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="Integer",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_refcodon")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_cds_strand")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testdbNSFPAnnotationWithAvgMatch(self):  # SNPs only
        """

        """
        self.logger.info("Initializing dbNSFP")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "dbNSFP_chr1_chr3_100vars_avg_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "dbNSFP_chr1_chr3_100vars_avg_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "1"
        m1.start = "35138"
        m1.end = "35139"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("dbNSFP_codonpos")
        cur_annotation = Annotation(value="2.5", datasourceName="dbNSFP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_refcodon")
        cur_annotation = Annotation(value="TAA|TAA|TAA|TAA", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_cds_strand")
        cur_annotation = Annotation(value="-|-|-|-", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testdbNSFPAnnotationWithMissingAvgMatch(self):  # SNPs only
        """

        """
        self.logger.info("Initializing dbNSFP")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "dbNSFP_chr1_chr3_100vars_avg_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "dbNSFP_chr1_chr3_100vars_avg_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "1"
        m1.start = "35137"
        m1.end = "35137"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("dbNSFP_codonpos")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="Float",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_refcodon")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_cds_strand")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testdbNSFPAnnotationWithOverlapMatch(self):  # SNPs only
        """

        """
        self.logger.info("Initializing dbNSFP")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "dbNSFP_chr1_chr3_100vars_overlap_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "dbNSFP_chr1_chr3_100vars_overlap_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "1"
        m1.start = "35138"
        m1.end = "35139"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("dbNSFP_codonpos")
        cur_annotation = Annotation(value="3|3|2|2", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_refcodon")
        cur_annotation = Annotation(value="TAA|TAA|TAA|TAA", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_cds_strand")
        cur_annotation = Annotation(value="-|-|-|-", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testdbNSFPAnnotationWithMissingOverlapMatch(self):  # SNPs only
        """

        """
        self.logger.info("Initializing dbNSFP")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "dbNSFP_chr1_chr3_100vars_overlap_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "dbNSFP_chr1_chr3_100vars_overlap_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "1"
        m1.start = "35136"
        m1.end = "35137"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("dbNSFP_codonpos")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_refcodon")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_cds_strand")
        cur_annotation = Annotation(value="", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testdbNSFPNoRefAltAnnotationWithExactMatch(self):
        """

        """
        self.logger.info("Initializing dbNSFP")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "dbNSFP_chr1_chr3_100vars_exact_no_ref_alt_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "dbNSFP_chr1_chr3_100vars_exact_no_ref_alt_ds.config"),
            tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "1"
        m1.start = "35140"
        m1.end = "35140"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("dbNSFP_codonpos")
        cur_annotation = Annotation(value="1|1|1", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_refcodon")
        cur_annotation = Annotation(value="TAA|TAA|TAA", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("dbNSFP_cds_strand")
        cur_annotation = Annotation(value="-|-|-", datasourceName="dbNSFP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")