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
        cur_annotation = Annotation(value="75.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithSNPOverlapMatch(self):
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

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingSNPAvgMatch(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = os.path.join(*["testdata", "small_esp_coverage_avg_ds", "hg19"])
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedTsvDirName, "small_esp_coverage_avg_ds.config"), tabixIndexedTsvDirName)

        m1 = MutationData()
        m1.chr = "X"
        m1.start = "100075333"
        m1.end = "100075333"

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String", description="",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingSNPOverlapMatch(self):
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

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithIndelExactMatch(self):
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
        cur_annotation = Annotation(value="83.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X|X|X|X|X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithIndelOverlapMatch(self):
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
        cur_annotation = Annotation(value="79.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="692.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Chromosome")
        cur_annotation = Annotation(value="X|X|X", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPCoverageAnnotationWithMissingIndelOverlapMatch(self):
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
        cur_annotation = Annotation(value="91.25", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")
