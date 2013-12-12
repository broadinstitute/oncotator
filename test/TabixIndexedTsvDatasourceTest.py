import unittest
import logging
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.Annotation import Annotation
from oncotator.MutationData import MutationData
from oncotator.utils.TagConstants import TagConstants
from TestUtils import TestUtils

TestUtils.setupLogging(__file__, __name__)


class TabixIndexedTsvDatasourceTest(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def testESPAnnotationWithMissingMutation(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "123456"
        m1.end = "123456"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalSamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AvgEAsampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPAnnotationWithExistingMutation(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "30348"
        m1.end = "30348"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="1.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalSamplesCovered")
        cur_annotation = Annotation(value="5", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        cur_annotation = Annotation(value="1.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AvgEAsampleReadDepth")
        cur_annotation = Annotation(value="1.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="2", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPAnnotationWithMissingIndel(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "10000"
        m1.end = "10005"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalSamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AvgEAsampleReadDepth")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testESPAnnotationWithExistingIndel(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceFactory.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "69428"
        m1.end = "69432"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        m1_annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        cur_annotation = Annotation(value="139.0|142.0|144.0|147.0|150.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalSamplesCovered")
        cur_annotation = Annotation(value="5335|5350|5360|5379|5400", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        cur_annotation = Annotation(value="110.0|113.0|115.0|117.0|119.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AvgEAsampleReadDepth")
        cur_annotation = Annotation(value="94.0|96.0|98.0|100.0|101.0", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        cur_annotation = Annotation(value="1911|1914|1917|1920|1925", datasourceName="ESP", dataType="String",
                                    description="", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")