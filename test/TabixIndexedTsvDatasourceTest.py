import unittest
import logging
from oncotator.DatasourceCreator import DatasourceCreator

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

    def _validateAnnotation(self, annotation, dataType, ds, desc, num, tags, val):
        self.assertEqual(annotation.getDataType(), dataType, "Expected data type is %s but was %s." %
                                                             (dataType, annotation.getDataType()))
        self.assertEqual(annotation.getDescription(), desc, "Expected description is %s but was %s." %
                                                            (desc, annotation.getDescription()))
        self.assertEqual(annotation.getDatasource(), ds, "Expected data type is %s but was %s." %
                                                         (ds, annotation.getDatasource()))
        self.assertEqual(annotation.getNumber(), num, "Expected num is %s but was %s." % (num, annotation.getNumber()))
        self.assertEqual(annotation.getValue(), val, "Expected val is %s but was %s." % (val, annotation.getValue()))
        self.assertEqual(annotation.getTags(), tags, "Expected tags is %s but was %s" % (tags, annotation.getTags()))

    def testCreateESPCoverageDataSource(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceCreator.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)

    def testESPAnnotationWithMissingMutation(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceCreator.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "123456"
        m1.end = "123456"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_TotalSamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_TotalEAsamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_AvgEAsampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

    def testESPAnnotationWithExistingMutation(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceCreator.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "30348"
        m1.end = "30348"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "1.0")

        annotation = m1_annotated.getAnnotation("ESP_TotalSamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "5")

        annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "1.0")

        annotation = m1_annotated.getAnnotation("ESP_TotalEAsamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "3")

        annotation = m1_annotated.getAnnotation("ESP_AvgEAsampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "1.0")

        annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "2")

    def testESPAnnotationWithMissingIndel(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceCreator.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "10000"
        m1.end = "10005"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_TotalSamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_TotalEAsamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_AvgEAsampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None, [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

    def testESPAnnotationWithExistingIndel(self):
        self.logger.info("Initializing ESP6500SI-V2 Coverage")
        tabixIndexedTsvDirName = "testdata/small_esp_coverage_ds"
        tabixIndexedTsvDatasource = DatasourceCreator.createDatasource(tabixIndexedTsvDirName + "/esp_coverage.config",
                                                                       tabixIndexedTsvDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "69428"
        m1.end = "69432"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedTsvDatasource.annotate_mutation(m1)
        annotation = m1_annotated.getAnnotation("ESP_AvgAAsampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "139.0|142.0|144.0|147.0|150.0")

        annotation = m1_annotated.getAnnotation("ESP_TotalSamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "5335|5350|5360|5379|5400")

        annotation = m1_annotated.getAnnotation("ESP_AvgSampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "110.0|113.0|115.0|117.0|119.0")

        annotation = m1_annotated.getAnnotation("ESP_TotalEAsamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "3424|3436|3443|3459|3475")

        annotation = m1_annotated.getAnnotation("ESP_AvgEAsampleReadDepth")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "94.0|96.0|98.0|100.0|101.0")

        annotation = m1_annotated.getAnnotation("ESP_TotalAAsamplesCovered")
        self._validateAnnotation(annotation, "String", "ESP", "", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "1911|1914|1917|1920|1925")

