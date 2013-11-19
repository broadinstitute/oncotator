import unittest
import logging
from oncotator.DatasourceCreator import DatasourceCreator

from oncotator.MutationData import MutationData
from TestUtils import TestUtils
import vcf
from oncotator.utils.TagConstants import TagConstants

TestUtils.setupLogging(__file__, __name__)


class TabixIndexedVcfDatasourceTest(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def testTags(self):
        self.logger.info("Initializing ESP6500SI-V2")
        tabixIndexedVcfDirName = "testdata/small_esp_ds"
        tabixIndexedVcfDatasource = DatasourceCreator.createDatasource(tabixIndexedVcfDirName + "/esp.config",
                                                                       tabixIndexedVcfDirName)
        tagsDict = tabixIndexedVcfDatasource._determine_tags()
        for ID in tagsDict:
            tags = tagsDict[ID]
            self.assertTrue(len(tags) == 2, "The length of tags is not 2 but %s." % len(tags))
            self.assertTrue(TagConstants.INFO in tags, "INFO tag is missing for %s." % ID)
            self.assertTrue(TagConstants.NOT_SPLIT in tags, "NOT_SPLIT tag is missing for %s." % ID)

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

    def testESPAnnotationWithMissingMutation(self):
        self.logger.info("Initializing ESP6500SI-V2")
        tabixIndexedVcfDirName = "testdata/small_esp_ds"
        tabixIndexedVcfDatasource = DatasourceCreator.createDatasource(tabixIndexedVcfDirName + "/esp.config",
                                                                       tabixIndexedVcfDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "123456"
        m1.end = "123456"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)
        annotation = m1_annotated.getAnnotation("ESP_PP")
        self._validateAnnotation(annotation, "String", "ESP", "proteinPosition", None, [TagConstants.INFO,
                                                                                        TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_CG")
        self._validateAnnotation(annotation, "Float", "ESP", "consScoreGERP", 1, [TagConstants.INFO,
                                                                                  TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_MAF")
        self._validateAnnotation(annotation, "Float", "ESP",
                                 "Minor Allele Frequency in percent in the order of EA,AA,All", 3,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], ",,")

    def testESPAnnotationWithExistingMutation(self):
        self.logger.info("Initializing ESP6500SI-V2")
        tabixIndexedVcfDirName = "testdata/small_esp_ds"
        tabixIndexedVcfDatasource = DatasourceCreator.createDatasource(tabixIndexedVcfDirName + "/esp.config",
                                                                       tabixIndexedVcfDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "69594"
        m1.end = "69594"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)
        annotation = m1_annotated.getAnnotation("ESP_PP")
        self._validateAnnotation(annotation, "String", "ESP", "proteinPosition", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "168/306")

        annotation = m1_annotated.getAnnotation("ESP_CG")
        self._validateAnnotation(annotation, "Float", "ESP", "consScoreGERP", 1,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "-0.2")

        annotation = m1_annotated.getAnnotation("ESP_MAF")
        self._validateAnnotation(annotation, "Float", "ESP",
                                 "Minor Allele Frequency in percent in the order of EA,AA,All", 3,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "0.0323,0.0,0.0205")

    def testESPAnnotationWithMissingIndel(self):
        self.logger.info("Initializing ESP6500SI-V2")
        tabixIndexedVcfDirName = "testdata/small_esp_ds"
        tabixIndexedVcfDatasource = DatasourceCreator.createDatasource(tabixIndexedVcfDirName + "/esp.config",
                                                                       tabixIndexedVcfDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "10000"
        m1.end = "10005"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)
        annotation = m1_annotated.getAnnotation("ESP_PP")
        self._validateAnnotation(annotation, "String", "ESP", "proteinPosition", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_CG")
        self._validateAnnotation(annotation, "Float", "ESP", "consScoreGERP", 1,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "")

        annotation = m1_annotated.getAnnotation("ESP_MAF")
        self._validateAnnotation(annotation, "Float", "ESP",
                                 "Minor Allele Frequency in percent in the order of EA,AA,All", 3,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], ",,")

    def testESPAnnotationWithExistingIndel(self):
        self.logger.info("Initializing ESP6500SI-V2")
        tabixIndexedVcfDirName = "testdata/small_esp_ds"
        tabixIndexedVcfDatasource = DatasourceCreator.createDatasource(tabixIndexedVcfDirName + "/esp.config",
                                                                       tabixIndexedVcfDirName)
        m1 = MutationData()
        m1.chr = "1"
        m1.start = "69428"
        m1.end = "69496"
        m1.ref_allele = "T"
        m1.alt_allele = "C"

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)
        annotation = m1_annotated.getAnnotation("ESP_PP")
        self._validateAnnotation(annotation, "String", "ESP", "proteinPosition", None,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "113/306|129/306|136/306")

        annotation = m1_annotated.getAnnotation("ESP_CG")
        self._validateAnnotation(annotation, "Float", "ESP", "consScoreGERP", 1,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT], "0.9|2.3|2.3")

        annotation = m1_annotated.getAnnotation("ESP_MAF")
        self._validateAnnotation(annotation, "Float", "ESP",
                                 "Minor Allele Frequency in percent in the order of EA,AA,All", 3,
                                 [TagConstants.INFO, TagConstants.NOT_SPLIT],
                                 "4.5707,0.3663,3.0647|0.0285,0.0,0.0183|0.0296,0.604,0.2364")