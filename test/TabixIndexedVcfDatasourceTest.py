# LICENSE_GOES_HERE

import unittest
import logging
import os
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.Annotation import Annotation
from TestUtils import TestUtils
from oncotator.utils.TagConstants import TagConstants
from oncotator.utils.MutUtils import MutUtils

TestUtils.setupLogging(__file__, __name__)


class TabixIndexedVcfDatasourceTest(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def testTags(self):
        """

        """
        self.logger.info("Initializing ESP6500SI-V2")
        tabixIndexedVcfDirName = os.path.join(*["testdata", "small_esp_ds"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "esp.config"), tabixIndexedVcfDirName)

        tagsDict = tabixIndexedVcfDatasource._determine_tags()
        for ID in tagsDict:
            tags = tagsDict[ID]
            self.assertTrue(len(tags) == 2, "The length of tags is not 2 but %s." % len(tags))
            self.assertTrue(TagConstants.INFO in tags, "INFO tag is missing for %s." % ID)
            self.assertTrue(TagConstants.NOT_SPLIT in tags, "NOT_SPLIT tag is missing for %s." % ID)

    def testExampleVcfDBAnnotationWithSNPExactMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_exact", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_exact.config"), tabixIndexedVcfDirName)

        chrom = "20"
        start = "1110696"
        end = "1110696"
        ref_allele = "A"
        alt_allele = "T"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.667", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="2,4", datasourceName="ESP", dataType="Integer",
                                    description="Allele Count", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False", datasourceName="ESP", dataType="Flag",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=0)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "20"
        start = "1230237"
        end = "1230237"
        ref_allele = "T"
        alt_allele = "A"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_NS")
        cur_annotation = Annotation(value="3", datasourceName="ESP", dataType="Integer",
                                    description="Number of Samples With Data",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithIndelExactMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_exact", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_exact.config"), tabixIndexedVcfDirName)

        chrom = "21"
        start = "1234567"
        end = "1234567"
        ref_allele = "GTC"
        alt_allele = "GTCT"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.167", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="3,1", datasourceName="ESP", dataType="Integer",
                                    description="Allele Count", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False", datasourceName="ESP", dataType="Flag",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=0)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AA")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Ancestral Allele", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="Float",
                                    description="A random variable, Z", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=3)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithMissingSNPExactMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_exact", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_exact.config"), tabixIndexedVcfDirName)

        chrom = "20"
        start = "17330"
        end = "17330"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Flag",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=0)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "11"
        start = "17330"
        end = "17330"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Flag",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=0)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="A random variable, Z", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=3)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithMissingIndelExactMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_exact", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_exact.config"), tabixIndexedVcfDirName)

        chrom = "21"
        start = "1234567"
        end = "1234567"
        ref_allele = "AGTC"
        alt_allele = "A"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Flag",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=0)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="A random variable, Z", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=3)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithSNPOverlapMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_overlap", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_overlap.config"), tabixIndexedVcfDirName)

        chrom = "20"
        start = "1110696"
        end = "1110696"
        ref_allele = "A"
        alt_allele = "T"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.333|0.667", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="2,4|2,4", datasourceName="ESP", dataType="String",
                                    description="Allele Count", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False|False", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value=",,|,,", datasourceName="ESP", dataType="String",
                                    description="A random variable, Z", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=3)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "20"
        start = "1230237"
        end = "1230237"
        ref_allele = "T"
        alt_allele = "A"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_NS")
        cur_annotation = Annotation(value="3", datasourceName="ESP", dataType="String",
                                    description="Number of Samples With Data",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "11"
        start = "1230237"
        end = "1230237"
        ref_allele = "T"
        alt_allele = "A"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_NS")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Number of Samples With Data",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithIndelOverlapMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_overlap", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_overlap.config"), tabixIndexedVcfDirName)

        chrom = "4"
        start = "1234567"
        end = "1234567"
        ref_allele = "GTC"
        alt_allele = "GTCTTA"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.5|0.5|0.5", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="3,3|3,3|3,3", datasourceName="ESP", dataType="String",
                                    description="Allele Count", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False|False|False|False", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AA")
        cur_annotation = Annotation(value="T", datasourceName="ESP", dataType="String",
                                    description="Ancestral Allele", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value="6.0,7.0,1.0|1.0,2.0,3.0|1.0,2.0,3.0|4.0,5.0,", datasourceName="ESP",
                                    dataType="String", description="A random variable, Z",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT], number=3)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithMissingSNPOverlapMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_overlap", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_overlap.config"), tabixIndexedVcfDirName)

        chrom = "20"
        start = "17329"
        end = "17329"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "11"
        start = "17330"
        end = "17330"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithMissingIndelOverlapMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_overlap", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_overlap.config"), tabixIndexedVcfDirName)

        chrom = "21"
        start = "1234570"
        end = "1234570"
        ref_allele = "AGTC"
        alt_allele = "A"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithSNPAvgMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_avg", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_avg.config"), tabixIndexedVcfDirName)

        chrom = "20"
        start = "1110696"
        end = "1110696"
        ref_allele = "A"
        alt_allele = "T"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.5", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="3.0", datasourceName="ESP", dataType="Float",
                                    description="Allele Count", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False|False", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="Float",
                                    description="A random variable, Z", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=3)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "20"
        start = "1230237"
        end = "1230237"
        ref_allele = "T"
        alt_allele = "A"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_NS")
        cur_annotation = Annotation(value="3.0", datasourceName="ESP", dataType="Float",
                                    description="Number of Samples With Data",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "11"
        start = "1230237"
        end = "1230237"
        ref_allele = "T"
        alt_allele = "A"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_NS")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Number of Samples With Data",
                                    tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithIndelAvgMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_avg", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_avg.config"), tabixIndexedVcfDirName)

        chrom = "4"
        start = "1234567"
        end = "1234567"
        ref_allele = "GTC"
        alt_allele = "GTCTTA"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.5", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="3.0", datasourceName="ESP", dataType="Float",
                                    description="Allele Count", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False|False|False", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AA")
        cur_annotation = Annotation(value="T", datasourceName="ESP", dataType="String",
                                    description="Ancestral Allele", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value="2.0,3.0,3.0", datasourceName="ESP", dataType="Float",
                                    description="A random variable, Z", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=3)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithMissingSNPAvgMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_avg", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_avg.config"), tabixIndexedVcfDirName)

        chrom = "20"
        start = "17329"
        end = "17329"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "11"
        start = "17330"
        end = "17330"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

    def testExampleVcfDBAnnotationWithMissingIndelAvgMatch(self):
        """

        """
        tabixIndexedVcfDirName = os.path.join(*["testdata", "vcf_db_avg", "hg19"])
        tabixIndexedVcfDatasource = DatasourceFactory.createDatasource(
            os.path.join(tabixIndexedVcfDirName, "vcf_db_avg.config"), tabixIndexedVcfDirName)

        chrom = "21"
        start = "1234570"
        end = "1234570"
        ref_allele = "AGTC"
        alt_allele = "A"
        build = "hg19"
        m1 = MutUtils.initializeMutFromAttributes(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")
