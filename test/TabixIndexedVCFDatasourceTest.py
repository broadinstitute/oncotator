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

import unittest
import logging
import os
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.Annotation import Annotation
from oncotator.MutationData import MutationData
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

    def _createMut(self, chrom, startPos, endPos, ref, alt, build):
        """

        :param chrom:
        :param startPos:
        :param endPos:
        :param ref:
        :param alt:
        :param build:
        :return:
        """
        mut = MutationData(chrom, int(startPos), int(endPos), ref, alt, build)

        varType = MutUtils.determineVariantType(mut)

        if varType == "snp":  # Snps
            mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue="")
        if varType == "del":  # deletion
            preceding_bases, updated_ref_allele, updated_start, updated_end =\
                MutUtils.retrievePrecedingBasesForDeletions(mut)
            mut.ref_allele = updated_ref_allele
            mut.alt_allele = "-"
            mut.start = updated_start
            mut.end = updated_end
            mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME,
                                 annotationValue=preceding_bases)
        elif varType == "ins":  # insertion
            preceding_bases, updated_alt_allele, updated_start, updated_end = \
                MutUtils.retrievePrecedingBasesForInsertions(mut)
            mut.ref_allele = "-"
            mut.alt_allele = updated_alt_allele
            mut.start = updated_start
            mut.end = updated_end
            mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME,
                                 annotationValue=preceding_bases)
        return mut

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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False", datasourceName="ESP", dataType="Flag",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=0)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "11"
        start = "17330"
        end = "17330"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False", datasourceName="ESP", dataType="Flag",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=0)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="Float",
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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False", datasourceName="ESP", dataType="Flag",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=0)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="Float",
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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.333,0.667", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="2,4", datasourceName="ESP", dataType="String",
                                    description="Allele Count", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="False", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Z")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
                                    description="A random variable, Z", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=3)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "20"
        start = "1230237"
        end = "1230237"
        ref_allele = "T"
        alt_allele = "A"
        build = "hg19"
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.5,0.5|0.017", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="3,3|1", datasourceName="ESP", dataType="String",
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
        cur_annotation = Annotation(value="1.0,2.0,3.0|4.0,5.0,|0.0,2.0,1.0", datasourceName="ESP", dataType="String",
                                    description="A random variable, Z", tags=[TagConstants.INFO,
                                                                              TagConstants.NOT_SPLIT], number=3)
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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "11"
        start = "17330"
        end = "17330"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

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
        cur_annotation = Annotation(value="False", datasourceName="ESP", dataType="String",
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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="0.339", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_AC")
        cur_annotation = Annotation(value="2.33333333333", datasourceName="ESP", dataType="Float",
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
        cur_annotation = Annotation(value="1.66666666667,2.33333333333,2.25", datasourceName="ESP", dataType="Float",
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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        chrom = "11"
        start = "17330"
        end = "17330"
        ref_allele = "T"
        alt_allele = "C"
        build = "hg19"
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
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
        m1 = self._createMut(chrom, start, end, ref_allele, alt_allele, build)

        m1_annotated = tabixIndexedVcfDatasource.annotate_mutation(m1)

        m1_annotation = m1_annotated.getAnnotation("ESP_AF")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="Float",
                                    description="Allele Frequency", tags=[TagConstants.INFO, TagConstants.SPLIT],
                                    number=-1)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_X")
        cur_annotation = Annotation(value=",", datasourceName="ESP", dataType="String",
                                    description="A random variable, X", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_H2")
        cur_annotation = Annotation(value="", datasourceName="ESP", dataType="String",
                                    description="HapMap2 membership", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=None)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")

        m1_annotation = m1_annotated.getAnnotation("ESP_Y")
        cur_annotation = Annotation(value=",,", datasourceName="ESP", dataType="String",
                                    description="A random variable, Y", tags=[TagConstants.INFO, TagConstants.NOT_SPLIT],
                                    number=-2)
        self.assertTrue(m1_annotation.isEqual(cur_annotation), "Annotations do not match.")
