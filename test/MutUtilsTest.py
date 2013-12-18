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

from oncotator.MutationData import MutationData
from oncotator.utils.MutUtils import MutUtils
from TestUtils import TestUtils
import vcf

TestUtils.setupLogging(__file__, __name__)


class MutUtilsTest(unittest.TestCase):
    def testProteinChange(self):
        """ Test that protein change parsing of start and end works.
        """
        # Each tuple is test, ground truth
        testInOuts = [
            ("p.K128_R130del", ['128','130']),
            ("p.W274G", ["274", "274"]),
            ("p.13_14AA>A", ["13", "14"]),
            ("p.G25_splice", ["25", "25"]),
            ("p.E813*", ["813", "813"]),
            ("p.SLPQPEQRPY59del", ["59", "59"])
        ]

        ctr = 1
        for test in testInOuts:
            result = MutUtils.extractProteinPosition(test[0])
            self.assertTrue(result != ['', ''], "Result was empty.  " + str(test[0]) + ".  ")
            self.assertTrue(result[0] == test[1][0] and result[1] == test[1][1], "Result did not match for " + str(test[0]) + ".  " + str(result) + "  GT: " + str(test[1]))
            ctr += 1
        self.assertTrue(MutUtils.extractProteinPosition("blahblah") == ['', ''])

    def testRetrieveMissingAnnotations(self):
        """ Test simple case.
        """
        m = MutationData()
        m.createAnnotation("a1", "1")
        m.createAnnotation("a2", "1")
        m.createAnnotation("a3", "1")
        m.createAnnotation("a4", "1")

        annotationNames = ["a3", "a2"]

        result = MutUtils.retrieveMissingAnnotations(m,annotationNames)

        self.assertIsNotNone(result)
        self.assertTrue(len(result) == 0, "Result was not empty: " + str(result))

        annotationNames = ["zztop", "a1", "blah", "dummy"]
        result = MutUtils.retrieveMissingAnnotations(m,annotationNames)
        self.assertTrue(result[0] == "blah", "Result was not sorted")
        self.assertTrue("blah" in result and "dummy" in result and "zztop" in result, "Incorrect elements (Truth: [zztop, blah, dummy]): " + str(result))

    def testChromosomeConversionHG19(self):
        """Test that an hg19 build with chrom = 23 or 24 gets converted to X or Y
        """
        self.assertEqual(MutUtils.convertChromosomeStringToMutationDataFormat("23", build="hg19"), "X", "chrom of 23 did not produce X: " + MutUtils.convertChromosomeStringToMutationDataFormat("23", build="hg19"))
        self.assertEqual(MutUtils.convertChromosomeStringToMutationDataFormat("24", build="hg19"), "Y", "chrom of 24 did not produce Y: " + MutUtils.convertChromosomeStringToMutationDataFormat("24", build="hg19"))

        self.assertEqual(MutUtils.convertChromosomeStringToMutationDataFormat("2", build="hg19"), "2", "chrom of 2 yielded different value: " + MutUtils.convertChromosomeStringToMutationDataFormat("2", build="hg19"))
        self.assertEqual(MutUtils.convertChromosomeStringToMutationDataFormat("4", build="hg19"), "4", "chrom of 4 yielded different value: " + MutUtils.convertChromosomeStringToMutationDataFormat("4", build="hg19"))

    def testChrom2HashCodeTable(self):
        chroms = ["1", "X", "3", "contig1", "Y", "25", "mt"]
        h = MutUtils.createChrom2HashCodeTable(chroms)
        self.assertTrue(h["1"] == 1, "For chrom 1, hash code should be 1 but it was %s." % h["1"])
        self.assertTrue(h["3"] == 3, "For chrom 3, hash code should be 3 but it was %s." % h["3"])
        self.assertTrue(h["25"] == 25, "For chrom 25, hash code should be 25 but it was %s." % h["25"])
        self.assertTrue(h["X"] == 26, "For chrom X, hash code should be 26 but it was %s." % h["X"])
        self.assertTrue(h["Y"] == 27, "For chrom Y, hash code should be 27 but it was %s." % h["Y"])
        self.assertTrue(h["mt"] == 28, "For chrom mt, hash code should be 28 but it was %s." % h["mt"])
        self.assertTrue(h["contig1"] == 29, "For chrom contig1, hash code should be 29 but it was %s." % h["contig1"])

        chroms = ["contig1", "mt"]
        h = MutUtils.createChrom2HashCodeTable(chroms)
        self.assertTrue(h["mt"] == 3, "For chrom mt, hash code should be 3 but it was %s." % h["mt"])
        self.assertTrue(h["contig1"] == 4, "For chrom contig1, hash code should be 4 but it was %s." % h["contig1"])

    def testRetrievePrecedingBasesForDeletions(self):
        chrom = "1"
        start = 1234567
        end = 1234567  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "GTC"
        alt_allele = "G"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_ref_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForDeletions(mut)
        mut.ref_allele = updated_ref_allele
        mut.alt_allele = "-"
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        self.assertTrue("_preceding_bases" in mut, "_preceding_bases is missing in the mutation data.")
        self.assertTrue(mut.start == 1234568, "Mut start should be 1234568 but was %s." % mut.start)
        self.assertTrue(mut.end == 1234569, "Mut end should be 1234569 but was %s." % mut.end)
        self.assertTrue(mut.ref_allele == "TC", "Ref allele should be TC but was %s." % mut.ref_allele)
        self.assertTrue(mut.alt_allele == "-", "Alt allele should be - but was %s." % mut.alt_allele)

        chrom = "1"
        start = 1234567
        end = 1234567  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "GTCT"
        alt_allele = "GTC"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_ref_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForDeletions(mut)
        mut.ref_allele = updated_ref_allele
        mut.alt_allele = "-"
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        self.assertTrue("_preceding_bases" in mut, "_preceding_bases is missing in the mutation data.")
        self.assertTrue(mut.start == 1234570, "Mut start should be 1234570 but was %s." % mut.start)
        self.assertTrue(mut.end == 1234570, "Mut end should be 1234570 but was %s." % mut.end)
        self.assertTrue(mut.ref_allele == "T", "Ref allele should be T but was %s." % mut.ref_allele)
        self.assertTrue(mut.alt_allele == "-", "Alt allele should be - but was %s." % mut.alt_allele)

        chrom = "1"
        start = 152497145
        end = 152497145  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "CCCGAGCTGCTTACGATAGCCTTCTT"
        alt_allele = "C"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_ref_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForDeletions(mut)
        mut.ref_allele = updated_ref_allele
        mut.alt_allele = "-"
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        self.assertTrue("_preceding_bases" in mut, "_preceding_bases is missing in the mutation data.")
        self.assertTrue(mut.start == 152497146, "Mut start should be 1234570 but was %s." % mut.start)
        self.assertTrue(mut.end == 152497170, "Mut end should be 1234570 but was %s." % mut.end)
        self.assertTrue(mut.ref_allele == "CCGAGCTGCTTACGATAGCCTTCTT", "Ref allele should be T but was %s."
                                                                       % mut.ref_allele)
        self.assertTrue(mut.alt_allele == "-", "Alt allele should be - but was %s." % mut.alt_allele)

    def testRetrievePrecedingBasesForInsertions(self):
        chrom = "1"
        start = 1234567
        end = 1234567  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "GTC"
        alt_allele = "GTCT"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_alt_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForInsertions(mut)
        mut.ref_allele = "-"
        mut.alt_allele = updated_alt_allele
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        self.assertTrue("_preceding_bases" in mut, "_preceding_bases is missing in the mutation data.")
        self.assertTrue(mut.start == 1234570, "Mut start should be 1234570 but was %s." % mut.start)
        self.assertTrue(mut.end == 1234570, "Mut end should be 1234570 but was %s." % mut.end)
        self.assertTrue(mut.ref_allele == "-", "Ref allele should be - but was %s." % mut.ref_allele)
        self.assertTrue(mut.alt_allele == "T", "Alt allele should be T but was %s." % mut.alt_allele)

        chrom = "1"
        start = 1234567
        end = 1234567  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "GTC"
        alt_allele = "GTCTT"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_alt_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForInsertions(mut)
        mut.ref_allele = "-"
        mut.alt_allele = updated_alt_allele
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        self.assertTrue("_preceding_bases" in mut, "_preceding_bases is missing in the mutation data.")
        self.assertTrue(mut.start == 1234570, "Mut start should be 1234570 but was %s." % mut.start)
        self.assertTrue(mut.end == 1234571, "Mut end should be 1234571 but was %s." % mut.end)
        self.assertTrue(mut.ref_allele == "-", "Ref allele should be - but was %s." % mut.ref_allele)
        self.assertTrue(mut.alt_allele == "TT", "Alt allele should be TT but was %s." % mut.alt_allele)

    def testRetrievePrecedingBaseFromAnnotationForDeletions(self):
        chrom = "1"
        start = 1234568
        end = 1234569  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "GTC"
        alt_allele = "G"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_ref_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForDeletions(mut)
        mut.ref_allele = updated_ref_allele
        mut.alt_allele = "-"
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        updated_ref_allele, updated_alt_allele, updated_start = \
            MutUtils.retrievePrecedingBaseFromAnnotationForDeletions(mut)
        self.assertTrue(updated_start == start, "Mut start should be %s but was %s." % (start, updated_start))
        self.assertTrue(updated_ref_allele == ref_allele, "Ref allele should be %s but was %s."
                                                          % (ref_allele, updated_ref_allele))
        self.assertTrue(updated_alt_allele == alt_allele, "Alt allele should be %s but was %s." % (alt_allele, updated_alt_allele))

        chrom = "1"
        start = 1234567
        end = 1234567  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "GTCT"
        alt_allele = "GTC"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_ref_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForDeletions(mut)
        mut.ref_allele = updated_ref_allele
        mut.alt_allele = "-"
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        updated_ref_allele, updated_alt_allele, updated_start = \
            MutUtils.retrievePrecedingBaseFromAnnotationForDeletions(mut)
        self.assertTrue(updated_start == start, "Mut start should be %s but was %s." % (start, updated_start))
        self.assertTrue(updated_ref_allele == ref_allele, "Ref allele should be %s but was %s."
                                                          % (ref_allele, updated_ref_allele))
        self.assertTrue(updated_alt_allele == alt_allele, "Alt allele should be %s but was %s."
                                                          % (alt_allele, updated_alt_allele))

        chrom = "1"
        start = 152497145
        end = 152497145  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "CCCGAGCTGCTTACGATAGCCTTCTT"
        alt_allele = "C"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_ref_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForDeletions(mut)
        mut.ref_allele = updated_ref_allele
        mut.alt_allele = "-"
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        updated_ref_allele, updated_alt_allele, updated_start = \
            MutUtils.retrievePrecedingBaseFromAnnotationForDeletions(mut)
        self.assertTrue(updated_start == start, "Mut start should be %s but was %s." % (start, updated_start))
        self.assertTrue(updated_ref_allele == ref_allele, "Ref allele should be %s but was %s."
                                                          % (ref_allele, updated_ref_allele))
        self.assertTrue(updated_alt_allele == alt_allele, "Alt allele should be %s but was %s."
                                                          % (alt_allele, updated_alt_allele))

    def testRetrievePrecedingBaseFromAnnotationForInsertions(self):
        chrom = "1"
        start = 1234567
        end = 1234567  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "GTC"
        alt_allele = "GTCT"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_alt_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForInsertions(mut)
        mut.ref_allele = "-"
        mut.alt_allele = updated_alt_allele
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        updated_ref_allele, updated_alt_allele, updated_start = \
            MutUtils.retrievePrecedingBaseFromAnnotationForInsertions(mut)
        self.assertTrue(updated_start == start, "Mut start should be %s but was %s." % (start, updated_start))
        self.assertTrue(updated_ref_allele == ref_allele, "Ref allele should be %s but was %s."
                                                          % (ref_allele, updated_ref_allele))
        self.assertTrue(updated_alt_allele == alt_allele, "Alt allele should be %s but was %s."
                                                          % (alt_allele, updated_alt_allele))

        chrom = "1"
        start = 1234567
        end = 1234567  # incorrect, but doesn't matter for the purposed of testing
        ref_allele = "GTC"
        alt_allele = "GTCTT"
        build = "19"
        mut = MutationData(chrom, start, end, ref_allele, alt_allele, build)
        preceding_bases, updated_alt_allele, updated_start, updated_end = \
            MutUtils.retrievePrecedingBasesForInsertions(mut)
        mut.ref_allele = "-"
        mut.alt_allele = updated_alt_allele
        mut.start = updated_start
        mut.end = updated_end
        mut.createAnnotation(annotationName=MutUtils.PRECEDING_BASES_ANNOTATION_NAME, annotationValue=preceding_bases)
        updated_ref_allele, updated_alt_allele, updated_start = \
            MutUtils.retrievePrecedingBaseFromAnnotationForInsertions(mut)
        self.assertTrue(updated_start == start, "Mut start should be %s but was %s." % (start, updated_start))
        self.assertTrue(updated_ref_allele == ref_allele, "Ref allele should be %s but was %s."
                                                          % (ref_allele, updated_ref_allele))
        self.assertTrue(updated_alt_allele == alt_allele, "Alt allele should be %s but was %s."
                                                          % (alt_allele, updated_alt_allele))

if __name__ == '__main__':
    unittest.main()
