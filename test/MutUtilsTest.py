# LICENSE_GOES_HERE

import unittest

from oncotator.MutationData import MutationData
from oncotator.utils.MutUtils import MutUtils
from TestUtils import TestUtils
import vcf

TestUtils.setupLogging(__file__, __name__)


class MutUtilsTest(unittest.TestCase):

    _multiprocess_can_split_ = True

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
        self.assertTrue(mut.start == 1234569, "Mut start should be 1234570 but was %s." % mut.start)
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
        self.assertTrue(mut.start == 1234569, "Mut start should be 1234570 but was %s." % mut.start)
        self.assertTrue(mut.end == 1234570, "Mut end should be 1234571 but was %s." % mut.end)
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
