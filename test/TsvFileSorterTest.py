# LICENSE_GOES_HERE

import os
import unittest
import hashlib
from oncotator.utils.TsvFileSorter import TsvFileSorter
from TestUtils import TestUtils
from oncotator.utils.CallbackException import CallbackException

TestUtils.setupLogging(__file__, __name__)


class TsvFileSorterTest(unittest.TestCase):

    def testSortFile(self):
        """
        Tests sorting a file on the filesystem.
        """
        inputFilename = os.path.join(*["testdata", "small_cosmic_gpp", "small_cosmic_gpp.tempForSorting.tsv"])
        outputFilename = os.path.join("out", "small_cosmic_gpp.tempForSorting.out.tsv")
        tsvFileSorter = TsvFileSorter(inputFilename)
        func = lambda val: ((val["Gene_name"]).lower(), int(val["startAA"]), int(val["endAA"]))
        tsvFileSorter.sortFile(outputFilename, func)

        self.assertTrue(os.path.exists(outputFilename), "No file was generated.")

    def testSortFileWithSpaces(self):
        """
        Tests sorting a file with spaces in the headers on the filesystem.
        """
        inputFilename = os.path.join(*["testdata", "small_cosmic_with_gp_and_gpp", "small_cosmic_trimmed_for_sorting.txt.tbi.byAA"])
        outputFilename = os.path.join("out", "small_cosmic_trimmed_for_sorting.txt.byAA.sorted.tsv")
        tsvFileSorter = TsvFileSorter(inputFilename)
        func = lambda val: ((val["Gene name"]).lower(), int(val["startAA"]), int(val["endAA"]))
        tsvFileSorter.sortFile(outputFilename, func)

        self.assertTrue(os.path.exists(outputFilename), "No file was generated.")

    def testSortMixedCaps(self):
        """
        Tests sorting a file with mixed capitalization in the reference column.
        """
        inputFilename = os.path.join(*["testdata", "sort_mixed_caps_tsv", "sort_mixed_caps.tsv"])
        outputFilename = os.path.join("out", "sort_mixed_caps.tsv.sorted.out.tsv")
        tsvFileSorter = TsvFileSorter(inputFilename)
        func = lambda val: ((val["Gene name"]).lower(), int(val["startAA"]), int(val["endAA"]))
        tsvFileSorter.sortFile(outputFilename, func)

        self.assertTrue(os.path.exists(outputFilename), "No file was generated.")

        guessmd5 = hashlib.md5(file(outputFilename, 'r').read()).hexdigest()
        gtmd5 = hashlib.md5(file(os.path.join(*["testdata", "sort_mixed_caps_tsv", "sort_mixed_caps_sorted.tsv"]),
                                 "r").read()).hexdigest()
        self.assertTrue(guessmd5 == gtmd5)

    def testMultiplePartitionSorting(self):
        """
        Tests that the sorting works when the partition size is small and input file must be broken into multiple
        partitions.
        """
        inputFilename = os.path.join(*["testdata", "sort_mixed_caps_tsv", "sort_mixed_caps.tsv"])
        outputFilename = os.path.join("out", "multiple_partitions_sort_mixed_caps.tsv.sorted.out.tsv")
        tsvFileSorter = TsvFileSorter(inputFilename)
        func = lambda val: ((val["Gene name"]).lower(), int(val["startAA"]), int(val["endAA"]))
        tsvFileSorter.sortFile(outputFilename, func, 3)
        self.assertTrue(os.path.exists(outputFilename), "No file was generated.")

        guessmd5 = hashlib.md5(file(outputFilename, "r").read()).hexdigest()
        gtmd5 = hashlib.md5(file(os.path.join(*["testdata", "sort_mixed_caps_tsv", "sort_mixed_caps_sorted.tsv"]),
                                 "r").read()).hexdigest()
        self.assertTrue(guessmd5 == gtmd5)

    def testCallbackExceptionIncorrectType(self):
        """
        Tests that the CallbackException is raised when the input anonymous function does not return a tuple given a
        row.
        """
        inputFilename = os.path.join(*["testdata", "sort_mixed_caps_tsv", "sort_mixed_caps.tsv"])
        outputFilename = os.path.join("out", "multiple_partitions_sort_mixed_caps.tsv.sorted.out.tsv")
        tsvFileSorter = TsvFileSorter(inputFilename)
        func = lambda val: (val["Gene name"]).lower()
        try:
            tsvFileSorter.sortFile(outputFilename, func, 3)
        except CallbackException as msg:
            self.assertTrue(msg.value == "The value returned by the callback must be a tuple. Instead, a value of "
                                         "<type 'str'> was returned.", "Error msg is different.")

if __name__ == '__main__':
    unittest.main()
