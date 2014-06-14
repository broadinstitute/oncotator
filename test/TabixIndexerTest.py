# LICENSE_GOES_HERE
from TestUtils import TestUtils


'''
Created on Jan 10, 2013

@author: lichtens
'''
import shutil
import unittest
import os
import logging
from oncotator.index.TabixIndexer import TabixIndexer
import pysam
import vcf

TestUtils.setupLogging(__file__, __name__)


class TabixIndexerTest(unittest.TestCase):

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()

    def tearDown(self):
        pass
    
    def testBasicGeneProteinPositionIndexCreationMixedCaps(self):
        """
        Test that the mixed caps (or underscore) in gene name does not throw off the indexing.
        """
        inputFilename = "testdata/sort_mixed_caps_tsv/sort_mixed_caps.tsv"
        # Copy the inputFilename into output dir, since the indexing rewrites the input file
        copyFilename = "out/sort_mixed_caps.tsv"
        shutil.copy(inputFilename, copyFilename)

        fileName, fileExtension = os.path.splitext(copyFilename)
        outputFilename = fileName + ".sorted" + fileExtension
        resultIndexedFile = TabixIndexer.indexGeneProteinPosition("Gene name", "Mutation AA",
                                                                  copyFilename, outputFilename)

        self.assertTrue(os.path.exists(resultIndexedFile), "No index file was generated.")

    def testBasicGeneProteinPositionIndexCreation(self):
        """
        Creates a tabix index based on a small tsv file.
        """
        inFile = "testdata/small_cosmic_gpp/small_cosmic_gpp.tsv"
        outFile = "out/small_cosmic_gpp.out.tsv"
        TabixIndexer.indexGeneProteinPosition("Gene name", "Mutation AA", inFile, outFile)
        self.assertTrue(os.path.exists("out/small_cosmic_gpp.out.tabix_indexed.tsv.gz.tbi"), "Index file was not created.")
        tabixFile = pysam.Tabixfile("out/small_cosmic_gpp.out.tabix_indexed.tsv.gz")

        gene = "BRAF"
        startAA = 599
        endAA = 600
        results = tabixFile.fetch(reference=gene, start=startAA, end=endAA)
        ctr = 0
        try:
            while results.next():
                ctr += 1
        except StopIteration:
            pass

        self.assertTrue(ctr == 2, "Returned wrong number of results (gt: 2): " + str(ctr))

        gene = "CDKN2A"
        startAA = 29
        endAA = 30

        results = tabixFile.fetch(reference=gene, start=startAA, end=endAA)
        ctr = 0
        try:
            while results.next():
                ctr += 1
        except StopIteration:
            pass

        # Should have one result
        self.assertTrue(ctr == 1, "Returned wrong number of results (gt: 1): " + str(ctr))

    def testTabixIndexedVcfCreation(self):
        """
        Test the creation of VCF based tabix index file.
        """
        inFile = "testdata/vcf/example.vcf"
        destDir = "out"

        resultIndexedFile = TabixIndexer.index(destDir=destDir, inputFilename=inFile, preset="vcf")
        self.assertTrue(os.path.exists(resultIndexedFile), "No index file was generated.")

        vcfReader = vcf.Reader(filename=resultIndexedFile, compressed=True, strict_whitespace=True)
        vcfRecords = vcfReader.fetch(chrom=20, start=1230237, end=1230237)
        for vcfRecord in vcfRecords:
            self.assertEqual(vcfRecord.INFO["NS"], 3, "Expected %s but got %s." % (3, vcfRecord.INFO["NS"]))
            self.assertEqual(vcfRecord.INFO["DP"], 13, "Expected %s but got %s." % (13, vcfRecord.INFO["DP"]))

        os.remove(resultIndexedFile)

    def testTabixIndexedTsvCreation(self):
        inFile = "testdata/ESP6500SI-V2.chr1.snps_indels.head.25.txt"
        destDir = "out"

        # chr, startPos, endPos
        resultIndexedFile = TabixIndexer.index(destDir=destDir, inputFilename=inFile, fileColumnNumList=[0, 1, 1])
        self.assertTrue(os.path.exists(resultIndexedFile), "No index file was generated.")

        chrom = "1"
        start = "69594"
        end = "69594"
        tsvRecords = None
        tsvReader = pysam.Tabixfile(filename=resultIndexedFile)  # initialize the tsv reader
        try:
            tsvRecords = tsvReader.fetch(chrom, int(start)-1, int(end), parser=pysam.asTuple())
        except ValueError:
            pass

        tsvRecord = None
        for tsvRecord in tsvRecords:
            self.assertEqual(tsvRecord[5], "2,6190", "Value in column sixth does not match the expected value.")

        self.assertIsNotNone(tsvRecord, "No record for %s:%s-%s was found." % (chrom, start, end))

        os.remove(resultIndexedFile)

    def testExistingTabixIndexedFile(self):
        inFile = "testdata/example.vcf.gz"
        destDir = "out"

        resultIndexedFile = TabixIndexer.index(destDir=destDir, inputFilename=inFile)

        self.assertTrue(os.path.exists(resultIndexedFile), "No index file was generated.")

        vcfReader = vcf.Reader(filename=resultIndexedFile, compressed=True, strict_whitespace=True)
        vcfRecords = vcfReader.fetch(chrom=20, start=1230237, end=1230237)
        for vcfRecord in vcfRecords:
            self.assertEqual(vcfRecord.INFO["NS"], 3, "Expected %s but got %s." % (3, vcfRecord.INFO["NS"]))
            self.assertEqual(vcfRecord.INFO["DP"], 13, "Expected %s but got %s." % (13, vcfRecord.INFO["DP"]))

        os.remove(resultIndexedFile)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testBasicGeneIndexCreation']
    unittest.main()