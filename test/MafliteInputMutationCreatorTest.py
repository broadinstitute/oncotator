# LICENSE_GOES_HERE



"""
Created on Nov 9, 2012

@author: lichtens
"""
import unittest

from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.input.MafliteMissingRequiredHeaderException import MafliteMissingRequiredHeaderException
from oncotator.output.TcgaMafOutputRenderer import TcgaMafOutputRenderer
from oncotator.utils.MutUtils import MutUtils
from oncotator.Annotator import Annotator
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.utils.GenericTsvReader import GenericTsvReader
import os
from TestUtils import TestUtils
TestUtils.setupLogging(__file__, __name__)

class MafliteInputMutationCreatorTest(unittest.TestCase):
    _multiprocess_can_split_ = True
    def setUp(self):
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def testMissingRequiredHeaders(self):
        try: 
            tmp = MafliteInputMutationCreator("testdata/maflite/brokenMaflite.tsv", 'configs/maflite_input.config')
            self.assertFalse(True, " Exception was not thrown")
        except MafliteMissingRequiredHeaderException as e:
            #str(e).find('alt_allele,end,ref_allele')<> -1
            missingCols = ['alt_allele', 'end', 'ref_allele']
            isMissingInData = []
            for c in missingCols:
                isMissingInData.append(str(e).find(c)<>-1)
            self.assertTrue(all(isMissingInData),"Incorrect columns identified as missing in maflite check: " + str(e))
            
    def testSimpleRead(self):
        """ Read a good maflite file and make sure that each mutation validates """
        tmp = MafliteInputMutationCreator("testdata/maflite/Patient0.indel.maf.txt", 'configs/maflite_input.config')
        muts = tmp.createMutations()
        
        # If no exception is thrown, then this test passes.
        for m in muts:
            MutUtils.validateMutation(m)
    
    def testNumberOfMuts(self):
        """ Make sure that the proper number of mutations were generated """
        inputFilename = "testdata/maflite/Patient0.snp.maf.txt"
        tmp = MafliteInputMutationCreator(inputFilename, 'configs/maflite_input.config')
        muts = tmp.createMutations()
        numMutsInput = len(file(inputFilename,'r').readlines()) - 1
        ctr = 0
        for m in muts:
            ctr += 1
        self.assertEqual(ctr, numMutsInput, "Did not see the proper number of mutations.")
    
    def testChromosomeM(self):
        """ Make sure that the chromosome created as M, rather than MT."""
        tmp = MafliteInputMutationCreator("testdata/maflite/chrM.maf.txt", 'configs/maflite_input.config')
        muts = tmp.createMutations()
        for m in muts:
            self.assertTrue(m.chr=="M", "mitochondria chromosome should be M, not " + m.chr)
            
    def testTCGAMAFAsInput(self):
        """ Test that we can take in a TCGA MAF (using MAFLITE), do no annotations, and still render it properly """
        tmp = MafliteInputMutationCreator("testdata/maf/Patient0.maf.annotated", 'configs/maflite_input.config')
        muts = tmp.createMutations()
        
        outputFilename = "out/testTCGAMAFAsInput.tsv"
        outputRenderer = TcgaMafOutputRenderer(outputFilename, 'configs/tcgaMAF2.4_output.config')
        outputRenderer.renderMutations(muts, tmp.getComments())
        
    def testTCGAMAFAsInputAndQuickAnnotate(self):
        """ Test that we can take in a TCGA MAF (using MAFLITE), do annotating, and still render it properly """
        inputFilename = "testdata/maf/Patient0.maf.annotated"
        tmp = MafliteInputMutationCreator(inputFilename, 'configs/maflite_input.config')
        outputFilename = "out/testTCGAMAFAsInputAndQuickAnnotate.tsv"
        outputRenderer = TcgaMafOutputRenderer(outputFilename, 'configs/tcgaMAF2.4_output.config')
        annotator = Annotator()
        
        annotator.setInputCreator(tmp)
        annotator.setOutputRenderer(outputRenderer)
        ds = DatasourceFactory.createDatasource("testdata/thaga_janakari_gene_ds/hg19/tj_data.config", "testdata/thaga_janakari_gene_ds/hg19/")
        annotator.addDatasource(ds)
        annotator.annotate()
        
        statinfo = os.stat(outputFilename)
        self.assertTrue(statinfo.st_size > 0, "Generated MAF file (" + outputFilename + ") is empty.")
        tsvReaderIn = GenericTsvReader(inputFilename)
        tsvReader = GenericTsvReader(outputFilename)
        
        self.assertTrue(tsvReader.getComments().find('#version') != -1, "First line did not specify a version number")
        self.assertTrue("i_TJ_Data_Why" in tsvReader.getFieldNames(), "New field missing (i_TJ_Data_Why) from header")
        self.assertTrue("i_TJ_Data_Who" in tsvReader.getFieldNames(), "New field missing (i_TJ_Data_Who) from header")
        
        ctrOut = 0
        for lineDict in tsvReader:
            ctrOut += 1
        ctrIn = 0
        for lineDict in tsvReaderIn:
            ctrIn += 1
        ctrIn += len(tsvReaderIn.getCommentsAsList())
        ctrOut += len(tsvReader.getCommentsAsList())

        self.assertTrue(ctrOut == (ctrIn + 2), "Output file should have same number of lines plus two (for maf version and Oncotator version comments) as input file.  (In,Out): " + str(ctrIn) + ", " + str(ctrOut))

    def testGetMetadata(self):
        """Make sure that we can retrieve metadata, even before createMutations has been called"""
        ic = MafliteInputMutationCreator("testdata/maflite/tiny_maflite.maf.txt")
        gtKeys = {'build', 'chr', 'start', 'end', 'ref_allele', 'alt_allele', 'tumor_barcode', 'normal_barcode',
                  'tumor_f', 'init_t_lod', 't_lod_fstar', 't_alt_count', 't_ref_count', 'judgement'}
        md = ic.getMetadata()
        ks = set(md.keys())
        diff = gtKeys.symmetric_difference(ks)
        self.assertTrue(len(diff) == 0, "Missing keys that should have been seen in the metadata: " + str(diff))

    def test_alt1_vs_alt2(self):
        """Test that we pick up the alternate that is different from the reference when both are specified"""
        ic = MafliteInputMutationCreator("testdata/maflite/alt1_vs_alt2.maflite")
        muts = ic.createMutations()
        ctr = 0
        for m in muts:
            ctr += 1
            self.assertTrue(m.alt_allele == "C", "Did not properly populate the alternate allele in line " + str(ctr) + "  " + m.alt_allele)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()