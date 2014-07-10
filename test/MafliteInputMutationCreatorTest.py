"""
By downloading the PROGRAM you agree to the following terms of use:

BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY

This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").

WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and

WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.

NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:

1. DEFINITIONS
1.1 "PROGRAM" shall mean the object code and source code known as Oncotator 1.0 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.  BROAD acknowledges that the PROGRAM employs one or more public domain code(s) that are freely available for public use.

2. LICENSE
2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.  LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.

3. OWNERSHIP OF INTELLECTUAL PROPERTY
LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.

Copyright 2014 Broad Institute, Inc.
Notice of attribution:  The Oncotator 1.0 program was made available through the generosity of the Broad Institute, Inc.

LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.

4. INDEMNIFICATION
LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorney fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.

5. NO REPRESENTATIONS OR WARRANTIES
THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.

6. ASSIGNMENT
This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.

7. MISCELLANEOUS
7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
"""
import shutil
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer


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

    def test_simple_seg_file_input(self):
        """Test that we can read in a seg file, do no annotation, and output as SIMPLE_TSV"""
        inputFilename = "testdata/seg/Patient0.seg.txt"
        output_filename = "out/test_simple_seg_file.tsv"
        if os.path.exists(output_filename):
            os.remove(output_filename)
        ic = MafliteInputMutationCreator(inputFilename, 'configs/seg_file_input.config')
        segs = ic.createMutations()

        i = 1
        for i,seg in enumerate(segs):
            pass

        self.assertTrue((i+1) == 27, "Found %d segments when there should have been 27." % (i+1) )

        ic = MafliteInputMutationCreator(inputFilename, 'configs/seg_file_input.config')
        segs = ic.createMutations()


        outputRenderer = SimpleOutputRenderer(output_filename, '')
        outputRenderer.renderMutations(segs)

        # Now check the output
        output_reader = GenericTsvReader(output_filename)

        required_cols = ["Sample", "Num_Probes", "Segment_Mean"]
        headers = output_reader.getFieldNames()
        for rcol in required_cols:
            self.assertTrue(rcol in headers)

        for line_dict in output_reader:
            self.assertTrue(line_dict['start'] is not None)
            self.assertTrue(line_dict['start'].strip() != "")
            self.assertTrue(line_dict['end'] is not None)
            self.assertTrue(line_dict['end'].strip() != "")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()