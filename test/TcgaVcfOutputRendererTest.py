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

import logging
import os
import unittest

from oncotator.Annotator import Annotator
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MutationData import MutationData
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.output.TcgaVcfOutputRenderer import TcgaVcfOutputRenderer
from TestUtils import TestUtils
from oncotator.utils.GenericTsvReader import GenericTsvReader
import vcf


TestUtils.setupLogging(__file__, __name__)
class TcgaVcfOutputRendererTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()

    def testHeaderCreation(self):
        """Test that a tcga vcf header can be generated, even from a blank mutation. """
        vcfOR = TcgaVcfOutputRenderer("out/TCGAVCFHeader.out.txt")
        m = MutationData()
        m.createAnnotation('center', "broad.mit.edu")
        hdr = vcfOR.createVcfHeader(m)
        self.assertTrue(hdr is not None)
        self.assertTrue(hdr <> "")
        self.assertTrue(hdr.find("broad.mit.edu") <> -1, "Could not find string that should have been in header.")

    def testChromRendering(self):
        """Make sure that the chromosome rendering in TCGA VCF is correct: "1" --> "1", "GLXXXX.Y" --> GLXXXX.Y, not <GLXXXX.Y>"""
        vcfOR = TcgaVcfOutputRenderer("out/TCGAVCF.empty.out.txt")
        testChrs = ["21", "MT", "GL1234.4", "1"]
        gt = ["21", "MT", "GL1234.4", "1"]
        ctr = 0
        for t in testChrs:
            val = vcfOR._renderChrom(t)
            self.assertTrue(val == gt[ctr], "Chrom value did not match ground truth: " + t + " --> " + val + "  GT: " + gt[ctr])
            ctr += 1

    def _createDatasourcesForTesting(self):
        dbDir = self.config.get('DEFAULT',"dbDir")
        return DatasourceFactory.createDatasources(dbDir, "hg19",isMulticore=False)

    def _createManualAnnotations(self):
        # These should be passed in to the oncotator via "-a"
        result = {"build":"37", 'center':"broad.mit.edu", 'individual_barcode':"TCGA-individual1",
                  'normal_accession':"accessionN", 'normal_barcode':"TCGA-ind1-N", 'normal_file':".", 'normal_uuid':"uuidN",
                  'platform':"Illumina", 'softwareName':"", 'softwareParams':"", 'softwareVersion':"", 'source':"",
                  'tumor_accession':"accessionT", 'tumor_barcode':"TCGA-ind1-T", 'tumor_file':".", 'tumor_subtype':"Primary",
                  'tumor_uuid':"uuidT", 'vcfProcessLog':"<InputVCF=<.>,InputVCFSource=<.>,InputVCFVer=<.>,InputVCFParam=<.>,InputVCFgeneAnno=<https://tcga-data.nci.nih.gov/docs/GAF/GAF3.0/>>",
                  'geneAnno':"https://tcga-data.nci.nih.gov/docs/GAF/GAF3.0/"}
        return result

    def testFullSnpVcf(self):
        """ Perform test of a SNP call stats (maflite) all the way through TCGA VCF creation.  Only checks that a file was created.
        """
        outputFilename = "out/TCGAVCFTest.snp.vcf"
        callStatsIn = MafliteInputMutationCreator("testdata/Test.call_stats.trim.txt")
        vcfOR = TcgaVcfOutputRenderer(outputFilename)
        datasources = self._createDatasourcesForTesting()

        annotator = Annotator()
        annotator.setInputCreator(callStatsIn)
        annotator.setOutputRenderer(vcfOR)
        annotator.setManualAnnotations(self._createManualAnnotations())
        for ds in datasources:
            annotator.addDatasource(ds)
        annotator.annotate()

        self.assertTrue(os.path.exists(outputFilename))

    def testFullIndelVcf(self):
        """ Perform test of a Indel maflite all the way through TCGA VCF creation
        """
        outputFilename = "out/TCGAVCFTest.indel.vcf"
        callStatsIn = MafliteInputMutationCreator("testdata/maflite/Patient0.indel.maf.txt")
        vcfOR = TcgaVcfOutputRenderer(outputFilename)
        datasources = self._createDatasourcesForTesting()

        annotator = Annotator()
        annotator.setInputCreator(callStatsIn)
        annotator.setOutputRenderer(vcfOR)
        annotator.setManualAnnotations(self._createManualAnnotations())
        for ds in datasources:
            annotator.addDatasource(ds)
        annotator.annotate()

        self.assertTrue(os.path.exists(outputFilename))

        # Check that the deletions have position decremented by one from what is present in the maflite
        #  Checking that 1	36643701 in the maflite (a deletion) becomes 1	36643700 in the vcf, but that the others are
        #  the same.
        maflite_ic = MafliteInputMutationCreator("testdata/maflite/Patient0.indel.maf.txt")
        muts = maflite_ic.createMutations()
        vcf_reader = vcf.Reader(open(outputFilename, 'r'))

        vcf_pos = [int(rec.POS) for rec in vcf_reader]
        for m in muts:
            # If the variant is a deletion, then the vcf position should be the same as maflite minus one.  Otherwise, the same.
            is_variant_deletion = (m.alt_allele == "") or (m.alt_allele == "-") or (m.alt_allele == ".")
            if is_variant_deletion:
                self.assertTrue((int(m.start) - 1) in vcf_pos, "Deletion was not correct for " + m.chr + ":" + m.start)
            else:
                self.assertTrue(int(m.start) in vcf_pos, "Insertion was not correct for " + m.chr + ":" + m.start)

    def _testInfoField(self, filter):
        outputFilename = "out/TCGAVCFTest.indel.vcf.dummy"
        vcfOR = TcgaVcfOutputRenderer(outputFilename)
        mq0 = "0"
        ss = "Somatic"
        m = MutationData()
        m.createAnnotation('t_ref_count', '20')
        m.createAnnotation('t_alt_count', '25')
        m.createAnnotation('n_ref_count', '100')
        m.createAnnotation('n_alt_count', '150')
        m.createAnnotation('dbSNP_RS', '')
        m.createAnnotation('gene', 'FAKE')
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', 'Missense')
        m.createAnnotation('transcript_id', 'tid001')
        infoData = vcfOR._generateInfoField(m, filter, mq0, ss)
        return infoData

    def testPassInfoFieldGeneration(self):
        """Test simple info field generation for pass"""
        filter='PASS'
        infoData = self._testInfoField(filter)
        self.assertIsNotNone(infoData)
        self.assertTrue(infoData <> "")
        self.assertTrue(infoData.find("SOMATIC") <> -1, "SOMATIC not found")
        self.assertTrue(infoData.find("Gene=FAKE") <> -1, "Gene not found")

    def testFailInfoFieldGeneration(self):
        """Test simple info field generation for fail"""
        filter='mf1'
        infoData = self._testInfoField(filter)
        self.assertIsNotNone(infoData)
        self.assertTrue(infoData <> "")
        self.assertTrue(infoData.find("SOMATIC") <> -1, "SOMATIC not found")
        self.assertTrue(infoData.find("Gene") == -1, "Gene was found when it should have been missing.")

    def testAnotherFullSNP(self):
        """Test SNP call stats .  Just make sure no exception is thrown."""
        inputFile = "testdata/maflite/Another.call_stats.txt"
        outputFilename = "out/Another.call_stats.out.vcf"
        callStatsIn = MafliteInputMutationCreator(inputFile)
        vcfOR = TcgaVcfOutputRenderer(outputFilename)
        datasources = self._createDatasourcesForTesting()

        annotator = Annotator()
        annotator.setInputCreator(callStatsIn)
        annotator.setOutputRenderer(vcfOR)
        annotator.setManualAnnotations(self._createManualAnnotations())
        for ds in datasources:
            annotator.addDatasource(ds)
        annotator.annotate()

        self.assertTrue(os.path.exists(outputFilename))
        statinfo = os.stat(outputFilename)
        self.assertTrue(statinfo.st_size > 0, "Generated VCF file (" + outputFilename + ") is empty.")

    def testPopulatedButNullValuesInInitNLod(self):
        """Test that if init_n_lod is "." or "", there is no error """
        m = MutationData()
        m.createAnnotation("init_n_lod", "")
        outputFilename = "out/blank.vcf"
        vcfOR = TcgaVcfOutputRenderer(outputFilename)
        lod = vcfOR._extract_lod(m,"init_n_lod")
        self.assertEqual(lod, 50)

        m["init_n_lod"] = '.'
        lod = vcfOR._extract_lod(m, "init_n_lod")
        self.assertEqual(lod, 50)

        m["init_n_lod"] = '6'
        lod = vcfOR._extract_lod(m, "init_n_lod")
        self.assertEqual(lod, 6)

        m["init_n_lod"] = '6.8'
        lod = vcfOR._extract_lod(m, "init_n_lod")
        self.assertEqual(lod, 6)

        m["init_n_lod"] = '-12.8'
        lod = vcfOR._extract_lod(m, "init_n_lod")
        self.assertEqual(lod, -12)

        m.createAnnotation("t_lod_fstar", "")
        lod = vcfOR._extract_lod(m, "t_lod_fstar")
        self.assertEqual(lod, 50)

        m["t_lod_fstar"] = '.'
        lod = vcfOR._extract_lod(m, "t_lod_fstar")
        self.assertEqual(lod, 50)

        m["t_lod_fstar"] = '6'
        lod = vcfOR._extract_lod(m, "t_lod_fstar")
        self.assertEqual(lod, 6)

        m["t_lod_fstar"] = '6.8'
        lod = vcfOR._extract_lod(m, "t_lod_fstar")
        self.assertEqual(lod, 6)

        m["t_lod_fstar"] = '-12.8'
        lod = vcfOR._extract_lod(m, "t_lod_fstar")
        self.assertEqual(lod, -12)

    def testEmptyInput(self):
        """Make sure that we can generate an empty vcf from an empty maflite"""
        inputFile = "testdata/maflite/empty.maflite"
        outputFilename = "out/empty.vcf"
        callStatsIn = MafliteInputMutationCreator(inputFile)
        vcfOR = TcgaVcfOutputRenderer(outputFilename)
        datasources = self._createDatasourcesForTesting()

        annotator = Annotator()
        annotator.setInputCreator(callStatsIn)
        annotator.setOutputRenderer(vcfOR)
        annotator.setManualAnnotations(self._createManualAnnotations())
        for ds in datasources:
            annotator.addDatasource(ds)
        annotator.annotate()

        self.assertTrue(os.path.exists(outputFilename))
        statinfo = os.stat(outputFilename)
        self.assertTrue(statinfo.st_size > 0, "Generated VCF file (" + outputFilename + ") is empty.")

    def testMafInput(self):
        """Make sure that we can render a TCGA VCF from a TCGA MAF -- using no datasources"""
        inputFile = "testdata/maf/Patient1.snp.maf.annotated"
        outputFilename = "out/maf2tcgavcf.vcf"
        mafIn = MafliteInputMutationCreator(inputFile)
        vcfOR = TcgaVcfOutputRenderer(outputFilename)

        annotator = Annotator()
        annotator.setInputCreator(mafIn)
        annotator.setOutputRenderer(vcfOR)
        annotator.setManualAnnotations(self._createManualAnnotations())
        annotator.annotate()
        self.assertTrue(os.path.exists(outputFilename))
        statinfo = os.stat(outputFilename)
        self.assertTrue(statinfo.st_size > 0, "Generated VCF file (" + outputFilename + ") is empty.")

if __name__ == '__main__':
    unittest.main()
