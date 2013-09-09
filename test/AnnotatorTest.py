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
from TestUtils import TestUtils
from oncotator.DatasourceCreator import DatasourceCreator
from oncotator.MutationData import MutationData
from oncotator.utils.OncotatorCLIUtils import OncotatorCLIUtils, RunSpecification

"""
Created on Nov 7, 2012

@author: lichtens
"""
import unittest
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.Annotator import Annotator
import logging
import os
from oncotator.utils.GenericTsvReader import GenericTsvReader

TestUtils.setupLogging(__file__, __name__)
class AnnotatorTest(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def testBlankAnnotatorInit(self):
        """ Test an extremely simple scenario, where no additional annotations are needed.  I.e. no data sources """
        self.logger.info("Starting Blank Annotator Init Test...")

        inputCreator = MafliteInputMutationCreator('testdata/maflite/tiny_maflite.maf.txt')
        outputRenderer = SimpleOutputRenderer("out/testBlankAnnotatorTestFile.tsv")

        # Assumed myIC and myOC have been initialized as the proper Input and Output Creators, respectively.
        # 1) Initialize the Annotator
        annotator = Annotator()
        annotator.setInputCreator(inputCreator)
        annotator.setOutputRenderer(outputRenderer)
        testOutputFilename = annotator.annotate()

        # Test that file exists and that it has correct # of mutations (+1 for header +1 for annotator comment line).
        numSamples = 1
        numExtraLines = 3  # one for header, two for comment lines
        numDoubleLines = 0  # Number of lines with two alt alleles
        numVariants = 9
        gt = numSamples * numVariants + numDoubleLines * numSamples + numExtraLines
        fp = file(testOutputFilename, 'r')
        ctr = 0
        for line in fp:
            ctr += 1
        fp.close()
        self.assertEqual(ctr, gt,
                         "Number of lines read was not correct: " + str(ctr) + " -- should have been: " + str(gt))

    def testVersionHeader(self):
        """ This method simply tests that the version string returned by the annotator does not cause an exception.
            Minimal checking that the returned sting is actually correct.
            Does not attempt to initialize input or output.  Only a gaf datasource.
         """
        annotator = Annotator()
        annotator.addDatasource(TestUtils.createGafDatasource(self.config))
        tmp = annotator.createHeaderString()
        self.assertTrue(tmp.find("Gaf ") != -1, "Could not find Gaf version in header string.")
        self.assertTrue(tmp.find("Oncotator") != -1, "Could not find the word Oncotator in header string.")

    def testManualAnnotations(self):
        """ Test that the manual annotation facility in the Annotator is working properly. """
        annotator = Annotator()
        overrides = {'source': 'Capture', 'status': 'Somatic', 'phase': 'Phase_I', 'sequencer': 'Illumina GAIIx'}
        annotator.setManualAnnotations(overrides)
        inputCreator = MafliteInputMutationCreator('testdata/maflite/Patient0.snp.maf.txt')
        outputRenderer = SimpleOutputRenderer("out/testManualAnnotationsFile.tsv")
        annotator.setInputCreator(inputCreator)
        annotator.setOutputRenderer(outputRenderer)

        testOutputFilename = annotator.annotate()

        keysOfInterest = overrides.keys()

        statinfo = os.stat(testOutputFilename)
        self.assertTrue(statinfo.st_size > 0, "Generated TSV file (" + testOutputFilename + ") is empty.")

        tsvReader = GenericTsvReader(testOutputFilename)

        ctr = 1
        for lineDict in tsvReader:
            for k in keysOfInterest:
                self.assertTrue(lineDict[k] != "__UNKNOWN__",
                                "__UNKNOWN__ value seen on line " + str(ctr) + ", when it should be populated: " + k)
                self.assertTrue(lineDict[k] != "",
                                "Blank value seen on line " + str(ctr) + ", when it should be populated: " + k)
                self.assertTrue(lineDict[k] == overrides[k],
                                "Value for " + k + " on line " + str(ctr) + " did not match override: " + str(
                                    lineDict[k]) + " <> " + str(overrides[k]))
            ctr += 1

    def testDefaultAnnotations(self):
        """Test that the default annotation values populate properly. """
        annotator = Annotator()
        default_annotations = {"test2": "foo2", "test3": "Should not be seen"}
        overrides = {'test3': 'foo3'}

        m1 = MutationData()
        m1.createAnnotation("test1", "foo1")
        m1.createAnnotation("test2", "")

        m2 = MutationData()
        m2.createAnnotation("test1", "")


        m3 = MutationData()
        m3.createAnnotation("test1", "")
        m3.createAnnotation("test2", "foo2-original")

        muts = [m1, m2, m3]

        muts2 = annotator._applyManualAnnotations(muts, overrides)
        muts_final_gen = annotator._applyDefaultAnnotations(muts2, default_annotations)

        muts_final = []
        for m in muts_final_gen:
            self.assertTrue(m['test3'] == "foo3", "Override did not work")
            muts_final.append(m)

        self.assertTrue(muts_final[0]['test1'] == "foo1")
        self.assertTrue(muts_final[0]['test2'] == "foo2")
        self.assertTrue(muts_final[0]['test3'] == "foo3")

        self.assertTrue(muts_final[1]['test1'] == "")
        self.assertTrue(muts_final[1]['test2'] == "foo2")
        self.assertTrue(muts_final[1]['test3'] == "foo3")

        self.assertTrue(muts_final[2]['test1'] == "")
        self.assertTrue(muts_final[2]['test2'] == "foo2-original")
        self.assertTrue(muts_final[2]['test3'] == "foo3")

    def testAnnotateListOfMutations(self):
        """Test that we can initialize an Annotator, without an input or output and then feed mutations, one at a time... using a runspec"""
        dbDir = self.config.get('DEFAULT',"dbDir")
        ds = DatasourceCreator.createDatasources(dbDir)
        runSpec = RunSpecification()
        runSpec.initialize(None, None, datasources=ds)
        annotator = Annotator()
        annotator.initialize(runSpec)

        m = MutationData()
        m.chr = "1"
        m.start = "12941796"
        m.end = "12941796"
        m.alt_allele = "G"
        m.ref_allele = "T"

        muts = [m]

        muts = annotator.annotate_mutations(muts)
        m2 = muts.next()
        self.assertTrue(m2.get("gene", None) is not None)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testBasicAnnotatorInit']
    unittest.main()