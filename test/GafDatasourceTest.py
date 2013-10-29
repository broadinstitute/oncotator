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
from itertools import chain
import unittest
import os
import logging
import traceback
from cPickle import dump

from oncotator.Metadata import Metadata
from oncotator.MutationData import MutationData
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.datasources import Gaf, TranscriptProvider
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.utils.MultiprocessingUtils import LoggingPool
from oncotator.utils.MutUtils import MutUtils
from TestUtils  import TestUtils


def annotate_mutation_global(t):
    """Annotate from any datasource given a tuple that is (datasource, mutation).
    Mutation is a single mutation, not a list/generator."""
    ds = t[0]
    m = t[1]
    return ds.annotate_mutation(m)

def annotate_mutations_global(t):
    """Annotate from any datasource given a tuple that is (datasource, mutations).
    Mutations is a list."""
    ds = t[1]
    ms = t[0]
    result = []
    for m in ms:
        result.append(ds.annotate_mutation(m))
    return result

TestUtils.setupLogging(__file__, __name__)


class GafDatasourceTest(unittest.TestCase):

    # HACK: Allow config to be viewed by unittest decorators.
    globalConfig = TestUtils.createUnitTestConfig()

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()

    def tearDown(self):
        pass

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testNoUnknownAnnotations(self):
        """ Make sure that the gaf 3.0 datasource does not annotate anything with source set to Unknown """
        inputCreator = MafliteInputMutationCreator('testdata/maflite/Patient0.snp.maf.txt')
        gafDatasource = TestUtils.createGafDatasource(self.config)
        mutations = inputCreator.createMutations()    
        for m in mutations:
            m = gafDatasource.annotate_mutation(m)
            MutUtils.validateMutation(m)
            unknownAnnotations = MutUtils.getUnknownAnnotations(m)
            self.assertTrue(len(unknownAnnotations) == 0, "Unknown annotations exist in mutation: " + str(unknownAnnotations))

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testNoLostMutations(self):
        """ Does a simple gaf datasource annotation run and makes sure that no mutations were lost """
        inputFilename = 'testdata/maflite/Patient0.snp.maf.txt'
        inputCreator = MafliteInputMutationCreator(inputFilename, "configs/maflite_input.config")
        gafDatasource = TestUtils.createGafDatasource(self.config)

        numMutsInput = len(file(inputFilename, 'r').readlines()) - 1
        mutations = inputCreator.createMutations()  
        ctr = 0  
        for m in mutations:
            m = gafDatasource.annotate_mutation(m)
            MutUtils.validateMutation(m)
            ctr += 1
        self.assertEqual(ctr, numMutsInput, "Gaf data source altered mutation count.")


    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testChrM(self):
        """ Test that mitochondrial mutations can be annotated properly. """
        inputCreator = MafliteInputMutationCreator('testdata/maflite/chrM.maf.txt', "configs/maflite_input.config")
        gafDatasource = TestUtils.createGafDatasource(self.config)
        mutations = inputCreator.createMutations() 
        for m in mutations:
            try:
                m = gafDatasource.annotate_mutation(m)
                MutUtils.validateMutation(m)
            except Exception as e:
                # Fail this test because an exception was thrown
                self.assertTrue(False, "Erroneous exception was thrown: " + str(e) + "\n" + traceback.format_exc())
            self.assertTrue(m['gene'] != '')

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testChrGLs(self):
        """ Test that mutations on unaligned transcripts can be annotated properly.  I.e. when chromosome = GL....."""
        inputCreator = MafliteInputMutationCreator('testdata/maflite/chrGLs.maf.tsv', "configs/maflite_input.config")
        gafDatasource = TestUtils.createGafDatasource(self.config)
        mutations = inputCreator.createMutations() 
        for m in mutations:
            try:
                m = gafDatasource.annotate_mutation(m)
                MutUtils.validateMutation(m)
            except Exception as e:
                # Fail this test because an exception was thrown
                self.assertTrue(False, "Erroneous exception was thrown: " + str(e) + "\n" + traceback.format_exc())
            self.assertTrue(m['gene'] != '')

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testMC1R(self):
        """Test that this version of the GAF produces a MC1R, instead of TUBB gene"""
        m = MutationData()
        m.chr = '16'
        m.start = '89985913'
        m.end = '89985913'
        m.ref_allele = 'G'
        m.alt_allele = 'A'
        gafDatasource = TestUtils.createGafDatasource(self.config)
        m = gafDatasource.annotate_mutation(m)

        # At some point, we would expect this to be MC1R, not TUBB3
        self.assertTrue(m['gene'] == "TUBB3", "Incorrect gene found: " + m['gene'] + "  If updating GAF, this may not be an error, but should be confirmed manually.")


    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testAKT1(self):
        """ Test that this version of the GAF produces the up to date gene for a position given from a website user.
        """
        m = MutationData()
        m.chr = '14'
        m.start = '105246407'
        m.end = '105246407'
        m.ref_allele = 'G'
        m.alt_allele = 'A'
        gafDatasource = TestUtils.createGafDatasource(self.config)
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['gene'] == "AKT1", "Incorrect gene found: " + m['gene'] + "  If updating GAF, this may not be an error, but should be confirmed manually.")

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def test_effect_tx_mode(self):
        gafDatasource = TestUtils.createGafDatasource(self.config)
        gafDatasource.set_tx_mode(TranscriptProvider.TX_MODE_BEST_EFFECT)

        # Canonical mutation was Intron
        m = MutationData()
        m.chr = '2'
        m.start = '219137340'
        m.end = '219137340'
        m.ref_allele = 'G'
        m.alt_allele = 'T'
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['gene'] == "PNKD")
        self.assertTrue(m['variant_classification'] == "Missense_Mutation")

        gafDatasource.set_tx_mode(TranscriptProvider.TX_MODE_CANONICAL)
        m = MutationData()
        m.chr = '2'
        m.start = '219137340'
        m.end = '219137340'
        m.ref_allele = 'G'
        m.alt_allele = 'T'
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['gene'] == "PNKD")
        self.assertTrue(m['variant_classification'] == "Intron", "Canonical no longer is Intron.  This test is no longer valid.  This failure can come up when changing the GAF datasource.")


    @unittest.skip('This test has been disabled, since the backing implementation is incomplete.')
    def testBasicInitWithENSEMBL(self):
        gaf_fname = self.config.get("ENSEMBL", "gaf_fname")
        gaf_transcripts_fname = self.config.get("ENSEMBL", "gaf_transcript_seqs_fname")

        gafDatasource = Gaf(gaf_fname, gaf_transcripts_fname)

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testSpliceSiteWithinNBases(self):
        """Test that a silent mutation is changed to splice site w/in 10 bases of a splice site """
        # chr21:10,998,326-10,998,346
        # 10,998,336 is a splice site.  (Junction between 10998335 and 336)
        # AGTTCTCCTT C TGGAAAAAAG
        refs = 'AGTTCTCCTTCTGGAAAAAAG'
        alts = 'TCAGACTGAAAATACCCCCCT'
        gafDatasource = TestUtils.createGafDatasource(self.config)
        vcs = []
        for s in range(10998326, 10998347):
            m = MutationData()
            m.start = str(s)
            m.end = str(s)
            m.chr = "21"
            m.ref_allele = refs[s - 10998326]
            m.alt_allele = alts[s - 10998326]

            m = gafDatasource.annotate_mutation(m)

            distanceFromSpliceSite = abs(10998336 - int(m.start))
            vc = m['variant_classification']
            self.assertTrue(vc != 'Silent', 'Silent mutation found when it should be a splice site.')

            vcs.append(vc)
            print vc + "  " + m.start

        self.assertTrue(all([tmp == "Splice_Site" for tmp in vcs[8:12]]), "Not all vcs within 2 bases were splice site: " + str(vcs[8:12]))
        self.assertTrue(all([tmp != "Splice_Site" for tmp in vcs[0:8]]), "No splice sites should be seen: " + str(vcs[0:8]))
        self.assertTrue(all([tmp != "Splice_Site" for tmp in vcs[12:20]]), "No splice sites should be seen: " + str(vcs[12:20]))

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testSilentMutationGoingToSpliceSite(self):
        """Test that a silent mutation within 10 bp of a splice junction should become a splice site"""
        #chr1:28,233,780-28,233,805 Junction is at chr1:28,233,793 & 94
        #

        refs = "TGGGCTCGGGCTCTCTGAAAAGAAAA"
        alts = "TGGGCTCAGGCTCGCTGAAAAGAAAA"
        vcs = []
        gafDatasource = TestUtils.createGafDatasource(self.config)
        numSpliceSites = 0
        numSilent = 0
        startWindow = 28233780
        for s in range(startWindow, 28233806):
            m = MutationData()
            m.start = str(s)
            m.end = str(s)
            m.chr = "1"
            m.ref_allele = refs[s - startWindow]
            m.alt_allele = alts[s - startWindow]

            m = gafDatasource.annotate_mutation(m)

            distanceFromSpliceSite = abs(28233793 - int(m.start))
            vc = m['variant_classification']
            vcs.append(vc)
            # self.assertTrue(vc <> 'Silent', 'Silent mutation found when it should be a splice site.')

            if vc.lower() == "splice_site":
                numSpliceSites += 1
            if vc.lower() == "silent":
                numSilent += 1
            print vc + "  " + m.start + "  " + str(distanceFromSpliceSite)

        self.assertTrue(numSpliceSites == 4, "Should have seen 4 splice site mutations, but saw: " + str(numSpliceSites))
        self.assertTrue(numSilent == 11, "Should have seen 11 Silent mutations, but saw: " + str(numSilent))

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testFlank(self):
        """Test that we can see a Flank mutation."""
        #chr1:28,233,780-28,233,805 Junction is at chr1:28,233,793 & 94
        #

        refs = "TGGGCTCGGGCTCTCTGAAAAGAAAA"
        alts = "TGGGCTCAGGCTCTCTGAAAAGAAAA"
        vcs = []
        gafDatasource = TestUtils.createGafDatasource(self.config)
        numSpliceSites = 0
        numSilent = 0
        startWindow = 11042200
        for s in range(startWindow, startWindow+len(refs)):
            m = MutationData()
            m.start = str(s)
            m.end = str(s)
            m.chr="1"
            m.ref_allele = refs[s-startWindow]
            m.alt_allele = alts[s-startWindow]

            m = gafDatasource.annotate_mutation(m)

            vc = m['variant_classification']
            vcs.append(vc)

            print vc + "  " + m.start

        pass

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testFlank2(self):
        """Test a second real-world flank scenario"""
        gafDatasource = TestUtils.createGafDatasource(self.config)

        # 1	228646357 nearest Gene=HIST3H2A C>T
        m = MutationData()
        m.start = str(228646357)
        m.end = str(228646357)
        m.chr="1"
        m.ref_allele = 'C'
        m.alt_allele = 'T'
        m = gafDatasource.annotate_mutation(m)

        self.assertTrue(m['gene'] == "HIST3H2A", "Wrong gene (GT: HIST3H2A): " + m['gene'] + "   -- if updating GAF, this test may fail as this gene may not be appropriate.")
        self.assertTrue(m['variant_classification'] == "5'Flank", "Should be 5'Flank, but was " + m['variant_classification'] + " -- if updating GAF, this test may fail as this test is data specific.  Also, this may fail if padding parameters are changed.")

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testMicroRNA(self):
        """Test proper annotation of miRNA
        """
        #uc021qwk.1	chr12:31379258-31379277:-	hsa-miR-3194-3p|?	chr12:31379258-31379277:-		Confidence=100
        gafDatasource = TestUtils.createGafDatasource(self.config)
        m = MutationData()
        m.start = 31379268
        m.end = 31379268
        m.chr= "12"
        m.alt_allele = 'G'

        # This is accurate
        m.ref_allele = 'A'
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['gene'].lower() == "hsa-mir-3194-3p", "Wrong gene (GT: hsa-mir-3194-3p): " + m['gene'] + "   -- if updating GAF, this test may fail as this result may not be appropriate.")

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testExonRetrievalForGene(self):
        """Make sure that the GAF datasource can retrieve exons, given a gene"""
        testGeneList = ['CEBPA', 'BRCA1', 'PIK3CA']
        gafDatasource = TestUtils.createGafDatasource(self.config)
        for gene in testGeneList:
            exons = gafDatasource.retrieveExons(gene, isCodingOnly=True)
            self.assertTrue(exons is not None)
            print(str(exons))

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testExonRetrievalForGenesFromFile(self):
        """Make sure that the GAF datasource can retrieve exons, given a list of genes from a simple file"""
        inputGeneList = file('testdata/testGeneList.txt', 'r')
        outputFileFP = file("out/testGeneListExons.txt", 'w')
        errorFileFP = file("out/testGeneListExons.err.txt", 'w')
        gafDatasource = TestUtils.createGafDatasource(self.config)
        for line in inputGeneList:
            gene = line.strip()
            exons = gafDatasource.retrieveExons(gene, isCodingOnly=True)
            if len(exons) == 0:
                errorFileFP.write("Could not locate " + gene + "\n")
            for e in exons:
                outputFileFP.write('%s\t%s\t%s\t%s\n' % (e[0], e[1], e[2], e[3]))

    @unittest.skip("The backing code is experimental and should not be run.")
    # @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testMulticoreAnnotate(self):
        """Test a (too) simple annotating exercise from GAF on 2 cores"""
        gafDatasource = TestUtils.createGafDatasourceProxy(self.config)

        # Test pickling
        dump(gafDatasource, file('out/testGAFPickle.pkl','w'))

        m1 = MutationData()
        m1.chr = '3'
        m1.start = '178866811'
        m1.end = '178866811'
        m1.ref_allele = "A"
        m1.alt_allele = "C"
        m1.build = "hg19"

        m2 = MutationData()
        m2.chr = '3'
        m2.start = '178866812'
        m2.end = '178866812'
        m2.ref_allele = "A"
        m2.alt_allele = "C"
        m2.build = "hg19"

        p = LoggingPool(processes=2)
        result = p.map(annotate_mutation_global, [(gafDatasource, m1), (gafDatasource, m2)])
        p.close()
        p.join()

        for r in result:
            self.assertTrue("transcript_id" in r.keys())
            self.assertTrue("gene" in r.keys())
            self.assertTrue(r["gene"] == "PIK3CA")
        self.assertTrue(result[0].start != result[1].start)

    @unittest.skip("The backing code is experimental and should not be run.")
    # @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testMulticoreAnnotateFromChunkedFile(self):
        #TODO: Add unit test that Mutation data is pickle-able
        inputFile = "testdata/maflite/Patient0.snp.maf.txt"
        outputFile = "out/testGAFMulticorePatient0.snp.maf.txt"
        chunkSize = 200
        numChunks = 4


        gafDatasource = TestUtils.createGafDatasourceProxy(self.config)
        ic = MafliteInputMutationCreator(inputFile)
        oc = SimpleOutputRenderer(outputFile)

        # createChunks
        muts = ic.createMutations()

        allAnnotatedChunksFlat = []
        are_mutations_remaining = True
        p = LoggingPool(processes=numChunks)
        while are_mutations_remaining:

            chunks = []
            for j in xrange(0, numChunks):
                chunk = []
                for i in xrange(0, chunkSize):
                    try:
                        chunk.append(muts.next())
                    except StopIteration:
                        are_mutations_remaining = False
                        break

                chunks.append((chunk, gafDatasource))

            annotatedChunks = p.map(annotate_mutations_global, chunks)
            annotatedChunksFlat = self._flattenChunks(annotatedChunks)
            allAnnotatedChunksFlat.append(annotatedChunksFlat)
        p.close()
        p.join()

        annotatedMuts = chain.from_iterable(allAnnotatedChunksFlat)

        ctr = 0
        oc.renderMutations(annotatedMuts, Metadata())
        tsvReader = GenericTsvReader(outputFile)
        for line in tsvReader:
            ctr += 1
        self.assertTrue(ctr == 730, "Should have read 730 variants, but read " + str(ctr))

    def testChangeInTxModeChangesHashcode(self):
        """Test that a change in the tx-mode will change the hashcode"""
        gafDatasource = TestUtils.createGafDatasource(self.config)

        gafDatasource.set_tx_mode(TranscriptProvider.TX_MODE_BEST_EFFECT)
        old_hashcode = gafDatasource.get_hashcode()
        gafDatasource.set_tx_mode(TranscriptProvider.TX_MODE_CANONICAL)
        new_hashcode = gafDatasource.get_hashcode()
        self.assertTrue(old_hashcode != new_hashcode)

    def _flattenChunks(self, chunks):
        [[(yield m) for m in c] for c in chunks]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()














