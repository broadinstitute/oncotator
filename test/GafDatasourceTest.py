# LICENSE_GOES_HERE
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
from oncotator.datasources.Gaf import Gaf
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.utils.MultiprocessingUtils import LoggingPool
from oncotator.utils.MutUtils import MutUtils
from TestUtils  import TestUtils
from oncotator.utils.VariantClassification import VariantClassification


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

globalConfig = TestUtils.createUnitTestConfig()

TestUtils.setupLogging(__file__, __name__)
# globalConfig = TestUtils.createUnitTestConfig()
# @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test.  GAF 3.0 will not be supported for much longer.")
@unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default datasource corpus, with GAF 3.0, is needed to run this test.")
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)

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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['gene'] == "AKT1", "Incorrect gene found: " + m['gene'] + "  If updating GAF, this may not be an error, but should be confirmed manually.")

    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def test_effect_tx_mode(self):
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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


    @unittest.skipIf(not os.path.exists(globalConfig.get("gaf3.0", "gafDir")), "Default Datasource, with GAF 3.0, corpus is needed to run this test")
    def testSpliceSiteWithinNBases(self):
        """Test that a silent mutation is changed to splice site w/in 10 bases of a splice site """
        # chr21:10,998,326-10,998,346
        # 10,998,336 is a splice site.  (Junction between 10998335 and 336)
        # AGTTCTCCTT C TGGAAAAAAG
        refs = 'AGTTCTCCTTCTGGAAAAAAG'
        alts = 'TCAGACTGAAAATACCCCCCT'
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)

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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)
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
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)

        gafDatasource.set_tx_mode(TranscriptProvider.TX_MODE_BEST_EFFECT)
        old_hashcode = gafDatasource.get_hashcode()
        gafDatasource.set_tx_mode(TranscriptProvider.TX_MODE_CANONICAL)
        new_hashcode = gafDatasource.get_hashcode()
        self.assertTrue(old_hashcode != new_hashcode)

    def test_start_codon(self):
        """Test a start codon hit in a GAF datasource"""
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)

        m = MutationData()
        m.start = str(22221729)
        m.end = str(22221729)
        m.chr="22"
        m.ref_allele = 'A'
        m.alt_allele = 'T'
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['variant_classification'] == VariantClassification.MISSENSE)

    @unittest.skip("GAF 3.0 datasources are not being supported much longer, but this test may have exposed a minor bug, so is being preserved if a bugfix is implemented.")
    def test_denovo(self):
        """GAF de novo test """
        gafDatasource = TestUtils.createTranscriptProviderDatasource(self.config)

        m = MutationData()
        m.start = str(22221735)
        m.end = str(22221737)
        m.chr="22"
        m.ref_allele = ''
        m.alt_allele = 'CAT'
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['variant_classification'] == 'De_novo_Start_OutOfFrame')

        m = MutationData()
        m.start = str(22221735)
        m.end = str(22221740)
        m.chr="22"
        m.ref_allele = ''
        m.alt_allele = 'AACATAA'
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['variant_classification'] == 'De_novo_Start_OutOfFrame')

        m = MutationData()
        m.start = str(22221735)
        m.end = str(22221739)
        m.chr="22"
        m.ref_allele = ''
        m.alt_allele = 'ACATAA'
        m = gafDatasource.annotate_mutation(m)
        self.assertTrue(m['variant_classification'] == 'De_novo_Start_InFrame')


    def _flattenChunks(self, chunks):
        [[(yield m) for m in c] for c in chunks]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()














