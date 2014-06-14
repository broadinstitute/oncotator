# LICENSE_GOES_HERE

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
        """Make sure that the chromosome rendering in TCGA VCF is correct: "1" --> "1"  ,  "GLXXXX.Y" --> <GLXXXX.Y>"""
        vcfOR = TcgaVcfOutputRenderer("out/TCGAVCF.empty.out.txt")
        testChrs = ["21", "MT", "GL1234.4", "1"]
        gt = ["21", "MT", "<GL1234.4>", "1"]
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
