# LICENSE_GOES_HERE
import unittest
import logging
import os

import pandas
import vcf

from oncotator.utils.MutUtils import MutUtils
from oncotator.input.VcfInputMutationCreator import VcfInputMutationCreator
from oncotator.Annotator import Annotator
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from TestUtils import TestUtils
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.utils.ConfigUtils import ConfigUtils
from oncotator.output.TcgaMafOutputRenderer import TcgaMafOutputRenderer
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.utils.TagConstants import TagConstants


TestUtils.setupLogging(__file__, __name__)


class VcfInputMutationCreatorTest(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def _createGafDataSource(self):
        self.logger.info("Initializing gaf 3.0")
        return TestUtils.createTranscriptProviderDatasource(self.config)

    def testBasicCreationWithExampleVcf(self):
        """
        Tests the ability to parse an input VCF file can be parsed without any errors.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])

        creator = VcfInputMutationCreator(inputFilename)
        muts = creator.createMutations()

        # You cannot use len(muts), since muts is a generator.
        ctr = 0
        for m in muts:
            ctr += 1
        self.assertTrue(ctr == 27, "Should have seen 27 (# REF alleles x # samples) mutations, but saw: " + str(ctr))
        self.assertTrue((m.chr == "21") and (m.start == 1234569), "Last mutation was not correct: " + str(m))

        # Reminder: muts is a generator, so it has to be reset
        creator.reset()
        muts = creator.createMutations()
        ctr = 0
        for m in muts:
            ctr += 1
        self.assertTrue(ctr == 27, "Should have seen 27 called mutations, but saw: " + str(ctr))

    def testSimpleAnnotationWithExampleVcf(self):
        """
        Tests the ability to do a simple Gaf 3.0 annotation.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])
        outputFilename = os.path.join("out", "simpleVCF.Gaf.annotated.out.tsv")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename, [])
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createTranscriptProviderDatasource(self.config))
        annotator.annotate()

    def testSimpleAnnotationWithAComplexVcf(self):
        """
        Tests the ability to parse a rather complex VCF file without any errors.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "random.vcf"])
        outputFilename = os.path.join("out", "random.tsv")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename, [])
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

    def _createTCGAMAFOverridesForVCF(self):
        """
        These are the default overrides for generating a TCGA MAF file.  These will appear on all mutations, but are here for a test.
        These were taken from version 0.5.25.0 of Oncotator.
        """
        #TODO: Remove the 'Match_Norm_Seq_Allele1' and 'Match_Norm_Seq_Allele2' from this list and populate properly, if possible.
        result = dict(source='Capture', status='Somatic', phase='Phase_I', sequencer='Illumina GAIIx',
                      Tumor_Validation_Allele1='', Tumor_Validation_Allele2='', Match_Norm_Validation_Allele1='',
                      Match_Norm_Validation_Allele2='', Verification_Status='', Validation_Status='',
                      Validation_Method='', Score='', BAM_file='', Match_Norm_Seq_Allele1='', Match_Norm_Seq_Allele2='',
                      Tumor_Sample_UUID='', Tumor_Sample_Barcode='', Strand="+", Center="broad.mit.edu",
                      NCBI_Build="37")
        return result

    def _createDatasourceCorpus(self):
        dbDir = self.config.get('DEFAULT', "dbDir")
        return DatasourceFactory.createDatasources(dbDir, "hg19", isMulticore=False)

    def testTCGAMAFRendering(self):
        """
        Tests the ability to render a germline VCF file as a TCGA MAF file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])
        outputFilename = os.path.join("out", "example.vcf.maf.annotated")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = TcgaMafOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.setManualAnnotations(self._createTCGAMAFOverridesForVCF())
        datasources = self._createDatasourceCorpus()
        for ds in datasources:
            annotator.addDatasource(ds)
        filename = annotator.annotate()

        self._validateTcgaMafContents(filename)

    def _validateTcgaMafContents(self, filename):
        """
        This is a utility, private method for unit tests to get a semblance that a valid maf file was created.
        
        Note: This method has nothing to do with the TCGA validator.
        
        TODO: This is code duplication from TCGA MAF Output RendererTest.  This should be refactored into a base class
        (to preserve self.assertTrue, etc).
        """
        statinfo = os.stat(filename)
        self.assertTrue(statinfo.st_size > 0, "Generated MAF file (" + filename + ") is empty.")

        tsvReader = GenericTsvReader(filename)

        self.assertTrue(tsvReader.getComments().find('#version') <> -1, "First line did not specify a version number")

        ctr = 1
        for lineDict in tsvReader:
            if lineDict['Entrez_Gene_Id'] == "0":
                self.assertTrue(lineDict['Hugo_Symbol'] == "Unknown",
                                "Entrez_Gene_Id was zero, but Hugo Symbol was not 'Unknown'.  Line: " + str(ctr))

            unknownKeys = []
            for k in lineDict.keys():
                if lineDict[k] == "__UNKNOWN__":
                    unknownKeys.append(k)

                self.assertTrue('\r' not in lineDict[k], "Carriage return character found in an annotation value.")

                configFile = ConfigUtils.createConfigParser('configs/tcgaMAF2.3_output.config')
                requiredColumns = configFile.get("general", "requiredColumns")
                optionalColumns = configFile.get("general", "optionalColumns")
                if (k not in requiredColumns) and (k not in optionalColumns):
                    self.assertTrue(k.startswith("i_"), "Internal column was not prepended with 'i_'")

            unknownKeys.sort()
            self.assertTrue(len(unknownKeys) == 0,
                            "__UNKNOWN__ values (" + str(len(unknownKeys)) + ") seen on line " + str(
                                ctr) + ", in fields: " + ", ".join(unknownKeys))

            ctr += 1

    def testSwitchedFieldsWithExampleVcf(self):
        """
        Tests whether the switched tags are ignored.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.bad.switched.fields.vcf"])
        outputFilename = os.path.join("out", "example.out.tsv")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename, [])
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)

    def testAnnotationWithExampleVcf(self):
        """
        Tests whether parsed annotations match the actual annotations.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])
        outputFilename = os.path.join("out", "example.out.tsv")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "example.expected.out.tsv"])

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        tsvReader = GenericTsvReader(outputFilename)

        current = pandas.read_csv(outputFilename, sep='\t', header=len(tsvReader.getCommentsAsList()))
        expected = pandas.read_csv(expectedOutputFilename, sep='\t')

        currentColNames = set()
        for i in range(len(current.columns)):
            currentColNames.add(current.columns[i])

        expectedColNames = set()
        for i in range(len(expected.columns)):
            expectedColNames.add(expected.columns[i])

        self.assertTrue(len(currentColNames.symmetric_difference(expectedColNames)) is 0,
                        "Should have the same columns")
        self.assertTrue(len(current.index) == len(expected.index), "Should have the same number of rows")

        for colName in currentColNames:
            self.assertTrue(sum((current[colName] == expected[colName]) | (pandas.isnull(current[colName]) &
                                                                           pandas.isnull(expected[colName]))) ==
                            len(current.index), "Should have the same values in column " + colName)

    def testAnnotationWithNoSampleNameExampleVcf(self):
        """
        Tests whether parsed annotations match the actual annotations when the input is a VCF file that has no samples.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.sampleName.removed.vcf"])
        outputFilename = os.path.join("out", "example.sampleName.removed.out.tsv")

        creator = VcfInputMutationCreator(inputFilename)
        renderer = SimpleOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

    def testGetMetaDataWithNoSampleNameExampleVcf(self):
        """
        Tests to ensure that the metadata can be retrieved even before createMutations has been called.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.sampleName.removed.vcf"])

        creator = VcfInputMutationCreator(inputFilename)
        gtKeys = {'genotype', 'read_depth', 'genotype_quality', 'haplotype_quality', 'q10', 's50', 'samples_number',
                  'depth_across_samples', 'allele_frequency', 'ancestral_allele', 'dbSNP_membership', 'id', 'qual',
                  'hapmap2_membership'}
        md = creator.getMetadata()
        ks = set(md.keys())
        diff = gtKeys.symmetric_difference(ks)
        self.assertTrue(len(diff) == 0, "Missing keys that should have been seen in the metadata: " + str(diff))

    def testSNPsAndIndelStartAndEndPos(self):
        """
        Tests that the start and end positions of SNPs and Indels are parsed as defined by the NCI's MAF specification
        (https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification).
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.snps.indels.vcf"])
        outputFilename = os.path.join("out", "example.snps.indels.out.tsv")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        tsvReader = GenericTsvReader(outputFilename)
        for row in tsvReader:
            if row['start'] == "16890445":
                self.assertEqual(row["end"], "16890445", "The value should be %s but it was %s." % ("16890445",
                                                                                                    row["end"]))
            elif row["start"] == "154524458":
                self.assertEqual(row["end"], "154524459", "The value should be %s but it was %s." % ("154524459",
                                                                                                     row["end"]))
            elif row["start"] == "114189432":
                self.assertEqual(row["end"], "114189433", "The value should be %s but it was %s." % ("114189433",
                                                                                                     row["end"]))

    def testSplitByNumberOfAltsWithFile(self):
        """
        Tests whether we properly determine that a field is split using an actual file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.split.tags.vcf"])
        creator = VcfInputMutationCreator(inputFilename)
        isSplit = dict()
        isSplit['read_depth'] = False
        isSplit['ESP_MAF'] = False
        isSplit['allele_frequency'] = True

        mapVcfFields2Tsv = dict()
        mapVcfFields2Tsv['read_depth'] = 'DP'
        mapVcfFields2Tsv['ESP_MAF'] = 'ESP_MAF'
        mapVcfFields2Tsv['allele_frequency'] = 'AF'

        muts = creator.createMutations()

        vcfReader = vcf.Reader(filename=inputFilename, strict_whitespace=True)

        chrom = None
        pos = None
        variant = None
        for m in muts:
            if (chrom != m['chr']) or (pos != m['start']):
                chrom = m['chr']
                pos = m['start']
                variant = vcfReader.next()

            for annotationName in isSplit.keys():
                if mapVcfFields2Tsv[annotationName] in variant.INFO:
                    a = m.getAnnotation(annotationName)
                    self.assertTrue((TagConstants.SPLIT in a.getTags()) == isSplit[annotationName],
                                    "Is " + annotationName + " split for chrom " + chrom + ", pos " + str(pos) +
                                    "? " + str(isSplit[annotationName]) + ", but saw: " +
                                    str(TagConstants.SPLIT in a.getTags()))

    def testGenotypeFieldIsHonored(self):
        """
        Tests that no issues arise with genotype values >1 when multiple variants appear on one line.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.severalGTs.vcf"])
        creator = VcfInputMutationCreator(inputFilename)
        muts = creator.createMutations()
        ctr = 0
        for mut in muts:

            if MutUtils.str2bool(mut["alt_allele_seen"]):
                self.assertTrue(mut['sample_name'] != "NA 00001")
                ctr += 1
        self.assertTrue(ctr == 7,
                        str(ctr) + " mutations with alt seen, but expected 7.  './.' should not show as a variant.")

    def testDuplicateAnnotation(self):
        """
        Tests that the duplicate annotations are parsed correctly.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.duplicate_annotation.vcf"])
        outputFilename = os.path.join("out", "example.duplicate_annotation.out.tsv")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        tsvReader = GenericTsvReader(outputFilename)
        fieldnames = tsvReader.getFieldNames()
        self.assertTrue("SS" in fieldnames, "SS field is missing in the header.")
        self.assertTrue("SS__FORMAT__" in fieldnames, "SS__FORMAT__ is missing in the header.")

        row = tsvReader.next()
        self.assertTrue("SS" in row, "SS field is missing in the row.")
        self.assertTrue("SS__FORMAT__" in row, "SS__FORMAT__ is missing in the row.")

        self.assertEqual("2", row["SS"], "Incorrect value of SS.")
        self.assertEqual("0", row["SS__FORMAT__"], "Incorrect value of SS__FORMAT__")

    def testDuplicateAnnotationMetaData(self):
        """
        Tests that the metadata is populated correctly in cases where duplicate annotations are present in the input VCF
        file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.duplicate_annotation.vcf"])

        creator = VcfInputMutationCreator(inputFilename)
        md = creator.getMetadata()

        self.assertTrue("SS" in md, "SS field is missing in metadata.")
        self.assertTrue("SS__FORMAT__" in md, "SS__FORMAT__ is missing in metadata.")

    def testMissingFilter(self):
        """
        Tests that the missing FILTER fields are parsed correctly.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.missing_filters.vcf"])
        outputFilename = os.path.join("out", "example.missing_filters.out.tsv")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "example.expected.missing_filters.out.tsv"])

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        tsvReader = GenericTsvReader(outputFilename)

        current = pandas.read_csv(outputFilename, sep='\t', header=len(tsvReader.getCommentsAsList()))
        expected = pandas.read_csv(expectedOutputFilename, sep='\t')

        currentColNames = set()
        for i in range(len(current.columns)):
            currentColNames.add(current.columns[i])

        expectedColNames = set()
        for i in range(len(expected.columns)):
            expectedColNames.add(expected.columns[i])

        self.assertTrue(len(currentColNames.symmetric_difference(expectedColNames)) is 0,
                        "Should have the same columns")
        self.assertTrue(len(current.index) == len(expected.index), "Should have the same number of rows")

        for colName in currentColNames:
            self.assertTrue(sum((current[colName] == expected[colName]) | (pandas.isnull(current[colName]) &
                                                                           pandas.isnull(expected[colName]))) ==
                            len(current.index), "Should have the same values in column " + colName)


if __name__ == "__main__":
    unittest.main()