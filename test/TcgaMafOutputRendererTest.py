# LICENSE_GOES_HERE
from TestUtils import TestUtils
from oncotator.input.VcfInputMutationCreator import VcfInputMutationCreator
from oncotator.utils.OptionConstants import OptionConstants


"""
Created on Nov 8, 2012

@author: lichtens
"""
import unittest

import logging 
from oncotator.utils.ConfigUtils import ConfigUtils
from oncotator.MutationData import MutationData
from oncotator.Annotator import Annotator
from oncotator.output.TcgaMafOutputRenderer import TcgaMafOutputRenderer
from oncotator.utils.OncotatorCLIUtils import OncotatorCLIUtils
import os
from oncotator.utils.GenericTsvReader import GenericTsvReader
from oncotator.DatasourceFactory import DatasourceFactory

TestUtils.setupLogging(__file__, __name__)


class TcgaMafOutputRendererTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    """ These are the default overrides for generating a TCGA MAF file.  These will appear on all mutations, but are here for a test.
        These were taken from version 0.5.25.0 of Oncotator.
    """
    TCGA_MAF_DEFAULTS = {'NCBI_Build': '37', 'Strand': "+", 'Center': 'broad.mit.edu', 'source': 'Capture',
                         'status': 'Somatic', 'phase': 'Phase_I', 'sequencer': 'Illumina GAIIx',
                         'Tumor_Validation_Allele1': '', 'Tumor_Validation_Allele2': '',
                         'Match_Norm_Validation_Allele1': '', 'Match_Norm_Validation_Allele2': '',
                         'Verification_Status': '', 'Validation_Status': '', 'Validation_Method': '', 'Score': '',
                         'BAM_file': '', 'Match_Norm_Seq_Allele1': '', 'Match_Norm_Seq_Allele2': ''}

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def _validateTcgaMafContents(self, filename):
        """ This is a utility, private method for unit tests to get a semblance that a valid maf file was created.  
        
        Note: This method has nothing to do with the TCGA validator.
        
        """
        configFile = ConfigUtils.createConfigParser(os.path.join("configs", "tcgaMAF2.4_output.config"))
        statinfo = os.stat(filename)
        self.assertTrue(statinfo.st_size > 0, "Generated MAF file (" + filename + ") is empty.")
        
        tsvReader = GenericTsvReader(filename)
        
        self.assertTrue(tsvReader.getComments().find('#version') != -1, "First line did not specify a version number")

        ctr = 1
        for lineDict in tsvReader:

            # TODO: Re-enable when GENCODE and HGNC datasources are concordant (or Entrez_Gene_ID is in the gencode gtf)
            # if lineDict['Entrez_Gene_Id'] == "0":
            #     self.assertTrue(lineDict['Hugo_Symbol'] == "Unknown", "Entrez_Gene_Id was zero, but Hugo Symbol was not 'Unknown'.  Line: " + str(ctr))
            
            unknownKeys = []
            self.assertTrue(lineDict["Tumor_Seq_Allele1"] != lineDict["Tumor_Seq_Allele2"], "Reference and alternate were equal in TCGA MAF output on line %d (%s)" % (ctr, lineDict["Tumor_Seq_Allele1"]))
            for k in lineDict.keys():
                if lineDict[k] == "__UNKNOWN__":
                    unknownKeys.append(k)

                self.assertTrue('\r' not in lineDict[k], "Carriage return character found in an annotation value.")

                requiredColumns = configFile.get("general", "requiredColumns")
                optionalColumns = configFile.get("general", "optionalColumns")
                exposedColumns = configFile.get("general", "exposedColumns")
                if (k not in requiredColumns) and (k not in optionalColumns) and (k not in exposedColumns):
                    self.assertTrue(k.startswith("i_"), "Internal column was not prepended with 'i_'")
                
            unknownKeys.sort()
            self.assertTrue(len(unknownKeys) == 0, "__UNKNOWN__ values (" + str(len(unknownKeys)) + ") seen on line " + str(ctr) + ", in fields: " + ", ".join(unknownKeys))
            
            ctr += 1

    def _determine_db_dir(self):
        return self.config.get('DEFAULT',"dbDir")

    def _annotateTest(self, inputFilename, outputFilename, datasource_dir, inputFormat="MAFLITE", outputFormat="TCGAMAF", default_annotations=TCGA_MAF_DEFAULTS, override_annotations=None, is_skip_no_alts=False):
        self.logger.info("Initializing Annotator...")

        if override_annotations is None:
            override_annotations = dict()

        annotator = Annotator()
        runSpec = OncotatorCLIUtils.create_run_spec(inputFormat, outputFormat, inputFilename, outputFilename, defaultAnnotations=default_annotations, datasourceDir=datasource_dir, globalAnnotations=override_annotations, is_skip_no_alts=is_skip_no_alts)
        annotator.initialize(runSpec)
        self.logger.info("Annotation starting...")
        return annotator.annotate()
    
    def testFullSNPOutput(self):
        """ Create a TCGA MAF from a SNP TSV file."""
        self.logger.info("Initializing Maflite SNP Test...")
        
        testOutputFilename = self._annotateTest('testdata/maflite/Patient0.snp.maf.txt', "out/testSNP_v2.4.maf.tsv", self._determine_db_dir())

        # Sanity checks to make sure that the generated maf file is not junk.
        self._validateTcgaMafContents(testOutputFilename)

    def testFullIndelOutput(self):
        """ Create a TCGA MAF from an Indel TSV file."""
        self.logger.info("Initializing Maflite indel Test...")
        
        testOutputFilename = self._annotateTest('testdata/maflite/Patient0.indel.maf.txt', "out/testIndel_v2.4.maf.tsv", self._determine_db_dir())
        
        # Sanity checks to make sure that the generated maf file is not junk.
        self._validateTcgaMafContents(testOutputFilename)

    def testEmptyInput(self):
        """ Create a TCGA MAF from an empty maflite file."""

        testOutputFilename = self._annotateTest('testdata/maflite/empty.maflite', "out/empty.maf.tsv", self._determine_db_dir())

        # Sanity checks to make sure that the generated maf file is not junk.
        self._validateTcgaMafContents(testOutputFilename)

    def testInternalFields(self):
        """ Test that an annotation that is not listed explicitly in the required or optional columns is rendered with i_ prepended """
        outputFilename = "out/testInternalFields_v2.4.maf.tsv"
        m = MutationData()
        m.createAnnotation("TEST", "THIS IS A TEST", "TESTING")
        
        # The next annotation is real and should not be considered internal.
        m.createAnnotation("gene", "EGFR")
        
        outputRenderer = TcgaMafOutputRenderer(outputFilename, configFile='configs/tcgaMAF2.4_output.config')
        outputRenderer.renderMutations(iter([m]), ['No comments'])
        
        configFile = ConfigUtils.createConfigParser('configs/tcgaMAF2.4_output.config')
        requiredColumns = configFile.get("general", "requiredColumns")
        self.assertTrue("Hugo_Symbol" in requiredColumns, " This test assumes that Hugo_Symbol is a required column in the TCGA MAF.  If not, the test must be modified.")

        statinfo = os.stat(outputFilename)
        self.assertTrue(statinfo.st_size > 0, "Generated MAF file (" + outputFilename + ") is empty.")
        
        tsvReader = GenericTsvReader(outputFilename)
        headers = tsvReader.getFieldNames()
        self.assertTrue("Hugo_Symbol" in headers, "Hugo_Symbol not found in output headers")
        self.assertTrue("TEST" not in headers, "TEST was found in output headers when it should have been renamed to i_TEST")
        self.assertTrue("i_TEST" in headers, "i_TEST not found in output headers")

    def testInternalFieldsSkipPrepend(self):
        """ Test that an annotation that is not listed explicitly in the required or optional columns is rendered with i_ prepended """
        outputFilename = "out/testInternalFields_v2.4.maf.tsv"
        m = MutationData()
        m.createAnnotation("TEST", "THIS IS A TEST", "TESTING")

        # The next annotation is real and should not be considered internal.
        m.createAnnotation("gene", "EGFR")

        outputRenderer = TcgaMafOutputRenderer(outputFilename, configFile='configs/tcgaMAF2.4_output.config', other_options={OptionConstants.NO_PREPEND:True})
        outputRenderer.renderMutations(iter([m]), ['No comments'])

        configFile = ConfigUtils.createConfigParser('configs/tcgaMAF2.4_output.config')
        requiredColumns = configFile.get("general", "requiredColumns")
        self.assertTrue("Hugo_Symbol" in requiredColumns, " This test assumes that Hugo_Symbol is a required column in the TCGA MAF.  If not, the test must be modified.")

        statinfo = os.stat(outputFilename)
        self.assertTrue(statinfo.st_size > 0, "Generated MAF file (" + outputFilename + ") is empty.")

        tsvReader = GenericTsvReader(outputFilename)
        headers = tsvReader.getFieldNames()
        self.assertTrue("Hugo_Symbol" in headers, "Hugo_Symbol not found in output headers")
        self.assertTrue("i_TEST" not in headers, "i_TEST was found in output headers when prepend was disabled.")
        self.assertTrue("TEST" in headers, "TEST was not found in output headers.")


    def testMutationDatasources(self):
        """ Test that we can create a simple TSV output from all of the current datasources (specified in the config file).  Note that no validation is done.  Simply that the output file was created. 
        TODO: This unit test needs to be moved."""
        testOutputFilename = self._annotateTest('testdata/maflite/Patient0.snp.maf.txt', "out/testsimpleSNP.maf.tsv", self._determine_db_dir(), outputFormat="SIMPLE_TSV")
        self.assertTrue(os.path.exists(testOutputFilename))

    def testExposedColumns(self):
        """Test that columns listed in the config file as exposed do not get the i_ prepend"""
        testOutputFilename = self._annotateTest('testdata/maflite/tiny_maflite.maf.txt', "out/testExposedCols.maf.tsv", self._determine_db_dir())

        # Sanity checks to make sure that the generated maf file is not junk.
        self._validateTcgaMafContents(testOutputFilename)

        # Check the columns, since the input has a couple of exposed columns.
        tsvReader = GenericTsvReader(testOutputFilename)
        headers = tsvReader.getFieldNames()
        headersToCheck = ['t_alt_count', 't_ref_count']
        for h in headersToCheck:
            self.assertFalse(("i_" + h) in headers, "i_ was prepended to " + h)
            self.assertTrue(h in headers, h + " not found.")

    def testProperConversionVcfToMaf(self):
        """Test that ref, alt, and positions are properly populated in a TCGA MAF generated from a VCF """

        # For this conversion, you must specify the barcodes manually
        override_annotations = TcgaMafOutputRendererTest.TCGA_MAF_DEFAULTS
        override_annotations.update({'tumor_barcode':'Patient0-Tumor', 'normal_barcode':'Patient0-Normal'})

        outputFilename = self._annotateTest('testdata/vcf/Patient0.somatic.strelka.indels.vcf', "out/testConversionFromVCF.maf.annotated", self._determine_db_dir(), inputFormat="VCF", outputFormat="TCGAMAF", override_annotations=override_annotations, is_skip_no_alts=True)

        # Sanity checks to make sure that the generated maf file is not junk.
        self._validateTcgaMafContents(outputFilename)

        # Check to make sure that the ref and alt are correct for a TCGA MAF.
        tsvReader = GenericTsvReader(outputFilename)

        ctr = 0

        for line_dict in tsvReader:
            ref = line_dict['Reference_Allele']
            alt = line_dict['Tumor_Seq_Allele2']

            # INS
            if len(alt) > len(ref):
                self.assertTrue(ref == "-", "Invalid insertion with " + ref + "  " + alt)

            # DEL
            if len(ref) > len(alt):
                self.assertTrue(alt == "-", "Invalid deletion with " + ref + "  " + alt)

            self.assertTrue(line_dict['Start_position'] in ["10089935", "57493929", "155301009", "64948169", "64948166",
                                                            "64948167", "64948168"])
            self.assertTrue(line_dict['Reference_Allele'] in ["-", "TC", "A", "TT", "TTT"])
            self.assertTrue(line_dict['Tumor_Seq_Allele2'] in ["-", "TC", "G", "T"])
            self.assertTrue(line_dict['Matched_Norm_Sample_Barcode'] == "Patient0-Normal")
            self.assertTrue(line_dict['Matched_Norm_Sample_UUID'] == "Patient0-Normal")
            self.assertTrue(line_dict['Tumor_Sample_Barcode'] == "Patient0-Tumor")
            self.assertTrue(line_dict['Tumor_Sample_UUID'] == "Patient0-Tumor")
            ctr += 1

        self.assertTrue(ctr == 8, str(ctr) + " mutations found, but should have been 8.")

    def testProperConversionVcfToMafWithThirdSample(self):
        """Test that ref, alt, and positions are properly populated in a TCGA MAF generated from a VCF, but that the NORMAL is treated as any other sample, since this VCF has three samples in it. """

        # For this conversion, you must specify the barcodes manually
        override_annotations = TcgaMafOutputRendererTest.TCGA_MAF_DEFAULTS
        override_annotations.update({'tumor_barcode': 'NA'})

        outputFilename = self._annotateTest(os.path.join(*["testdata", "vcf", "Patient0.somatic.strelka.indels.met.vcf"]),
                                            os.path.join("out", "testConversionFromVCFv2.maf.annotated"),
                                            self._determine_db_dir(), inputFormat="VCF", outputFormat="TCGAMAF",
                                            override_annotations=override_annotations)

        # Sanity checks to make sure that the generated maf file is not junk.
        self._validateTcgaMafContents(outputFilename)

        # Check to make sure that the ref and alt are correct for a TCGA MAF.
        tsvReader = GenericTsvReader(outputFilename)

        ctr = 0

        for line_dict in tsvReader:
            ctr += 1

        self.assertTrue(ctr == 24, str(ctr) + " mutations found, but should have been 24." )

    def test_validation_correction(self):
        """ Test that the validation allele fields are determined automatically when not specified by the user for invalid mutation.
        """
        m = MutationData()
        m.chr = "3"
        m.start = "178948145"
        m.end = "178948145"
        m.alt_allele = "A"
        m.ref_allele = "G"
        m['validation_status'] = "Invalid"
        m['Match_Norm_Validation_Allele1'] = ""
        m['Match_Norm_Validation_Allele2'] = ""
        m['Tumor_Validation_Allele1'] = ""
        m['Tumor_Validation_Allele2'] = ""
        m['Mutation_Status'] = "Somatic"

        output_filename = os.path.join("out", "test_validation_correction1.maf.tsv")

        outputRenderer = TcgaMafOutputRenderer(output_filename,
                                               configFile=os.path.join("configs", "tcgaMAF2.4_output.config"))
        outputRenderer.renderMutations([m].__iter__())

        tsv_reader = GenericTsvReader(output_filename)

        for line_dict in tsv_reader:
            self.assertTrue(line_dict['Match_Norm_Validation_Allele1'] == line_dict['Match_Norm_Validation_Allele2'], "Matched norm alleles did not match.")
            self.assertTrue(line_dict['Tumor_Validation_Allele1'] == line_dict['Tumor_Validation_Allele2'], "Tumor alleles did not match for an invalid validation result.")
            self.assertTrue(line_dict['Match_Norm_Validation_Allele1'] == line_dict['Tumor_Validation_Allele2'], "Tumor alleles did not match normal alleles for an invalid validation result.")
            self.assertTrue(line_dict['Match_Norm_Validation_Allele1'] == line_dict['Reference_Allele'], "Norm validation alleles did not match reference (norm, reference): (%s, %s)" %(line_dict['Match_Norm_Validation_Allele1'] ,line_dict['Reference_Allele']) )
            self.assertTrue("G" == line_dict['Reference_Allele'], "Reference allele should have been G, but was " + line_dict['Reference_Allele'])
            self.assertTrue("None" == line_dict['Mutation_Status'], "Mutation Status must be None when Validation Status is Invalid: " + line_dict['Mutation_Status'])

    def test_validation_correction_valid(self):
        """ Test that the validation allele fields are determined automatically when not specified by the user for a valid mutation.
        """
        m = MutationData()
        m.chr = "3"
        m.start = "178948145"
        m.end = "178948145"
        m.alt_allele = "A"
        m.ref_allele = "G"
        m['validation_status'] = "Valid"
        m['Match_Norm_Validation_Allele1'] = ""
        m['Match_Norm_Validation_Allele2'] = ""
        m['Tumor_Validation_Allele1'] = ""
        m['Tumor_Validation_Allele2'] = ""
        m['Mutation_Status'] = "Somatic"

        output_filename = os.path.join("out", "test_validation_correction2.maf.tsv")

        outputRenderer = TcgaMafOutputRenderer(output_filename,
                                               configFile=os.path.join("configs", "tcgaMAF2.4_output.config"))
        outputRenderer.renderMutations([m].__iter__())

        tsv_reader = GenericTsvReader(output_filename)

        for line_dict in tsv_reader:
            self.assertTrue(line_dict['Match_Norm_Validation_Allele1'] == line_dict['Match_Norm_Validation_Allele2'], "Matched norm alleles did not match.")
            self.assertTrue(line_dict['Tumor_Validation_Allele1'] == line_dict['Reference_Allele'], "Tumor validation allele 1 did not match reference for a valid validation result.")
            self.assertTrue(line_dict['Tumor_Validation_Allele2'] == line_dict['Tumor_Seq_Allele2'], "Tumor validation allele 2 did not match Tumor_Seq_Allele2 for a valid validation result.")
            self.assertTrue(line_dict['Match_Norm_Validation_Allele1'] == line_dict['Tumor_Validation_Allele1'], "Tumor allele 1 did not match normal alleles for a valid validation result.")
            self.assertTrue(line_dict['Match_Norm_Validation_Allele1'] == line_dict['Reference_Allele'], "Norm validation alleles did not match reference (norm, reference): (%s, %s)" %(line_dict['Match_Norm_Validation_Allele1'] ,line_dict['Reference_Allele']) )
            self.assertTrue("G" == line_dict['Reference_Allele'], "Reference allele should have been G, but was " + line_dict['Reference_Allele'])
            self.assertTrue("A" == line_dict['Tumor_Seq_Allele2'], "Alt allele should have been A, but was " + line_dict['Tumor_Seq_Allele2'])


if __name__ == "__main__":
    unittest.main()