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
from oncotator.DatasourceCreator import DatasourceCreator
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
        return TestUtils.createGafDatasource(self.config)

    def testBasicCreationWithExampleVcf(self):
        inputFilename = 'testdata/vcf/example.vcf'
        
        creator = VcfInputMutationCreator(inputFilename)
        muts = creator.createMutations()
        
        # You cannot use len(muts), since muts is a generator.
        ctr = 0
        for m in muts:
            ctr += 1
        self.assertTrue(ctr == 27, "Should have seen 27 (# REF alleles x # samples) mutations, but saw: " + str(ctr))
        self.assertTrue((m.chr == "21") and (m.start == 1234567), "Last mutation was not correct: " + str(m))
        
        # Reminder:  muts is a generator, so it has to be reset
        creator.reset()
        muts = creator.createMutations()
        ctr = 0
        for m in muts:
            ctr += 1
        self.assertTrue(ctr == 27, "Should have seen 27 called mutations, but saw: " + str(ctr))

    def testSimpleAnnotationWithExampleVcf(self):
        ''' Tests the ability to do a simple Gaf 3.0 annotation. '''
        inputFilename = 'testdata/vcf/example.vcf'
        outputFilename = 'out/simpleVCF.Gaf.annotated.out.tsv'
        
        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename, [])
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createGafDatasource(self.config))
        annotator.annotate()

    def testSimpleAnnotationWithGermlineVcf(self):
        ''' Tests the ability to parse Germline vcf. '''
        inputFilename = 'testdata/vcf/random.vcf'
        outputFilename = 'out/random.tsv'
        
        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename, [])
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

    def _createTCGAMAFOverridesForVCF(self):
        ''' These are the default overrides for generating a TCGA MAF file.  These will appear on all mutations, but are here for a test.
        These were taken from version 0.5.25.0 of Oncotator.
        '''
        #TODO: Remove the 'Match_Norm_Seq_Allele1' and 'Match_Norm_Seq_Allele2' from this list and populate properly, if possible.
        result = {'source':'Capture', 'status':'Somatic', 'phase':'Phase_I', 'sequencer':'Illumina GAIIx',
                'Tumor_Validation_Allele1': '', 'Tumor_Validation_Allele2': '', 'Match_Norm_Validation_Allele1': '', 'Match_Norm_Validation_Allele2': '',
                'Verification_Status': '','Validation_Status': '', 'Validation_Method': '', 'Score': '', 'BAM_file': '',
                'Match_Norm_Seq_Allele1':'', 'Match_Norm_Seq_Allele2':'','Tumor_Sample_UUID':'','Tumor_Sample_Barcode':'',
                'Strand':"+", 'Center':"broad.mit.edu", "NCBI_Build":"37"}
        return result

    def _createDatasourceCorpus(self):
        dbDir = self.config.get('DEFAULT',"dbDir")
        return DatasourceCreator.createDatasources(dbDir, "hg19",isMulticore=False)

    def testTCGAMAFRendering(self):
        ''' Test the ability to render a germline vcf as a TCGA MAF '''
        creator = VcfInputMutationCreator('testdata/vcf/example.vcf')
        creator.createMutations()
        renderer = TcgaMafOutputRenderer('out/example.vcf.maf.annotated')
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
        ''' This is a utility, private method for unit tests to get a semblance that a valid maf file was created.  
        
        Note: This method has nothing to do with the TCGA validator.
        
        TODO: This is code duplication from TCGA MAF Output RendererTest.  This should be refactored into a base class
        (to preserve self.assertTrue, etc).
        
        '''
        
        statinfo = os.stat(filename)
        self.assertTrue(statinfo.st_size > 0, "Generated MAF file (" + filename + ") is empty.")
        
        tsvReader = GenericTsvReader(filename)
        
        self.assertTrue(tsvReader.getComments().find('#version') <> -1, "First line did not specify a version number") 
        
        ctr = 1
        for lineDict in tsvReader:
            if lineDict['Entrez_Gene_Id'] == "0":
                self.assertTrue(lineDict['Hugo_Symbol'] == "Unknown", "Entrez_Gene_Id was zero, but Hugo Symbol was not 'Unknown'.  Line: " + str(ctr))
            
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
            self.assertTrue(len(unknownKeys) == 0, "__UNKNOWN__ values (" + str(len(unknownKeys)) + ") seen on line " + str(ctr) + ", in fields: " + ", ".join(unknownKeys))
            
            ctr = ctr + 1

    def testSwitchedFieldsWithExampleVcf(self):
        '''Test whether the switched tags are ignored.'''
        inputFilename = 'testdata/vcf/example.bad.switched.fields.vcf'
        outputFilename = 'out/example.out.tsv'
        
        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = SimpleOutputRenderer(outputFilename, [])
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)


    def testAnnotationWithExampleVcf(self):
        ''' Test whether parsed annotations match the actual annotations. '''
        inputFilename = 'testdata/vcf/example.vcf'
        outputFilename = 'out/example.out.tsv'
        expectedOutputFilename = 'testdata/vcf/example.expected.out.tsv'

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
        
        self.assertTrue(len(currentColNames.symmetric_difference(expectedColNames)) is 0, "Should have the same columns")
        self.assertTrue(len(current.index) == len(expected.index), "Should have the same number of rows")

        for colName in currentColNames:
            self.assertTrue(sum((current[colName] == expected[colName]) | (pandas.isnull(current[colName]) &
                                                                           pandas.isnull(expected[colName]))) ==
                            len(current.index), "Should have the same values in column " + colName)

    def testAnnotationWithNoSampleNameExampleVcf(self):
        """ Test whether parsed annotations match the actual annotations. """
        inputFilename = 'testdata/vcf/example.sampleName.removed.vcf'
        outputFilename = 'out/example.sampleName.removed.out.tsv'

        creator = VcfInputMutationCreator(inputFilename)
        renderer = SimpleOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

    def testGetMetaDataWithNoSampleNameExampleVcf(self):
        """ Make sure that we can retrieve metadata, even before createMutations has been called """
        inputFilename = 'testdata/vcf/example.sampleName.removed.vcf'

        creator = VcfInputMutationCreator(inputFilename)
        gtKeys = set(['genotype', 'read_depth', 'genotype_quality', 'haplotype_quality', 'q10', 's50', 'samples_number',
                      'depth_across_samples', 'allele_frequency', 'ancestral_allele', 'dbSNP_membership', 'id', 'qual',
                      'hapmap2_membership'])
        md = creator.getMetadata()
        ks = set(md.keys())
        diff = gtKeys.symmetric_difference(ks)
        self.assertTrue(len(diff) == 0, "Missing keys that should have been seen in the metadata: " + str(diff))

    def testSnpsAndIndelStartAndEndPos(self):
        inputFilename = "testdata/vcf/example.snps.indels.vcf"
        outputFilename = 'out/example.snps.indels.out.tsv'

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
        """ Test whether we properly determine that a field is split ... using an actual file"""
        inputFilename = 'testdata/vcf/example.split.tags.vcf'
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
                                    annotationName + " is split? " + str(isSplit[annotationName]) + ", but saw: " +
                                    str(TagConstants.SPLIT in a.getTags()))

    def testGenotypeFieldIsHonored(self):
        """Test that Oncotator does not have issues with genotype values >1 when multiple variants appear on one line"""
        inputFilename = 'testdata/vcf/example.severalGTs.vcf'
        creator = VcfInputMutationCreator(inputFilename)
        muts = creator.createMutations()
        ctr = 0
        for mut in muts:

            if MutUtils.str2bool(mut["altAlleleSeen"]):
                self.assertTrue(mut['sampleName'] != "NA 00001")
                ctr += 1
        self.assertTrue(ctr == 7, str(ctr) + " mutations with alt seen, but expected 7.  './.' should not show as a variant.")


if __name__ == "__main__":
    unittest.main()