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


import unittest
from oncotator.input.VcfInputMutationCreator import VcfInputMutationCreator
from oncotator.Annotator import Annotator
from oncotator.output.VcfOutputRenderer import VcfOutputRenderer
from TestUtils import TestUtils
from oncotator.utils.version import VERSION
import string
import logging
import vcf
import os
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.utils.OptionConstants import OptionConstants
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.utils.GenericTsvReader import GenericTsvReader


TestUtils.setupLogging(__file__, __name__)


class VcfOutputRendererTest(unittest.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass

    def testRemoveAnnotationsThatBeginWithUnderscore(self):
        """
        Tests that any annotation that begins with "_" does not appear in the output.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])
        outputFilename = os.path.join("out", "example.out.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        vcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        for key in vcfReader.infos.keys():
            self.assertTrue(not key.startswith("_"), "INFO tag %s begins with _." % key)

        for key in vcfReader.formats.keys():
            self.assertTrue(not key.startswith("_"), "FORMAT tag %s begins with _." % key)

        for key in vcfReader.filters.keys():
            self.assertTrue(not key.startswith("_"), "FILTER tag %s begins with _." % key)

    def testHeaderWithExampleVcf(self):
        """
        Tests that the output header is rendered correctly when the input is a VCF file.
        """
        headerFilename = os.path.join(*["testdata", "vcf", "example.header.txt"])
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])
        outputFilename = os.path.join("out", "example.header.vcf")

        ver_str = string.join(["##oncotator_version", string.replace(VERSION, " ", "_")], "=")
        expected = set()
        with open(headerFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith("##oncotator_version"):
                    expected.add(ver_str)
                else:
                    expected.add(line.rstrip('\n'))
        
        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()
        
        current = set()
        with open(outputFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith('##'):
                    if line.startswith("##oncotator_version"):
                        current.add(ver_str)
                    else:
                        current.add(line.rstrip('\n'))

        self.assertTrue(len(current) == len(expected), "Number of lines is not the same.")
        self.assertTrue(len(current.symmetric_difference(expected)) == 0, "Headers do not match.")

    def testHeaderWithExampleVcfWithoutAnySamples(self):
        """
        Tests that the output header is rendered correctly when no samples exist in the input VCF file.
        """
        headerFilename = os.path.join(*["testdata", "vcf", "example.sampleName.removed.header.txt"])
        inputFilename = os.path.join(*["testdata", "vcf", "example.sampleName.removed.vcf"])
        outputFilename = os.path.join("out", "example.sampleName.removed.header.vcf")

        ver_str = string.join(["##oncotator_version", string.replace(VERSION, " ", "_")], "=")
        expected = set()
        with open(headerFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith("##oncotator_version"):
                    expected.add(ver_str)
                elif line.startswith("##"):
                    expected.add(line.rstrip('\n'))

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        current = set()
        with open(outputFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith("##oncotator_version"):
                    current.add(ver_str)
                elif line.startswith("##"):
                    current.add(line.rstrip('\n'))

        self.assertTrue(len(current) == len(expected), "Number of lines is not the same.")
        self.assertTrue(len(current.symmetric_difference(expected)) == 0, "Headers do not match.")

    def testHeaderWithExampleVcfWithoutAnySamplesOrVariants(self):
        """
        Tests that the output header is rendered correctly when no variants exist in the input VCF file.
        """
        headerFilename = os.path.join(*["testdata", "vcf", "example.sampleName.variants.removed.header.txt"])
        inputFilename = os.path.join(*["testdata", "vcf", "example.sampleName.variants.removed.vcf"])
        outputFilename = os.path.join("out", "example.sampleName.variants.removed.header.vcf")

        expected = set()
        with open(headerFilename, 'r') as fp:
            for line in iter(fp):
                expected.add(line.rstrip('\n'))

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        current = set()
        with open(outputFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith('#'):
                    current.add(line.rstrip('\n'))

    def testContentofExampleVcfwith1Variant(self):
        """
        Tests that the output VCF file is rendered correctly when the input is a VCF file with only one variant.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.1row.vcf"])
        outputFilename = os.path.join("out", "example.variants.1row.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=inputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)

        self.assertEquals(expectedVcfReader.samples, currentVcfReader.samples, "Sample names do not match.")
        self.assertTrue(len(set(expectedVcfReader.formats.keys()) - set(currentVcfReader.formats.keys())) == 0)
        for k in expectedVcfReader.formats.keys():
            self.assertTrue(expectedVcfReader.formats[k] == currentVcfReader.formats[k], "Value is not the same for format of " + k + " : " + str(expectedVcfReader.formats[k]) + " to " + str(currentVcfReader.formats[k]))
        self.assertTrue(len(set(expectedVcfReader.infos.keys()) - set(currentVcfReader.infos.keys())) == 0)
        for k in expectedVcfReader.infos.keys():
            self.assertTrue(expectedVcfReader.infos[k] == currentVcfReader.infos[k], "Value is not the same for infos of " + k + " : " + str(expectedVcfReader.infos[k]) + " to " + str(currentVcfReader.infos[k]))

        for expectedRecord, currentRecord in zip(expectedVcfReader, currentVcfReader):
            self.assertEqual(dict(expectedRecord.INFO), dict(currentRecord.INFO))
            self.assertEquals(expectedRecord.samples, currentRecord.samples)

            for expectedCall, currentCall in zip(expectedRecord.samples, currentRecord.samples):
                for field in expectedCall.data._fields:
                    self.assertIn(field, currentCall.data._fields, "Format field %s is missing in output vcf." % field)
                    self.assertEqual(getattr(expectedCall.data, field), getattr(currentCall.data, field),
                                     "Format field %s values do not match." % field)

                for field in currentCall.data._fields:
                    if not hasattr(expectedCall.data, field):
                        val = getattr(currentCall.data, field)
                        if not isinstance(val, list):
                            val = [val]
                        self.assertEqual(filter(None, val), [], "Format field %s values do not match." % field)

    def testMetaDataofExampleVcf(self):
        """
        Tests that the meta-data (sample names, INFO and FORMAT sections, etc.) of the output VCF file is rendered
        correctly when the input is a VCF file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])
        outputFilename = os.path.join("out", "example.test_metadata.variants.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=inputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)

        self.assertEquals(expectedVcfReader.samples, currentVcfReader.samples, "Sample names do not match.")
        self.assertEquals(dict(expectedVcfReader.formats), dict(currentVcfReader.formats),
                          "Format meta-information does not match.")
        self.assertEquals(dict(expectedVcfReader.infos), dict(currentVcfReader.infos),
                          "Info meta-information does not match.")

        for expectedRecord, currentRecord in zip(expectedVcfReader, currentVcfReader):
            self.assertEqual(dict(expectedRecord.INFO), dict(currentRecord.INFO))

    def testContentofExampleVcf(self):
        """
        Tests that the output VCF file is rendered correctly when the input is a VCF file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])
        outputFilename = os.path.join("out", "example.content_of_example.variants.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=inputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testContentofExampleWithESP_MAFVcf(self):
        """
        Tests that the output VCF file is rendered correctly when the input is a VCF file with ESP MAF annotations.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.withESP_MAF.vcf"])
        outputFilename = os.path.join("out", "example.variants.withESP_MAF.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=inputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testContentofExampleWithVcfThatHasNoFormat(self):
        """
        Tests that the output VCF file is rendered correctly when the input is a vcf that is missing a FORMAT section.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.no.format.vcf"])
        outputFilename = os.path.join("out", "example.no.format.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=inputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    @TestUtils.requiresDefaultDB()
    def testGafAnnotatedContentofExampleWithESP_MAFVcf(self):
        """
        Tests that the output VCF file is rendered correctly when the input is a VCF file with ESP MAF and Gaf
        annotations.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.withESP_MAF.vcf"])
        outputFilename = os.path.join("out", "example.variants.gaf_annotated.withESP_MAF.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createTranscriptProviderDatasource(self.config))
        annotator.annotate()

    def _compareGenotypeFields(self, currentSampleFields, expectedSampleFields):
        """
        Tests the equality of rendered and expected genotype fields.

        :param currentSampleFields: current genotype fields
        :param expectedSampleFields: expected genotype fields
        """
        self.assertTrue("GT" in currentSampleFields, "Rendered vcf should have the field, GT.")
        self.assertTrue("GT" in expectedSampleFields, "Input vcf should have the field, GT")
        self.assertTrue("GT" == currentSampleFields[0], "Rendered vcf should have the GT field in the front of FORMAT.")
        self.assertTrue("GT" == expectedSampleFields[0], "Input vcf should have the GT field in the front of FORMAT.")
        self.assertTrue(len(set(currentSampleFields).symmetric_difference(expectedSampleFields)) >= 0,
                        "Should at least have all common fields")

    def _compareVcfs(self, expectedVcfReader, currentVcfReader):
        """
        Tests the equality of rendered and expected VCF files.

        :param expectedVcfReader: expected vcf
        :param currentVcfReader: rendered vcf
        """
        for expectedRecord in expectedVcfReader:
            currentRecord = currentVcfReader.next()
            self.assertTrue(expectedRecord.CHROM == currentRecord.CHROM, "Should have the same chromosome")
            self.assertTrue(expectedRecord.POS == currentRecord.POS,
                            "Should have the same position; expected: %s, saw: %s"
                            % (expectedRecord.POS, currentRecord.POS))
            self.assertTrue(expectedRecord.ID == currentRecord.ID, "Should have the same ID")
            self.assertTrue(expectedRecord.REF == currentRecord.REF, "Should have the same reference allele")

            expectedAlts = [alt.sequence if alt is not None else None for alt in expectedRecord.ALT]
            currentAlts = [alt.sequence if alt is not None else None for alt in currentRecord.ALT]
            self.assertTrue(len(expectedAlts) == len(currentAlts), "Should have the same number of alternate alleles")
            self.assertTrue(len(set(expectedAlts).symmetric_difference(currentAlts)) == 0,
                            "Should have the same alternate alleles")

            self.assertTrue(expectedRecord.QUAL == currentRecord.QUAL, "Should have the same qual")

            self.assertTrue(len(expectedRecord.FILTER) == len(currentRecord.FILTER),
                            "Should have the same number of FILTER tags")
            self.assertTrue(len(set(expectedRecord.FILTER).symmetric_difference(currentRecord.FILTER)) == 0,
                            "Should have the same FILTER tags")
            self.assertTrue(len(expectedRecord.INFO.keys()) == len(currentRecord.INFO.keys()),
                            "Should have the same number of INFO keys")
            self.assertTrue(len(set(expectedRecord.INFO.keys()).symmetric_difference(currentRecord.INFO.keys())) == 0,
                            "Should have the same INFO keys")
            keys = currentRecord.INFO.keys()
            for key in keys:
                expectedVal = expectedRecord.INFO[key]
                currentVal = currentRecord.INFO[key]
                if not isinstance(expectedVal, list):
                    expectedVal = [expectedVal]
                if not isinstance(currentVal, list):
                    currentVal = [currentVal]
                self.assertTrue(len(expectedVal) == len(currentVal),
                                "Should have the same number of value for INFO key, %s." % key)
                self.assertTrue(len(set(expectedVal).symmetric_difference(currentVal)) == 0,
                                "Should have the same value for INFO key, %s." % key)

            expectedSampleNames = [sample.sample for sample in expectedRecord.samples]
            currentSampleNames = [sample.sample for sample in currentRecord.samples]
            self.assertTrue(len(expectedSampleNames) == len(currentSampleNames),
                            "Should have the same number of sample names")
            self.assertTrue(len(set(expectedSampleNames).symmetric_difference(currentSampleNames)) == 0,
                            "Should have the sample names")
            self.assertTrue(len(expectedSampleNames) ==
                            sum([1 for i, j in zip(expectedSampleNames, currentSampleNames) if i == j]),
                            "Should have the sample names in the same order")

            # Current, as a consequence of the way VCF is rendered, will have more fields than expected
            for i in xrange(len(currentSampleNames)):
                currentSample = currentRecord.samples[i]
                expectedSample = expectedRecord.samples[i]

                currentSampleFields = currentSample.data._fields
                expectedSampleFields = expectedSample.data._fields

                self._compareGenotypeFields(currentSampleFields, expectedSampleFields)

                currentGenotypeData = currentRecord.genotype(currentSampleNames[i])
                expectedGenotypeData = expectedRecord.genotype(currentSampleNames[i])

                for currentSampleField in currentSampleFields:
                    currentGenotypeVal = currentGenotypeData[currentSampleField]
                    if not isinstance(currentGenotypeVal, list):
                        currentGenotypeVal = [currentGenotypeVal]

                    if currentSampleField in expectedSampleFields:
                        expectedGenotypeVal = expectedGenotypeData[currentSampleField]
                        if not isinstance(expectedGenotypeVal, list):
                            expectedGenotypeVal = [expectedGenotypeVal]
                        self.assertTrue(len(set(expectedGenotypeVal).symmetric_difference(currentGenotypeVal)) == 0,
                                        "Should have the same value for genotype field, %s. Should be %s but was %s." %
                                        (currentSampleField, string.join(map(str, expectedGenotypeVal), ","),
                                        string.join(map(str, currentGenotypeVal), ",")))
                    else:
                        self.assertTrue(len(filter(None, currentGenotypeVal)) == 0,
                                        "Rendered vcf should have missing value for genotype field, %s."
                                        % currentSampleField)

            currentSampleFields = currentRecord.FORMAT.split(":") if currentRecord.FORMAT is not None else None
            expectedSampleFields = expectedRecord.FORMAT.split(":") if expectedRecord.FORMAT is not None else None

            if currentSampleFields is None and expectedSampleFields is not None:
                self._compareGenotypeFields(currentSampleFields, expectedSampleFields)

    def testMissingGenotypeTag(self):
        """
        Tests that the VCF file created from a MAF file contains missing genotype information in the FORMAT fields.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.pair_name.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.missing_gt.pair_name.vcf")

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        vcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self.assertTrue("GT" in vcfReader.formats, "Vcf is missing FORMAT the following field: GT")
        for record in vcfReader:
            for sampleName in vcfReader.samples:
                genotypeData = record.genotype(sampleName)
                self.assertTrue(genotypeData["GT"] is None, "%s must have a missing genotype." % sampleName)

    def testMaf2Vcf_PairNameAnnnotationExist(self):
        """
        Tests MAF to VCF file conversion when "PairName" annotation exists in the input MAF file.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.pair_name.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.pair_name.vcf")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "maf2vcf.example.pair_name.vcf"])

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=expectedOutputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testMaf2Vcf_OnlyNormalAndTumorSampleBarcodeExist(self):
        """
        Tests MAF to VCF file conversion when "PairName" annotation is missing but only normal and tumor sample names
        exist.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.normal_tumor_sample_name.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.normal_tumor_sample_name.vcf")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "maf2vcf.example.normal_tumor_sample_name.vcf"])

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=expectedOutputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testMaf2Vcf_OnlyNormalSampleBarcodeExist(self):
        """
        Tests MAF to VCF file conversion when "PairName" annotation is missing but only normal sample names exist.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.normal_sample_name.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.normal_sample_name.vcf")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "maf2vcf.example.normal_sample_name.vcf"])

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=expectedOutputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testMaf2Vcf_OnlyTumorSampleBarcodeExist(self):
        """
        Tests MAF to VCF file conversion when "PairName" annotation is missing but only tumor sample names exist.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.tumor_sample_name.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.tumor_sample_name.vcf")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "maf2vcf.example.tumor_sample_name.vcf"])

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=expectedOutputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testContigsInVcf2Vcf(self):
        """
        Tests contig information in the header is rendered correctly.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.contigs.vcf"])
        outputFilename = os.path.join("out", "example.contigs.variants.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expected = set()
        with open(inputFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith("##contig=<ID="):
                    expected.add(line.rstrip('\n'))

        current = set()
        with open(outputFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith("##contig=<ID="):
                    current.add(line.rstrip('\n'))

        self.assertTrue(len(current) == len(expected), "Number of lines of contig information is not the same.")
        self.assertTrue(len(current.symmetric_difference(expected)) == 0, "Lines of contig information do not match.")

    def testAltsInVcf2Vcf(self):
        """
        Tests alts information in the header is rendered correctly.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.contigs.alts.vcf"])
        outputFilename = os.path.join("out", "example.contigs.alts.variants.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expected = set()
        with open(inputFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith("##ALT=<"):
                    expected.add(line.rstrip('\n'))

        current = set()
        with open(outputFilename, 'r') as fp:
            for line in iter(fp):
                if line.startswith("##ALT=<"):
                    current.add(line.rstrip('\n'))

        self.assertTrue(len(current) == len(expected), "Number of lines of alts information is not the same.")
        self.assertTrue(len(current.symmetric_difference(expected)) == 0, "Lines of alts information do not match.")

    @TestUtils.requiresDefaultDB()
    def testINSMaf2Vcf(self):
        """
        Tests that the insertions are rendered correctly when the input is a MAF file.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.normal_tumor_sample_name.ins.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.normal_tumor_sample_name.ins.vcf")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "maf2vcf.example.normal_tumor_sample_name.ins.vcf"])

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createReferenceDatasource(self.config))
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=expectedOutputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)


    @TestUtils.requiresDefaultDB()
    def testDELMaf2Vcf(self):
        """
        Tests that the deletions are rendered correctly when the input is a MAF file.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.normal_tumor_sample_name.del.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.normal_tumor_sample_name.del.vcf")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "maf2vcf.example.normal_tumor_sample_name.del.vcf"])

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createReferenceDatasource(self.config))
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=expectedOutputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    @TestUtils.requiresDefaultDB()
    def testINSVcf2Vcf(self):
        """
        Tests that the insertions are rendered correctly when the input is a VCF file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.normal_tumor_sample_name.ins.vcf"])
        outputFilename = os.path.join("out", "vcf2vcf.example.normal_tumor_sample_name.ins.vcf")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "vcf2vcf.example.normal_tumor_sample_name.ins.vcf"])

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createReferenceDatasource(self.config))
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=expectedOutputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    @TestUtils.requiresDefaultDB()
    def testDELVcf2Vcf(self):
        """
        Tests that the deletions are rendered correctly when the input is a VCF file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.normal_tumor_sample_name.del.vcf"])
        outputFilename = os.path.join("out", "vcf2vcf.example.normal_tumor_sample_name.del.vcf")
        expectedOutputFilename = os.path.join(*["testdata", "vcf", "vcf2vcf.example.normal_tumor_sample_name.del.vcf"])

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createReferenceDatasource(self.config))
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=expectedOutputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    @TestUtils.requiresDefaultDB()
    def testWhitespaceInAnnotationName(self):
        """
        Tests that the INFO keys contain no white-space, semi-colons, or equals-signs.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.with_space.tsv"])
        outputFilename = os.path.join("out", "example.with_space.out.vcf")

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createReferenceDatasource(self.config))
        annotator.annotate()

        vcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        chars = ["=", ";", " ", ":"]
        for record in vcfReader:
            for key in record.INFO.keys():
                for char in chars:
                    msg = "The %s INFO key in record (chr:%s, pos:%s) contains %s character." \
                          % (key, record.CHROM, record.POS, char)
                    self.assertTrue(key.find(char) == -1, msg)

    def testDuplicateAnnotationException(self):
        """
        Tests that the DuplicateAnnotationException is not thrown when INFO and FORMAT sections have the same key.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.duplicate_annotation.vcf"])
        outputFilename = os.path.join("out", "example.duplicate_annotation.out.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        vcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self.assertTrue("SS" in vcfReader.infos, "SS is missing in INFO.")
        self.assertTrue("SS" in vcfReader.formats, "SS is missing in FORMAT.")

    @TestUtils.requiresDefaultDB()
    def testMaf2VcfMissingFilterAnnotations(self):
        """
        Tests that the missing FILTER fields are rendered correctly when the input is a MAF file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.expected.out.tsv"])
        outputFilename = os.path.join("out", "example.expected.out.vcf")

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createReferenceDatasource(self.config))
        annotator.annotate()

        vcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self.assertTrue("s50" in vcfReader.filters, "s50 is missing in FILTER.")
        self.assertTrue("q10" in vcfReader.filters, "q10 is missing in FILTER.")

    def testMissingFilters(self):
        """
        Tests that the missing FILTER fields are rendered correctly when the input is a VCF file.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.missing_filters.vcf"])
        outputFilename = os.path.join("out", "example.missing_filters.out.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        vcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        for record in vcfReader:

            if record.CHROM == "20" and record.POS == 14370:
                self.assertEqual(record.FILTER, None, "")

            if record.CHROM == "20" and record.POS == 14370:
                self.assertEqual(record.FILTER, None, "")

    def testMaf2VcfCommentInHeader(self):
        """
        Tests MAF file to VCF file conversion when the MAF header contains comments.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.normal_sample_name.with_comments.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.normal_sample_name.with_comments.vcf")

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        vcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self.assertTrue("comment" in vcfReader.metadata, "comment section is missing in the meta-information lines.")
        self.assertItemsEqual(["_muTector_v1.0.33440", "_Oncotator_v0.5.23.0|GAF_2.1_hg19_Jun2011|dbSNP_build_134"],
                              vcfReader.metadata["comment"],
                              "comment section in the meta-information lines is incorrect.")

    def testMaf2VcfInferGenotypes(self):
        """
        Tests that the output VCF file has correctly inferred genotypes when the input is a MAF file.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "example.normal_sample_name.maf"])
        outputFilename = os.path.join("out", "maf2vcf.example.normal_sample_name.infer_genotypes.vcf")

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename, otherOptions={OptionConstants.VCF_OUT_INFER_GENOTYPES: "true"})
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        vcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        for rec in vcfReader:
            for sample in rec.samples:
                if sample.sample == "NA0001-NT":
                    genotypeData = rec.genotype(sample.sample)
                    if rec.CHROM == "20" and rec.POS == 14370:
                        msg = "Genotype for %s at chrom:%s and pos:%s was %s but should have beeen %s." % \
                              ("NA0001-NT", rec.CHROM, str(rec.POS), genotypeData["GT"], "0/2")
                        self.assertEqual(genotypeData["GT"], "0/2", msg)

                    if rec.CHROM == "21" and rec.POS == 123090:
                        msg = "Genotype for %s at chrom:%s and pos:%s was %s but should have beeen %s." % \
                              ("NA0001-NT", rec.CHROM, str(rec.POS), genotypeData["GT"], "0/1")
                        self.assertEqual(genotypeData["GT"], "0/1", msg)

                    if rec.CHROM == "22" and rec.POS == 3239880:
                        msg = "Genotype for %s at chrom:%s and pos:%s was %s but should have beeen %s." % \
                              ("NA0001-NT", rec.CHROM, str(rec.POS), genotypeData["GT"], "None")
                        self.assertEqual(genotypeData["GT"], None, msg)

    def testVcf2VcfInferGenotypes(self):
        """
        Tests that the output VCF file is rendered correctly when the input is a VCF file and infer genotypes is true.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.vcf"])
        outputFilename = os.path.join("out", "example.variants.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename, otherOptions={"infer_genotypes": True})
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=inputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testMaf2VcfInferGenotypesNoFormatFields(self):
        """
        Tests that MAF to VCF infer genotypes works when there are no format fields in the input maf
        """
        inputFilename = os.path.join(*["testdata", "maflite", "onp_combination.maf.txt"])
        outputFilename = os.path.join("out", "testMaf2VcfInferGenotypesNoFormatFields.vcf")
        creator = MafliteInputMutationCreator(inputFilename)
        renderer = VcfOutputRenderer(outputFilename, otherOptions={OptionConstants.VCF_OUT_INFER_GENOTYPES: True})
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        with open(outputFilename) as file:
            outfile = "\n".join(file.readlines())
            self.assertIn("\t0/1\t", outfile)

    @TestUtils.requiresDefaultDB()
    def testMaf2Vcf2MafConversions(self):
        """
        Tests that the MAF to VCF to MAF conversion works for insertions and deletions.
        """
        inputFilename = os.path.join(*["testdata", "maflite", "indels.verify.maf"])
        outputFilename = os.path.join("out", "indels.verify.out.vcf")

        creator = MafliteInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createReferenceDatasource(self.config))
        annotator.annotate()

        # Output Vcf file is the new input file
        creator = VcfInputMutationCreator(outputFilename)
        creator.createMutations()
        outputFilename = os.path.join(*["out", "indels.verify.maf"])
        renderer = SimpleOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        inputFileReader = GenericTsvReader(inputFilename)
        outputFileReader = GenericTsvReader(outputFilename)

        for inputRow in inputFileReader:
            outputRow = outputFileReader.next()
            self.assertEqual(inputRow["Chromosome"], outputRow["chr"], "Chromosomes do not match.")
            self.assertEqual(inputRow["Start_position"], outputRow["start"], "Start positions do not match.")
            self.assertEqual(inputRow["End_position"], outputRow["end"], "End positions do not match.")
            self.assertEqual(inputRow["Reference_Allele"], outputRow["ref_allele"], "Reference alleles do not match.")
            self.assertEqual(inputRow["Tumor_Seq_Allele2"], outputRow["alt_allele"], "Alternate alleles do not match.")

    def testVcf2VcfWithMultipleRecordsWithSamePosition(self):
        """
        Tests that the output VCF file is rendered correctly when the input is a VCF file that has multiple records
        with the same position.
        """
        inputFilename = os.path.join(*["testdata", "vcf", "example.multiple_records_same_position.vcf"])
        outputFilename = os.path.join("out", "example.multiple_records_same_position.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        creator.createMutations()
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename=inputFilename, strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    @TestUtils.requiresDefaultDB()
    def test_protein_rendering_exception_for_non_primitive(self):
        """Test problematic variants from the GATK, even with VariantToAllelicPrimitive.  Tests no exception is thrown."""
        inputFilename = os.path.join(*["testdata", "vcf", "protein_rendering_exception.vcf"])
        outputFilename = os.path.join("out", "protein_rendering_exception.annotated.vcf")

        creator = VcfInputMutationCreator(inputFilename)
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        annotator.addDatasource(transcript_ds)
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        reader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        for record in reader:
            self.assertEqual(len(record.INFO["variant_classification"]), 2,
                "Length of must be 2 but was %s." % len(record.INFO["variant_classification"]))
            self.assertEqual(len(record.INFO["variant_type"]), 2,
                "Length of must be 2 but was %s." % len(record.INFO["variant_type"]))

    @TestUtils.requiresDefaultDB()
    def test_rendering_vcf_output_format_field(self):
        """Test that we do not drop FORMAT fields for a given maflite (issue 201) """
        inputFilename = os.path.join(*["testdata", "maflite", "example_input_for_vcf_out.call_stats.txt"])
        outputFilename = os.path.join("out", "example_input_for_vcf_out.call_stats.txt.vcf")
        creator = MafliteInputMutationCreator(inputFilename)
        renderer = VcfOutputRenderer(outputFilename)
        annotator = Annotator()
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        annotator.addDatasource(transcript_ds)
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        reader = vcf.Reader(filename=outputFilename, strict_whitespace=True)
        for record in reader:
            self.assertTrue(len(record.FORMAT.split(":")) == 4, "There should be four FORMAT fields, but saw %d: %s" % (len(record.FORMAT.split(":")),record.FORMAT))

if __name__ == "__main__":
    unittest.main()