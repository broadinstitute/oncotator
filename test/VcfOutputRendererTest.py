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
from oncotator.input.VcfInputMutationCreator import VcfInputMutationCreator
from oncotator.Annotator import Annotator
from oncotator.output.VcfOutputRenderer import VcfOutputRenderer
from TestUtils import TestUtils
from oncotator.utils.version import VERSION
import string
import logging
import vcf
TestUtils.setupLogging(__file__, __name__)


class VcfOutputRendererTest(unittest.TestCase):
    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        pass

    def tearDown(self):
        pass
    
    def _createGafDataSource(self):   
        self.logger.info("Initializing gaf 3.0")
        return TestUtils.createGafDatasource(self.config)

    def testHeaderWithExampleVcf(self):
        ver_str = string.join(["##oncotator_version", string.replace(VERSION, " ", "_")], "=")
        expected = set()
        with open('testdata/vcf/example.header.txt', 'r') as fp:
            for line in iter(fp):
                if line.startswith("##oncotator_version"):
                    expected.add(ver_str)
                else:
                    expected.add(line.rstrip('\n'))
        
        creator = VcfInputMutationCreator('testdata/vcf/example.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.header.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()
        
        current = set()
        with open('out/example.header.vcf', 'r') as fp:
            for line in iter(fp):
                if line.startswith('##'):
                    if line.startswith("##oncotator_version"):
                        current.add(ver_str)
                    else:
                        current.add(line.rstrip('\n'))

        self.assertTrue(len(current) == len(expected), "Number of lines are not the same.")
        self.assertTrue(len(current.symmetric_difference(expected)) == 0, "Headers do not match.")

    def testHeaderWithExampleVcfWithoutAnySamples(self):
        ver_str = string.join(["##oncotator_version", string.replace(VERSION, " ", "_")], "=")
        expected = set()
        with open('testdata/vcf/example.sampleName.removed.header.txt', 'r') as fp:
            for line in iter(fp):
                if line.startswith("##oncotator_version"):
                    expected.add(ver_str)
                elif line.startswith("##"):
                    expected.add(line.rstrip('\n'))

        creator = VcfInputMutationCreator('testdata/vcf/example.sampleName.removed.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.sampleName.removed.header.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        current = set()
        with open('out/example.sampleName.removed.header.vcf', 'r') as fp:
            for line in iter(fp):
                if line.startswith("##oncotator_version"):
                    current.add(ver_str)
                elif line.startswith("##"):
                    current.add(line.rstrip('\n'))

        self.assertTrue(len(current) == len(expected), "Number of lines are not the same.")
        self.assertTrue(len(current.symmetric_difference(expected)) == 0, "Headers do not match.")

    def testHeaderWithExampleVcfWithoutAnySamplesOrVariants(self):
        expected = set()
        with open('testdata/vcf/example.sampleName.variants.removed.header.txt', 'r') as fp:
            for line in iter(fp):
                    expected.add(line.rstrip('\n'))

        creator = VcfInputMutationCreator('testdata/vcf/example.sampleName.variants.removed.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.sampleName.variants.removed.header.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        current = set()
        with open('out/example.sampleName.variants.removed.header.vcf', 'r') as fp:
            for line in iter(fp):
                if line.startswith('#'):
                    current.add(line.rstrip('\n'))

    def testContentofExampleVcfwith1Variant(self):
        creator = VcfInputMutationCreator('testdata/vcf/example.1row.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.variants.1row.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename='testdata/vcf/example.1row.vcf', strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename='out/example.variants.1row.vcf', strict_whitespace=True)

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

    def testInfoContentofExampleVcf(self):
        creator = VcfInputMutationCreator('testdata/vcf/example.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.variants.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename='testdata/vcf/example.vcf', strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename='out/example.variants.vcf', strict_whitespace=True)

        self.assertEquals(expectedVcfReader.samples, currentVcfReader.samples, "Sample names do not match.")
        self.assertEquals(dict(expectedVcfReader.formats), dict(currentVcfReader.formats),
                          "Format meta-information does not match.")
        self.assertEquals(dict(expectedVcfReader.infos), dict(currentVcfReader.infos),
                          "Info meta-information does not match.")

        for expectedRecord, currentRecord in zip(expectedVcfReader, currentVcfReader):
            self.assertEqual(dict(expectedRecord.INFO), dict(currentRecord.INFO))

    def testContentofExampleVcf(self):
        creator = VcfInputMutationCreator('testdata/vcf/example.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.variants.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename='testdata/vcf/example.vcf', strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename='out/example.variants.vcf', strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testContentofExampleWithESP_MAFVcf(self):
        creator = VcfInputMutationCreator('testdata/vcf/example.withESP_MAF.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.variants.withESP_MAF.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename='testdata/vcf/example.withESP_MAF.vcf', strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename='out/example.variants.withESP_MAF.vcf', strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testContentofExampleWithVcfThatHasNoFormat(self):
        creator = VcfInputMutationCreator('testdata/vcf/example.no.format.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.no.format.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.annotate()

        expectedVcfReader = vcf.Reader(filename='testdata/vcf/example.no.format.vcf', strict_whitespace=True)
        currentVcfReader = vcf.Reader(filename='out/example.no.format.vcf', strict_whitespace=True)
        self._compareVcfs(expectedVcfReader, currentVcfReader)

    def testGafAnnotatedContentofExampleWithESP_MAFVcf(self):
        creator = VcfInputMutationCreator('testdata/vcf/example.withESP_MAF.vcf')
        creator.createMutations()
        renderer = VcfOutputRenderer('out/example.variants.gaf_annotated.withESP_MAF.vcf')
        annotator = Annotator()
        annotator.setInputCreator(creator)
        annotator.setOutputRenderer(renderer)
        annotator.addDatasource(TestUtils.createGafDatasource(self.config))
        annotator.annotate()

    def _compareGenotypeFields(self, currentSampleFields, expectedSampleFields):
        self.assertTrue("GT" in currentSampleFields, "Rendered vcf should have the field, GT.")
        self.assertTrue("GT" in expectedSampleFields, "Input vcf should have the field, GT")
        self.assertTrue("GT" == currentSampleFields[0], "Rendered vcf should have the GT field in the front of FORMAT.")
        self.assertTrue("GT" == expectedSampleFields[0], "Input vcf should have the GT field in the front of FORMAT.")
        self.assertTrue(len(set(currentSampleFields).symmetric_difference(expectedSampleFields)) >= 0,
                        "Should at least have all common fields")

    def _compareVcfs(self, expectedVcfReader, currentVcfReader):
        for expectedRecord in expectedVcfReader:
            currentRecord = currentVcfReader.next()
            self.assertTrue(expectedRecord.CHROM == currentRecord.CHROM, "Should have the same chromosome")
            self.assertTrue(expectedRecord.POS == currentRecord.POS, "Should have the same position")
            self.assertTrue(expectedRecord.ID == currentRecord.ID, "Should have the same ID")
            self.assertTrue(expectedRecord.REF == currentRecord.REF, "Should have the same reference allele")

            expectedAlts = [alt.sequence if alt is not None else None for alt in expectedRecord.ALT]
            currentAlts = [alt.sequence if alt is not None else None for alt in currentRecord.ALT]
            self.assertTrue(len(set(expectedAlts).symmetric_difference(currentAlts)) == 0,
                            "Should have the same alternate alleles")

            self.assertTrue(expectedRecord.QUAL == currentRecord.QUAL, "Should have the same qual")
            self.assertTrue(len(set(expectedRecord.FILTER).symmetric_difference(currentRecord.FILTER)) == 0,
                            "Should have the same FILTER tags")
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
                self.assertTrue(len(set(expectedVal).symmetric_difference(currentVal)) == 0,
                                "Should have the same value for INFO key, %s." % key)

            expectedSampleNames = [sample.sample for sample in expectedRecord.samples]
            currentSampleNames = [sample.sample for sample in currentRecord.samples]
            self.assertTrue(len(set(expectedSampleNames).symmetric_difference(currentSampleNames)) == 0,
                            "Should have the same number of sample names")
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
                                        "Should have the same value for genotype field, %s." % currentSampleField)
                    else:
                        self.assertTrue(len(filter(None, currentGenotypeVal)) == 0,
                                        "Rendered vcf should have missing value for genotype field, %s."
                                        % currentSampleField)

            currentSampleFields = currentRecord.FORMAT.split(":") if currentRecord.FORMAT is not None else None
            expectedSampleFields = expectedRecord.FORMAT.split(":") if expectedRecord.FORMAT is not None else None

            if currentSampleFields is None and expectedSampleFields is not None:
                self._compareGenotypeFields(currentSampleFields, expectedSampleFields)

    def testChrom2HashCodeTable(self):
        chroms = ["1", "X", "3", "contig1", "Y", "25", "mt"]
        renderer = VcfOutputRenderer("")
        h = renderer._createChrom2HashCodeTable(chroms)
        self.assertTrue(h["1"] == 1, "For chrom 1, hash code should be 1 but it was %s." % h["1"])
        self.assertTrue(h["3"] == 3, "For chrom 3, hash code should be 3 but it was %s." % h["3"])
        self.assertTrue(h["25"] == 25, "For chrom 25, hash code should be 25 but it was %s." % h["25"])
        self.assertTrue(h["X"] == 26, "For chrom X, hash code should be 26 but it was %s." % h["X"])
        self.assertTrue(h["Y"] == 27, "For chrom Y, hash code should be 27 but it was %s." % h["Y"])
        self.assertTrue(h["mt"] == 28, "For chrom mt, hash code should be 28 but it was %s." % h["mt"])
        self.assertTrue(h["contig1"] == 29, "For chrom contig1, hash code should be 29 but it was %s." % h["contig1"])

        chroms = ["contig1", "mt"]
        h = renderer._createChrom2HashCodeTable(chroms)
        self.assertTrue(h["mt"] == 3, "For chrom mt, hash code should be 3 but it was %s." % h["mt"])
        self.assertTrue(h["contig1"] == 4, "For chrom contig1, hash code should be 4 but it was %s." % h["contig1"])



if __name__ == "__main__":
    unittest.main()