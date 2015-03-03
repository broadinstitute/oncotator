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
import os
from oncotator.input.OnpQueue import OnpQueue

from oncotator.utils.RunSpecificationFactory import RunSpecificationFactory
from test.TestUtils import TestUtils
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.input.OnpCombiner import OnpCombiner
from oncotator.MutationData import MutationData
from oncotator.Annotator import Annotator
from oncotator.utils.OptionConstants import OptionConstants

__author__ = 'louisb'

TestUtils.setupLogging(__file__, __name__)


class OnpCombinerTest(unittest.TestCase):
    _multiprocess_can_split_ = True
    def test_output_order(self):
        """Test that indels are not output out of order"""
        inputs = [(1, 1, 1, "C", "G", "hg19"),
                    (1, 2, 2, "-", "T", "hg19"),
                    (1, 3, 3, "A", "G", "hg19")]
        self._onp_ordered_combiner_test(inputs, inputs)

    def test_output_order_chrom_boundry(self):
        """Test that indels are not output out of order"""
        inputs = [(1, 1, 1, "C", "G", "hg19"),
                    (2, 2, 2, "A", "G", "hg19"),
                    (2, 3, 3, "T", "-", "hg19"),
                    (2, 3, 3, "G", "T", "hg19")]
        expected = [(1, 1, 1, "C", "G", "hg19"),
                    (2, 2, 3, "AG", "GT", "hg19"),
                    (2, 3, 3, "T", "-", "hg19")]
        self._onp_ordered_combiner_test(inputs, expected)

    def test_more_indel_ordering(self):
        inputs = [(2,3,3,'T','A',19),
                  (2,4,4,'C','T',19),
                  (2,5,6,'-','CA',19),
                  (2,7,7,'G','A',19)]
        expected = [(2,3,4,'TC','AT',19),
                    (2,5,6,'-','CA',19),
                    (2,7,7,'G','A',19)]
        self._onp_ordered_combiner_test(inputs,expected)

    def test_complex_single_sample(self):
        '''test with multiple chromosome, overlapping indels, and existing onps'''
        # chr, start, end, ref, alt, build
        inputs = [("1", "1", "1", "G", "C", "hg19"),
                  ("2", "3", "3", "T", "A", "hg19"),
                  ("2", "4", "4", "T", "A", "hg19"),
                  ("2", "5", "6", "-", "CA", "hg19"),
                  ("2", "7", "7", "G", "A", "hg19"),
                  ("3", "12", "12", "T", "-", "hg19"),
                  ("3", "12", "12", "G", "A", "hg19"),
                  ("3", "13", "13", "A", "T", "hg19"),
                  ("4", "100", "102", "AAA", "TTT", "hg19"),
                  ("4", "101", "101", "A", "G", "hg19"),
                  ("4", "102", "102", "C", "T", "hg19"),
                  ("5", "10","11","AT","GA", "hg19"),
                  ("5", "11","11","T","C","hg19"),
                  ("5", "12","14","CTT","AAA",'hg19')]
        expected = [("1", "1", "1", "G", "C", "hg19"),
                  ("2", "3", "4", "TT", "AA", "hg19"),
                  ("2", "5", "6", "-", "CA", "hg19"),
                  ("2", "7", "7", "G", "A", "hg19"),
                  ("3", "12", "12", "T", "-", "hg19"),
                  ("3", "12", "13", "GA", "AT", "hg19"),
                  ("4", "100", "102", "AAA", "TTT", "hg19"),
                  ("4", "101", "102", "AC", "GT", "hg19"),
                  ("5", "10","14","ATCTT","GAAAA", "hg19"),
                  ("5", "11","14","TCTT","CAAA","hg19")]
        self._onp_ordered_combiner_test(inputs, expected)

    def test_tnp_with_overlapping_snps(self):
        '''test that an onp with overlapping snps is dealt with correctly'''
        inputs = [(4,100,102,"AAA","TTT"),
            (4,101,101,'A','G'),
            (4,102,102,'C','T')]
        expected = [(4,100,102,"AAA","TTT"),
            (4,101,102,'AC',"GT")]
        self._onp_ordered_combiner_test(inputs, expected)




    def _tuples_to_MutationData(self, tuples):
        tuples = [map(str, tuple) for tuple in tuples]
        return map(lambda x: MutationData(*x), tuples)

    def _onp_unordered_combiner_test(self,inputs, expected):
        """Convert input and expected tuples into MutationData objects, then run the inputs through the ONP combiner on
        the inputs and compare to the expected"""
        input_muts = iter(self._tuples_to_MutationData(inputs))
        expected = self._tuples_to_MutationData(expected)
        combiner = OnpQueue(input_muts)
        results = list(combiner.get_combined_mutations())
        self.assert_mutations_match_expected(expected, results)

    def _assert_mutation_lists_equal(self, expected_muts, results):
        for expected, result in zip(expected_muts, results):
            self.assertTrue(expected.attributesEqual(result),
                            "\n%s != %s\nExpected:\n%s\n!=\nSeen:\n%s" % (result.positionStr(), expected.positionStr(),
                                                                          "\n".join(
                                                                              [m.positionStr() for m in expected_muts]),
                                                                          ("\n".join(
                                                                              [m.positionStr() for m in results]))))

    def _onp_ordered_combiner_test(self,inputs, expected):
        input_muts = iter(self._tuples_to_MutationData(inputs))
        expected_muts = self._tuples_to_MutationData(expected)
        combiner = OnpQueue(input_muts)
        results = list(combiner.get_combined_mutations())
        self._assert_mutation_lists_equal(expected_muts, results)

    def test_indels_only_onp_combiner(self):
        """test indels get output if there are no snps"""
        inputs = [(12, 102, 103, "C", "-", "hg19"),
                    (13, 102, 102, "-", "A", "hg19")]
        self._onp_ordered_combiner_test(inputs, inputs)

    def test_indels_multi_chromosome_ordering(self):
        """test that indels from multiple chromosome don't get sorted by position together"""
        inputs = [(12, 102, 102, "C", "-", "hg19"),
                    (13, 1, 1, "-", "A", "hg19"),
            (13,2,2,"-","A","hg19"),
            (14,1,1,"-","G","hg19")]
        self._onp_ordered_combiner_test(inputs, inputs)

    def test_indels_snps_adjacency_ordering(self):
        """test that indels and snps don't shift order when next to eachother"""
        inputs = [(1, 1, 2, "C", "T", "hg19"),
                    (1, 2, 3, "T", "-", "hg19")]
        self._onp_ordered_combiner_test(inputs, inputs)

    def test_combinatorial_onp_create(self):
        """test that every combination of dnps is created for multiple snps at the same site"""
        inputs = [(12, 102, 102, "C", "G", "hg19"),
                    (12, 102, 102, "C", "A", "hg19"),
                    (12, 103, 103, "A", "T", "hg19"),
                    (12, 103, 103, "A", "G", "hg19")]
        expected = [(12,102,103,"CA","AT","hg19"),
                    (12,102,103,"CA","AG", "hg19"),
                    (12,102,103,"CA","GG", "hg19"),
                    (12,102,103,"CA","GT","hg19")]
        self._onp_unordered_combiner_test(inputs, expected)

    def test_merge_existing_dnp(self):
        """test merging an existing DNP with a SNP"""
        inputs = [(1,1,2, "AA", "TT", "hg19"),
                  (1,3,3, "G", "C", "hg19")]
        output = [(1,1,3, "AAG", "TTC", "hg19")]
        self._onp_unordered_combiner_test(inputs, output)

    def test_combine_simple_no_sample(self):
        """test that combination still happens if no sample name is specified"""
        inputs = [(12, 101, 101, "C", "T", "hg19"),
                       (12, 102, 102, "G", "A", "hg19")]
        expected = [(12,101,102, "CG", "TA", "hg19")]
        self._onp_unordered_combiner_test(inputs, expected)

    def test_combine_different_chrom_no_sample(self):
        """be sure that chromosome is checked when combining snps"""
        input = [(12, 101, 101, "C", "T", "hg19"),
                 (13, 102, 102, "G", "A", "hg19")]
        self._onp_unordered_combiner_test(input, input)

    def test_onp_ignore_indels(self):
        """make sure indels aren't being combined with onps"""
        file = 'testdata/maflite/onp.indel.maf.txt'
        input = OnpCombiner(MafliteInputMutationCreator(file))
        expected = list(MafliteInputMutationCreator(file).createMutations())
        onp_muts = [mut for mut in input.createMutations()]
        self.assert_mutations_match_expected(expected=expected, result=onp_muts)

    def assert_mutations_match_expected(self, expected, result):
        """test that results match expected, ignoring order"""
        for mut in result:
            self.assertTrue( any(map(mut.attributesEqual, expected)), "result %s not in expected \n%s" % (mut.positionStr(), "\n".join(map(MutationData.positionStr, expected))))
        for mut in expected:
            self.assertTrue( any(map(mut.attributesEqual, result)), "expected %s not in result \n%s" % (mut.positionStr(), "\n".join(map(MutationData.positionStr, result))))

    def test_multi_sample_maflite(self):
        """Tests a multi sample maf with several unaligned onps"""
        input = OnpCombiner(MafliteInputMutationCreator('testdata/maflite/onp_combination.maf.txt'))
        onp_muts =list(input.createMutations())
        expected = self._tuples_to_MutationData([(1, 1, 1, "G", "A", "hg19"),
                                                 (1, 1, 4, "GTTT", "CGAA", "hg19"),
                                                 (1, 2, 3, "TT", "CC", "hg19"),
                                                 (1, 12, 12, "T", "G", "hg19"),
                                                 (2, 13, 13, "A", "G", "hg19")])
        expected_pair_names = ['P2','P1','P3','P1','P1']
        self._assert_mutation_lists_equal(onp_muts, expected)
        for mut, pair in zip(onp_muts, expected_pair_names):
            self.assertEqual(mut["Pair_Name"], pair)


    def test_rendering_combined_to_vcf(self):
        """Test that we produce a merged ONP vcf file without crashing """
        input_filename = os.path.join(*["testdata", "maflite", "onp_combination.maf.txt"])
        output_filename = os.path.join("out", "onp_combination.vcf")
        spec = RunSpecificationFactory.create_run_spec("MAFLITE","VCF", input_filename, output_filename,
                                                other_opts={OptionConstants.VCF_OUT_INFER_GENOTYPES: True,
                                                            OptionConstants.INFER_ONPS: True})
        annotator = Annotator()
        annotator.initialize(spec)
        annotator.annotate()

    def test_rendering_combined_to_tsv(self):
        """Test that we produce a merged ONP simple tsv file without crashing """
        input_filename = os.path.join(*["testdata", "maflite", "onp_combination.maf.txt"])
        output_filename = os.path.join("out", "onp_combination.tsv")
        spec = RunSpecificationFactory.create_run_spec("MAFLITE","SIMPLE_TSV",input_filename, output_filename,
                                                other_opts={OptionConstants.INFER_ONPS: True})
        annotator = Annotator()
        annotator.initialize(spec)
        annotator.annotate()

    @TestUtils.requiresDefaultDB()
    def test_single_sample_onp_combiner(self):
        """test that we can create an onp combined TCGA maf without crashing"""
        input_filename = 'testdata/maflite/onp.singlesample.maf.txt'
        output_filename = 'out/testSingleSampleOnpCombiner.maf'
        config = TestUtils.createUnitTestConfig()
        defaultdb = config.get('DEFAULT',"dbDir")
        spec = RunSpecificationFactory.create_run_spec("MAFLITE","TCGAMAF", input_filename, output_filename,datasourceDir=defaultdb,
                                                other_opts={OptionConstants.INFER_ONPS: True})
        annotator = Annotator()
        annotator.initialize(spec)
        annotator.annotate()

    def test_mutation_combiner(self):
        """Test that attributes and annotations are set properly with combine mutations"""
        mut1 = MutationData(chr=1,start=100, end=100, ref_allele="G", alt_allele="A")
        mut1.createAnnotation("SomeValue", "value1", "INPUT", "STRING", "a value")
        mut2 = MutationData(chr=1,start=101, end=101, ref_allele="C", alt_allele="T")
        mut2.createAnnotation("SomeValue", "value2", tags=["IT"])
        mut2.createAnnotation("AnotherValue","5")
        result = OnpQueue._combine_mutations([mut1, mut2])

        expected = MutationData(chr=1, start=100, end=101, ref_allele="GC", alt_allele="AT")
        expected.createAnnotation("SomeValue", "value1|value2", "INPUT", "STRING", "a value", tags=["IT"])
        expected.createAnnotation("AnotherValue", "5")
        self.assertTrue(result.attributesEqual(expected))
        self.assertEqual(result, expected)

    def test_mutation_combiner_no_mut(self):
        """Combining no mutations should return None"""
        result = OnpQueue._combine_mutations([])
        self.assertIsNone(result)
