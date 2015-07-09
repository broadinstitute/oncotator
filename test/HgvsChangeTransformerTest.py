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
import logging

from oncotator.MutationData import MutationData
from oncotator.MutationDataFactory import MutationDataFactory
from test.TestUtils import TestUtils
from oncotator.utils.HgvsChangeTransformer import HgvsChangeTransformer


TestUtils.setupLogging(__file__, __name__)


class HgvsChangeTransformerTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        self.hgvs_datasource = HgvsChangeTransformer()

    ### TODO need test to assert that all necessary fields are present

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_missense(self):
        #rs80358866
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '13')
        m.createAnnotation('start', 32914782)
        m.createAnnotation('end', 32914782)
        m.createAnnotation('ref_allele', 'C')
        m.createAnnotation('alt_allele', 'T')
        m.createAnnotation('variant_classification', 'Missense_Mutation')
        m.createAnnotation('annotation_transcript', 'ENST00000380152.3')
        m.createAnnotation('genome_change', 'g.chr13:32914782C>T')
        m.createAnnotation('transcript_change', 'c.6290C>T')
        m.createAnnotation('protein_change', 'p.T2097M')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr13.hg19:g.32914782C>T')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000380152.3:c.6290C>T')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000369497:p.Thr2097Met')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_nonsense(self):
        #rs35229491
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '5')
        m.createAnnotation('start', 45303809)
        m.createAnnotation('end', 45303809)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'A')
        m.createAnnotation('variant_classification', 'Nonsense_Mutation')
        m.createAnnotation('annotation_transcript', 'ENST00000303230.4')
        m.createAnnotation('genome_change', 'g.chr5:45303809G>A')
        m.createAnnotation('transcript_change', 'c.1510C>T')
        m.createAnnotation('protein_change', 'p.R504*')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr5.hg19:g.45303809G>A')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000303230.4:c.1510C>T')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000307342:p.Arg504*')

    @TestUtils.requiresDefaultDB()
    def test_annotate_renders_with_no_build(self):
        #rs148119501
        """If mutation instance being annotated does not have a build value or is '', annotate should
        return a genome_change value with just chr. i.e. chr2:g.80529551A>C vs. chr2.hg19:g.80529551A>C"""
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('chr', 'chr2')
        m.createAnnotation('start', 80529551)
        m.createAnnotation('end', 80529551)
        m.createAnnotation('ref_allele', 'A')
        m.createAnnotation('alt_allele', 'C')
        m.createAnnotation('variant_classification', 'Intron')
        m.createAnnotation('annotation_transcript', 'ENST00000402739.4')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr2:80529551A>C')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr2:g.80529551A>C')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_intron(self):
        #rs148119501
        #+ strand transcript
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '2')
        m.createAnnotation('start', 80529551)
        m.createAnnotation('end', 80529551)
        m.createAnnotation('ref_allele', 'A')
        m.createAnnotation('alt_allele', 'C')
        m.createAnnotation('variant_classification', 'Intron')
        m.createAnnotation('annotation_transcript', 'ENST00000402739.4')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr2:80529551A>C')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr2.hg19:g.80529551A>C')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000402739.4:c.1057-90785A>C')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

        #- strand transcript
        #rs78420771
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '10')
        m.createAnnotation('start', 118891993)
        m.createAnnotation('end', 118891993)
        m.createAnnotation('ref_allele', 'A')
        m.createAnnotation('alt_allele', 'G')
        m.createAnnotation('variant_classification', 'Intron')
        m.createAnnotation('annotation_transcript', 'ENST00000277905.2')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr10:118891993A>G')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr10.hg19:g.118891993A>G')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000277905.2:c.430-5T>C')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_5_utr(self):
        #rs141173433
        #negative strand
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '7')
        m.createAnnotation('start', 6865862)
        m.createAnnotation('end', 6865862)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'C')
        m.createAnnotation('variant_classification', "5'UTR")
        m.createAnnotation('annotation_transcript', 'ENST00000316731.8')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr7:6865862G>C')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr7.hg19:g.6865862G>C')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000316731.8:c.-34C>G')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

        #positive strand
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '7')
        m.createAnnotation('start', 55086964)
        m.createAnnotation('end', 55086964)
        m.createAnnotation('ref_allele', 'A')
        m.createAnnotation('alt_allele', 'T')
        m.createAnnotation('variant_classification', "5'UTR")
        m.createAnnotation('annotation_transcript', 'ENST00000275493.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr7:55086964A>T')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr7.hg19:g.55086964A>T')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000275493.2:c.-7A>T')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_3_utr(self):
        #rs143436239
        #negative strand
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '8')
        m.createAnnotation('start', 27145409)
        m.createAnnotation('end', 27145409)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'A')
        m.createAnnotation('variant_classification', "3'UTR")
        m.createAnnotation('annotation_transcript', 'ENST00000521253.1')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr8:27145409G>A')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr8.hg19:g.27145409G>A')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000521253.1:c.*220C>T')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

        #positive strand
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '7')
        m.createAnnotation('start', 55273314)
        m.createAnnotation('end', 55273314)
        m.createAnnotation('ref_allele', 'C')
        m.createAnnotation('alt_allele', 'T')
        m.createAnnotation('variant_classification', "3'UTR")
        m.createAnnotation('annotation_transcript', 'ENST00000275493.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr7:55273314C>T')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)


        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr7.hg19:g.55273314C>T')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000275493.2:c.*4C>T')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_igr(self):
        #rs112615235
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('variant_classification', 'IGR')
        m.createAnnotation('chr', '15')
        m.createAnnotation('start', 30938316)
        m.createAnnotation('end', 30938316)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'A')
        m.createAnnotation('genome_change', '')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)

        m = transcript_ds.annotate_mutation(m)
        tx = transcript_ds.get_transcript(m.get('annotation_transcript', None))
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr15.hg19:g.30938316G>A')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), '')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_silent(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('variant_classification', 'Silent')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 19549914)
        m.createAnnotation('end', 19549914)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'A')
        m.createAnnotation('annotation_transcript', 'ENST00000477853.1')
        m.createAnnotation('genome_change', 'g.chr1:19549914G>A')
        m.createAnnotation('transcript_change', 'c.2352C>T')
        m.createAnnotation('protein_change', 'p.I784I')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.19549914G>A')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000477853.1:c.2352C>T')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_splice_site(self):
        #splice site mutation occuring in intron prior to coding start position
        #rs61191258
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('variant_classification', 'Splice_Site')
        m.createAnnotation('chr', '19')
        m.createAnnotation('start', 52994576)
        m.createAnnotation('end', 52994576)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'A')
        m.createAnnotation('annotation_transcript', 'ENST00000421239.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr19:52994576G>A')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr19.hg19:g.52994576G>A')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000421239.2:c.-121-1G>A')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

        #splice site mutation occuring in intron after coding start position
        #rs144524702
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('variant_classification', 'Splice_Site')
        m.createAnnotation('chr', '5')
        m.createAnnotation('start', 484634)
        m.createAnnotation('end', 484634)
        m.createAnnotation('ref_allele', 'C')
        m.createAnnotation('alt_allele', 'T')
        m.createAnnotation('annotation_transcript', 'ENST00000264938.3')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr5:484634C>T')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr5.hg19:g.484634C>T')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000264938.3:c.932+1G>A')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_de_novo_start_OutOfFrame(self):
        #rs114472931
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 45140082)
        m.createAnnotation('end', 45140082)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'T')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('variant_classification', 'De_novo_Start_OutOfFrame')
        m.createAnnotation('annotation_transcript', 'ENST00000372237.3')
        m.createAnnotation('genome_change', 'g.chr1:45140082G>T')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.45140082G>T')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000372237.3:c.-19C>A')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_ONP_missense(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '22')
        m.createAnnotation('start', 27003913)
        m.createAnnotation('end', 27003914)
        m.createAnnotation('ref_allele', 'CC')
        m.createAnnotation('alt_allele', 'AT')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('variant_classification', 'Missense_Mutation')
        m.createAnnotation('annotation_transcript', 'ENST00000215939.2')
        m.createAnnotation('genome_change', 'g.chr22:27003913_27003914CC>AT')
        m.createAnnotation('transcript_change', 'c.371_372GG>AT')
        m.createAnnotation('protein_change', 'p.W124Y')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr22.hg19:g.27003913_27003914delinsAT')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000215939.2:c.371_372delinsAT')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000215939:p.Trp124Tyr')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_inframe_1(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '5')
        m.createAnnotation('start', 113698631)
        m.createAnnotation('end', 113698632)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'GCC')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('variant_classification', 'In_Frame_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000512097.3')
        m.createAnnotation('genome_change', 'g.chr5:113698631_113698632insGCC')
        m.createAnnotation('transcript_change', 'c.159_160insGCC')
        m.createAnnotation('protein_change', 'p.54_54A>AA')
        #m.createAnnotation('ref_context', 'CTGCAGCCGCTGCCGCCGCCGC')
        m.createAnnotation('ref_context', 'TCCTCCCCGTCTGCAGCCGCTGCCGCCGCCGCCGCTGTTTCG') # need a larger ref_context size to get the correct mapping
        
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        #this ins of GCC occurs in a GCC-repeat region and thus need to 3' adjust position for HGVS compliance
        # it is technically a duplication
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr5.hg19:g.113698641_113698643dupGCC')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000512097.3:c.169_171dupGCC')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000427120:p.Ala58dup')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_inframe_2(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '7')
        m.createAnnotation('start', 11871469)
        m.createAnnotation('end', 11871470)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'GCAGCG')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('variant_classification', 'In_Frame_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000423059.4')
        m.createAnnotation('genome_change', 'g.chr7:11871469_11871470insGCAGCG')
        m.createAnnotation('transcript_change', 'c.103_104insCGCTGC')
        m.createAnnotation('protein_change', 'p.34_35insPL')
        #m.createAnnotation('ref_context', 'cagcagcaggagcagcggcagc')
        m.createAnnotation('ref_context', 'CGCAGCCCTGCCGGCGCCCGGGCGTAGCAGCAGCAGCAGGAGCAGCGGCAGCGGCAGCGGCAGCGGCAGCAGCTGCAGGACG') # need a larger ref_context size to get the correct mapping
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr7.hg19:g.11871488_11871493dupGCAGCG')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000423059.4:c.98_103dupCGCTGC')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000406482:p.Pro33_Leu34dup')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_inframe_3(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '8')
        m.createAnnotation('start', 10467629)
        m.createAnnotation('end', 10467630)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'TTC')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('variant_classification', 'In_Frame_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000382483.3')
        m.createAnnotation('genome_change', 'g.chr8:10467629_10467630insTTC')
        m.createAnnotation('transcript_change', 'c.3978_3979insGAA')
        m.createAnnotation('protein_change', 'p.1326_1327KT>KET')
        m.createAnnotation('ref_context', 'ccttcttctgttttagtttcct')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr8.hg19:g.10467629_10467630insTTC')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000382483.3:c.3978_3979insGAA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000371923:p.Lys1326_Thr1327insGlu')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_inframe_4(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '8')
        m.createAnnotation('start', 10467628)
        m.createAnnotation('end', 10467629)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'CCC')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('variant_classification', 'In_Frame_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000382483.3')
        m.createAnnotation('genome_change', 'g.chr8:10467628_10467629insCCC')
        m.createAnnotation('transcript_change', 'c.3979_3980insGGG')
        m.createAnnotation('protein_change', 'p.1327_1327T>RA')
        m.createAnnotation('ref_context', 'cccttcttctgttttagtttcc')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr8.hg19:g.10467628_10467629insCCC')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000382483.3:c.3979_3980insGGG')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000371923:p.Thr1327delinsArgAla')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_inframe_5(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '2')
        m.createAnnotation('start', 3197914)
        m.createAnnotation('end', 3197915)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'CAT')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('variant_classification', 'In_Frame_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000398659.4')
        m.createAnnotation('genome_change', 'g.chr2:3197914_3197915insCAT')
        m.createAnnotation('transcript_change', 'c.757_758insATG')
        m.createAnnotation('protein_change', 'p.252_253insD')
        m.createAnnotation('ref_context', 'CTGTCCGTGGGCATTCTCTATG')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr2.hg19:g.3197915_3197917dupCAT')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000398659.4:c.755_757dupATG')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000381652:p.Asn252_Ala253insAsp')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_inframe_6(self):
        #This is an insertion of a STOP in between two amino acids
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 248637602)
        m.createAnnotation('end', 248637603)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'TGA')
        m.createAnnotation('variant_classification', 'In_Frame_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000359594.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('ref_context', 'AAGAAAAGTAGTAAAGGGCAAGC')
        m.createAnnotation('protein_change', 'p.317_318ins*')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.248637602_248637603insTGA')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000359594.2:c.951_952insTGA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000352604:p.Lys318*')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_inframe_7(self):
        #This is an insertion of a STOP right before a stop
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 248637605)
        m.createAnnotation('end', 248637606)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'TGA')
        m.createAnnotation('variant_classification', 'In_Frame_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000359594.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('ref_context', 'AAGAAAAGTAGTAAAGGGCAAGC')
        m.createAnnotation('protein_change', 'p.319_319*>**')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.248637605_248637606insTGA')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000359594.2:c.954_955insTGA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_frameshift(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '4')
        m.createAnnotation('start', 1388441)
        m.createAnnotation('end', 1388442)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'CG')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('variant_classification', 'Frame_Shift_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000324803.4')
        m.createAnnotation('genome_change', 'g.chr4:1388441_1388442insCG')
        m.createAnnotation('transcript_change', 'c.142_143insCG')
        m.createAnnotation('protein_change', 'p.M48fs')
        m.createAnnotation('ref_context', 'CTGCTCACACATGCCCATGTGG')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        #this ins of CG does NOT occurs next to a CG and does not need to be position adjusted
        # it is technically an insertion
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr4.hg19:g.1388441_1388442insCG')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000324803.4:c.142_143insCG')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000323978:p.Met48fs')

    @TestUtils.requiresDefaultDB()
    def test_annotate_INS_frameshift_2(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '9')
        m.createAnnotation('start', 135977871)
        m.createAnnotation('end', 135977872)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'CGCT')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('variant_classification', 'Frame_Shift_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000393160.3')
        m.createAnnotation('genome_change', 'g.chr9:135977871_135977872insCGCT')
        m.createAnnotation('transcript_change', 'c.1835_1836insAGCG')
        m.createAnnotation('protein_change', 'p.-612fs')
        m.createAnnotation('ref_context', 'ACTCGCTCCAGCGCTTGACAAT')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr9.hg19:g.135977872_135977875dupCGCT')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000393160.3:c.1832_1835dupAGCG')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000376867:p.Arg612fs')

    @TestUtils.requiresDefaultDB()
    def test_annotate_DEL_inframe(self):
        #rs141326765
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '14')
        m.createAnnotation('start', 70924869)
        m.createAnnotation('end', 70924871)
        m.createAnnotation('ref_allele', 'ATG')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'In_Frame_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000603540.1')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr14:70924869_70924871delATG')
        m.createAnnotation('transcript_change', 'c.653_655delATG')
        m.createAnnotation('protein_change', 'p.D219del')
        m.createAnnotation('ref_context', 'GTGGTGAACCATGATTTCTTCAT')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        #this deletion is straightforward, no position adjustments necessary
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr14.hg19:g.70924869_70924871delATG')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000603540.1:c.653_655delATG')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000474385:p.Asp219del')

    @TestUtils.requiresDefaultDB()
    def test_annotate_DEL_inframe_2(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '12')
        m.createAnnotation('start', 50156659)
        m.createAnnotation('end', 50156667)
        m.createAnnotation('ref_allele', 'AAGAAGAAA')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'In_Frame_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000552699.1')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr12:50156659_50156667delAAGAAGAAA')
        m.createAnnotation('transcript_change', 'c.868_876delAAGAAGAAA')
        m.createAnnotation('protein_change', 'p.KKK290del')
        m.createAnnotation('ref_context', 'TTTCTAGGATAAGAAGAAAGAGAAGAAAT')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        #this deletion is straightforward, no position adjustments necessary
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr12.hg19:g.50156659_50156667delAAGAAGAAA')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000552699.1:c.868_876delAAGAAGAAA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000446734:p.Lys290_Lys292del')

    @TestUtils.requiresDefaultDB()
    def test_annotate_DEL_inframe_3(self):
        #rs141326765
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '19')
        m.createAnnotation('start', 40900180)
        m.createAnnotation('end', 40900182)
        m.createAnnotation('ref_allele', 'TCC')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'In_Frame_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000324001.7')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr19:40900180_40900182delTCC')
        m.createAnnotation('transcript_change', 'c.4077_4079delGGA')
        m.createAnnotation('protein_change', 'p.1359_1360EE>E')
        m.createAnnotation('ref_context', 'ACTGCcctcttcctcctcctcct')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        #this deletion is straightforward, no position adjustments necessary
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr19.hg19:g.40900189_40900191delTCC')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000324001.7:c.4077_4079delGGA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000326018:p.Glu1361del')

    @TestUtils.requiresDefaultDB()
    def test_annotate_DEL_frameshift(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '19')
        m.createAnnotation('start', 11348960)
        m.createAnnotation('end', 11348960)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'Frame_Shift_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000294618.7')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr19:11348960delG')
        m.createAnnotation('transcript_change', 'c.1664delC')
        m.createAnnotation('protein_change', 'p.P555fs')
        m.createAnnotation('ref_context', 'GAGGCTGTGCGGGTACACGTA')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        #Here only the genomic change needs to get '3 shifted because the transcript is negative strand 
        #and the coding postion is already the most 3'
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr19.hg19:g.11348960delG')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000294618.7:c.1664delC')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000294618:p.Pro555fs')

    @TestUtils.requiresDefaultDB()
    def test_annotate_SNP_nonstop(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('variant_classification', 'Nonstop_Mutation')
        m.createAnnotation('chr', '7')
        m.createAnnotation('start', 55273310)
        m.createAnnotation('end', 55273310)
        m.createAnnotation('ref_allele', 'A')
        m.createAnnotation('alt_allele', 'G')
        m.createAnnotation('annotation_transcript', 'ENST00000275493.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr7:55273310A>G')
        m.createAnnotation('transcript_change', 'c.3633A>G')
        m.createAnnotation('protein_change', 'p.*1211W')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr7.hg19:g.55273310A>G')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000275493.2:c.3633A>G')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000275493:p.*1211Trpext*6') #6 new amino acids added until another stop codon is encountered
        # "p.*1211Trpext?" would describe a variant in the stop codon at position 1211 changing it to a codon for Tryptophan (Trp, W) and adding a tail of new amino acids of unknown length since the shifted frame does not contain a new stop codon.

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_DEL_1(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 248637607)
        m.createAnnotation('end', 248637607)
        m.createAnnotation('ref_allele', 'A')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'Stop_Codon_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000359594.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr1:248637607delA')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        m.createAnnotation('ref_context', 'CAAGAAAAGTAGTAAAGGGCA')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        #Here only the genomic change needs to get '3 shifted because the transcript is negative strand 
        #and the coding postion is already the most 3'
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.248637607delA')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000359594.2:c.956delA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000352604:p.*319Cysext*?')

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_DEL_2(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 248637605)
        m.createAnnotation('end', 248637606)
        m.createAnnotation('ref_allele', 'GT')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'Stop_Codon_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000359594.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr1:248637605_248637606delGT')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        m.createAnnotation('ref_context', 'CAAGAAAAGTAGTAAAGGGCA')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        #Here only the genomic change needs to get '3 shifted because the transcript is negative strand 
        #and the coding postion is already the most 3'
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.248637605_248637606delGT')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000359594.2:c.954_955delGT')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000352604:p.*319Valext*?')

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_DEL_3(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 248637608)
        m.createAnnotation('end', 248637610)
        m.createAnnotation('ref_allele', 'GTA')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'Stop_Codon_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000359594.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('ref_context', 'AAGAAAAGTAGTAAAGGGCAAGC')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.248637608_248637610delGTA')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000359594.2:c.957_*2delGTA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_DEL_4(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 248637602)
        m.createAnnotation('end', 248637610)
        m.createAnnotation('ref_allele', 'AAAGTAGTA')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'Stop_Codon_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000359594.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('ref_context', 'AAGAAAAGTAGTAAAGGGCAAGC')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.248637602_248637610delAAAGTAGTA')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000359594.2:c.951_*2delAAAGTAGTA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000352604:p.Glu317_*319delinsGluext*?')

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_DEL_5(self):
        #negative strand transcript
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '2')
        m.createAnnotation('start', 29416090)
        m.createAnnotation('end', 29416091)
        m.createAnnotation('ref_allele', 'TC')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'Stop_Codon_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000389048.3')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('ref_context', 'GCGACCGAGCTCAGGGCCCAGG')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr2.hg19:g.29416090_29416091delTC')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000389048.3:c.4862_4863delGA')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000373700:p.*1621Cysext*53')

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_DEL_6(self):
        #negative strand transcript
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '2')
        m.createAnnotation('start', 29416092)
        m.createAnnotation('end', 29416094)
        m.createAnnotation('ref_allele', 'AGG')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'Stop_Codon_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000389048.3')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('ref_context', 'GACCGAGCTCAGGGCCCAGGCTG')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr2.hg19:g.29416092_29416094delAGG')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000389048.3:c.4859_4861delCCT')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000373700:p.Pro1620_*1621delinsArgext*41')

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_DEL_7(self):
        #negative strand transcript
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DEL')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '2')
        m.createAnnotation('start', 29416085)
        m.createAnnotation('end', 29416090)
        m.createAnnotation('ref_allele', 'CGAGCT')
        m.createAnnotation('alt_allele', '-')
        m.createAnnotation('variant_classification', 'Stop_Codon_Del')
        m.createAnnotation('annotation_transcript', 'ENST00000389048.3')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('ref_context', 'AGTGTGCGACCGAGCTCAGGGCCCAG')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr2.hg19:g.29416085_29416090delCGAGCT')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000389048.3:c.4863_*5delAGCTCG')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000373700:p.*1621Trpext*39')

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_INS(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 248637606)
        m.createAnnotation('end', 248637607)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'CAT')
        m.createAnnotation('variant_classification', 'Stop_Codon_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000359594.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('ref_context', 'AAGAAAAGTAGTAAAGGGCAAGC')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.248637606_248637607insCAT')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000359594.2:c.955_956insCAT')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000352604:p.*319Serext*1')

    @TestUtils.requiresDefaultDB()
    def test_annotate_stop_codon_ONP(self):
        m = MutationDataFactory.default_create()
        m.createAnnotation('variant_type', 'DNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '1')
        m.createAnnotation('start', 248637605)
        m.createAnnotation('end', 248637606)
        m.createAnnotation('ref_allele', 'GT')
        m.createAnnotation('alt_allele', 'CC')
        m.createAnnotation('variant_classification', 'Nonstop_Mutation')
        m.createAnnotation('annotation_transcript', 'ENST00000359594.2')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('ref_context', 'ACCAAGAAAAGTAGTAAAGGGC')
        m.createAnnotation('protein_change', 'p.318_319K*>NQ')
        
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)

        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr1.hg19:g.248637605_248637606delinsCC')
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000359594.2:c.954_955delinsCC')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), 'ENSP00000352604:p.Lys318_*319delinsAsnGlnext*1')

    @TestUtils.requiresDefaultDB()
    def test_annotate_DEL_ref_hg(self):
        """Make sure that a simple HGVS annotation run can actually see ref_hg. """

        m = MutationDataFactory.default_create()
        m.chr = "2"
        m.start = "201722365"
        m.end = "201722366"
        # m.createAnnotation('variant_type', VariantClassification.VT_DEL)
        m.ref_allele = "AC"
        m.alt_allele = "-"
        m.createAnnotation('build', 'hg19')
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)

        # This test should still pass even without ref context annotation being populated.
        ref_hg_ds = TestUtils.createReferenceDatasource(self.config)
        m = ref_hg_ds.annotate_mutation(m)
        m = transcript_ds.annotate_mutation(m)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        hgvs_dict = self.hgvs_datasource.hgvs_annotate_mutation_given_tx(m, tx)


        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(hgvs_dict.get('HGVS_genomic_change', None), 'chr2.hg19:g.201722369_201722370delAC') # NOTE: This is right-shifted in HGVS
        self.assertEqual(hgvs_dict.get('HGVS_coding_DNA_change', None), 'ENST00000321356.4:c.832+74GT>-')
        self.assertEqual(hgvs_dict.get('HGVS_protein_change', None), '')


