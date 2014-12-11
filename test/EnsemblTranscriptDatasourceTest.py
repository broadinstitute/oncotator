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
import ConfigParser
import logging
import shutil

import unittest
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MutationData import MutationData
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.utils.ConfigUtils import ConfigUtils
from oncotator.utils.VariantClassification import VariantClassification
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from test.TestUtils import TestUtils, data_provider_decorator

TestUtils.setupLogging(__file__, __name__)
class EnsemblTranscriptDatasourceTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.config = TestUtils.createUnitTestConfig()

    def tearDown(self):
        pass

    def test_intitialize(self):
        """Test a simple initialization of an ensembl datasource """
        base_config_location = "testdata/ensembl/saccer/"
        config_parser = ConfigUtils.createConfigParser(base_config_location + "ensembl.config")
        title = config_parser.get("general", "title")
        version = config_parser.get("general", "version")
        src_file = config_parser.get("general", "src_file")

        ensembl_ds = EnsemblTranscriptDatasource(title=title, version=version, src_file=src_file)
        self.assertIsNotNone(ensembl_ds)
        ensembl_ds.set_tx_mode(TranscriptProvider.TX_MODE_BEST_EFFECT)
        self.assertTrue(TranscriptProvider.TX_MODE_BEST_EFFECT == ensembl_ds.get_tx_mode())

    def test_overlapping_single_transcripts(self):
        base_config_location = "testdata/ensembl/saccer/"

        ensembl_ds = DatasourceFactory.createDatasource(base_config_location + "ensembl.config", base_config_location)
        recs = ensembl_ds.get_overlapping_transcripts("I", "500", "500")
        self.assertTrue(len(recs) == 1)
        self.assertTrue(recs[0].get_gene() == 'YAL069W')

    def test_overlapping_multiple_transcripts_snp(self):
        base_config_location = "testdata/ensembl/saccer/"

        ensembl_ds = DatasourceFactory.createDatasource(base_config_location + "ensembl.config", base_config_location)
        recs = ensembl_ds.get_overlapping_transcripts("I", "550", "550")
        self.assertTrue(len(recs) == 2)
        ids = set()
        for r in recs:
            ids.add(r.get_transcript_id())

        self.assertTrue(len(ids - {'YAL069W', 'YAL068W-A'}) == 0)

    def test_overlapping_multiple_transcripts_indel(self):
        base_config_location = "testdata/ensembl/saccer/"

        ensembl_ds = DatasourceFactory.createDatasource(base_config_location + "ensembl.config", base_config_location)
        recs = ensembl_ds.get_overlapping_transcripts("I", "2500", "8000")
        self.assertTrue(len(recs) == 2)
        ids = set()
        for r in recs:
            ids.add(r.get_transcript_id())

        self.assertTrue(len(ids - {'YAL067W-A', 'YAL067C'}) == 0)

    def _create_ensembl_ds_from_saccer(self):
        gencode_input_gtf = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf"
        gencode_input_fasta = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa"
        base_output_filename = "out/test_saccer_ds"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([gencode_input_gtf], [gencode_input_fasta], base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, title="ensembl", version="71")
        return ensembl_ds

    def test_simple_annotate_with_nonhuman(self):
        """Test a very simple annotation with a nonhuman genome (saccer)"""
        ensembl_ds = self._create_ensembl_ds_from_saccer()

        m = MutationData()
        m.chr = "I"
        m.start = "500"
        m.end = "500"
        m.ref_allele = "C"
        m.alt_allele = "A"

        m2 = ensembl_ds.annotate_mutation(m)

        self.assertTrue(m2['annotation_transcript'] == "YAL069W")
        self.assertTrue(m2['gene'] == "YAL069W")

    def test_simple_annotate(self):
        """ Annotate a simple example.
        """
        base_config_location = "testdata/ensembl/saccer/"
        config_parser = ConfigUtils.createConfigParser(base_config_location + "ensembl.config")
        title = config_parser.get("general", "title")
        version = config_parser.get("general", "version")
        src_file = config_parser.get("general", "src_file")

        ensembl_ds = EnsemblTranscriptDatasource(title=title, version=version, src_file=src_file)

        m = MutationData()
        m.chr = "22"
        m.start = "22161963"
        m.end = "22161963"
        m.ref_allele = "C"
        m.alt_allele = "A"

        m2 = ensembl_ds.annotate_mutation(m)

    def test_overlapping_multiple_genes(self):
        """Test that we can collect multiple overlapping genes """
        ds = TestUtils._create_test_gencode_v19_ds("out/overlapping_genes_multiple_")
        genes = ds.get_overlapping_genes("22", 22080000, 22120000)
        self.assertTrue(len({"MAPK1", "YPEL1"} - genes) ==0 )

    def test_overlapping_gene(self):
        """Test that we can collect an overlapping gene """
        ds = TestUtils._create_test_gencode_v19_ds("out/overlapping_genes_")
        genes = ds.get_overlapping_genes("22", 22115000, 22120000)
        self.assertTrue(len({"MAPK1"} - genes) == 0)

    def test_check_for_appris_tag(self):
        """Test that a transcript with an appris tag returns the right rank"""
        ds = TestUtils._create_test_gencode_v19_ds("out/appris_tag",)
        txs = ds.get_overlapping_transcripts("22", 22222050, 22222050, padding=100)
        self.assertTrue( len(txs) == 1)
        self.assertEquals(ds._get_appris_rank(txs[0]),0)

    def test_check_for_missing_appris_tag(self):
        """Check that the correct value is returned for a site with no appris tag """
        ds = TestUtils._create_test_gencode_v19_ds("out/appris_no_tag",)
        txs = ds.get_overlapping_transcripts("16", 61556, 61556, padding=100)
        self.assertTrue( len(txs) > 0)
        self.assertEquals(ds._get_appris_rank(txs[0]), TranscriptProviderUtils.NO_APPRIS_VALUE)

    @TestUtils.requiresDefaultDB()
    def test_appris_ccds_tag(self):
        m = MutationData(chr="1", start="200818757", end="200818757", ref_allele="C", alt_allele="A", build="hg19")
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        m = transcript_ds.annotate_mutation(m)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(tx._transcript_id,'ENST00000358823.2')

    @TestUtils.requiresDefaultDB()
    def test_appris_selects_transcript(self):
        m = MutationData(chr="2", start="201722365", end="201722366", ref_allele="AC", alt_allele="-", build="hg19")
        transcript_ds = TestUtils.createTranscriptProviderDatasource(self.config)
        m = transcript_ds.annotate_mutation(m)
        tx = transcript_ds.get_transcript(m['annotation_transcript'])
        self.assertTrue(tx is not None, "Transcript was None when it should have been found.  Does the ground truth transcript above need to be updated?")
        self.assertEqual(tx._transcript_id,'ENST00000321356.4')

    def test_overlapping_gene_5flank(self):
        """Test that we can collect an overlapping gene on its 5' Flank """
        ds = TestUtils._create_test_gencode_v19_ds("out/overlapping_genes_flank")
        txs = ds.get_overlapping_transcripts("22", 22222050, 22222050, padding=100)
        self.assertTrue( len(txs) == 1)
        self.assertTrue(txs[0].get_transcript_id() == "ENST00000398822.3")

        txs = ds.get_overlapping_transcripts("22", 22224920, 22224920)
        self.assertTrue(len(txs) == 0)


    def test_small_positive_strand_transcript_change(self):
        """Test one location on a transcript and make sure that the transcript change rendered properly """
        ds = TestUtils._create_test_gencode_v19_ds("out/small_positive_strand_")

        # Now for a negative strand
        m = MutationData()
        m.chr = "22"
        m.start = "22221730"
        m.end = "22221730"
        m.ref_allele = "T"
        m.alt_allele = "G"
        m2 = ds.annotate_mutation(m)
        self.assertTrue(m2['transcript_change'] == "c.1A>C", "Incorrect transcript change: " + m2['transcript_change'])

        # positive strand
        m = MutationData()
        m.chr = "3"
        m.start = "178916614"
        m.end = "178916614"
        m.ref_allele = "G"
        m.alt_allele = "T"
        m2 = ds.annotate_mutation(m)
        self.assertTrue(m2['transcript_change'] == "c.1G>T", "Incorrect transcript change: " + m2['transcript_change'])

    def test_hgvs_annotations_simple_SNP(self):
        """Test that HGVS annotations appear (incl. protein change) in a mutation, so we believe that the Transcript objects are populated properly."""
        ds = TestUtils._create_test_gencode_v19_ds("out/test_hgvs_annotations_SNP_")

        # Now for a negative strand
        m = MutationData()
        m.chr = "22"
        m.start = "22221730"
        m.end = "22221730"
        m.ref_allele = "T"
        m.alt_allele = "G"
        m.build = "hg19"
        m2 = ds.annotate_mutation(m)
        self.assertEqual(m2.get('HGVS_genomic_change', None), 'chr22.hg19:g.22221730T>G')
        self.assertEqual(m2.get('HGVS_coding_DNA_change', None), 'ENST00000215832.6:c.1A>C')
        self.assertEqual(m2.get('HGVS_protein_change', None), 'ENSP00000215832:p.Met1Leu')

    def test_hgvs_annotations_IGR(self):
        """Test that the HGVS annotations appear for IGR"""
        ds = TestUtils._create_test_gencode_v19_ds("out/test_hgvs_annotations_IGR_")
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('variant_classification', 'IGR')
        m.createAnnotation('chr', '15')
        m.createAnnotation('start', 30938316)
        m.createAnnotation('end', 30938316)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'A')
        m2 = ds.annotate_mutation(m)
        self.assertEqual(m2.get('HGVS_genomic_change', None), 'chr15.hg19:g.30938316G>A')
        self.assertEqual(m2.get('HGVS_coding_DNA_change', None), '')
        self.assertEqual(m2.get('HGVS_protein_change', None), '')

    def test_no_mapping_file(self):
        """Test that we can still create (from scratch) and instantiate a EnsemblDatasource when no protein mapping is specified (i.e. limited HGVS support)"""
        """Test that HGVS annotations appear (incl. protein change) in a mutation, so we believe that the Transcript objects are populated properly."""
        ds = TestUtils._create_test_gencode_v19_ds("out/test_hgvs_annotations_no_mapping_file_", protein_id_mapping_file=None)

        # Now for a negative strand
        m = MutationData()
        m.chr = "22"
        m.start = "22221730"
        m.end = "22221730"
        m.ref_allele = "T"
        m.alt_allele = "G"
        m.build = "hg19"
        m2 = ds.annotate_mutation(m)
        self.assertEqual(m2.get('HGVS_genomic_change', None), 'chr22.hg19:g.22221730T>G')
        self.assertEqual(m2.get('HGVS_coding_DNA_change', None), 'ENST00000215832.6:c.1A>C')
        self.assertEqual(m2.get('HGVS_protein_change', None), 'unknown_prot_seq_id:p.Met1Leu')

    @TestUtils.requiresDefaultDB()
    def test_retrieve_transcripts_from_region(self):
        """Test that we can retrieve a large number of transcripts.  Requires a full gencode datasource."""
        config = TestUtils.createUnitTestConfig()
        transcript_ds = TestUtils.createTranscriptProviderDatasource(config)
        filtered_txs = transcript_ds.get_transcripts_by_pos(chr="1", start="1", end="100000000")

        self.assertTrue(len(filtered_txs) > 4000)
        gene_set = set([tx.get_gene() for tx in filtered_txs])
        self.assertTrue(len(gene_set) > 1500)


    segment_start_data_negative_strand = lambda: (
        # The start is between exon 0 and exon 1.  Since the start and end are in genomic space, we expect 0,"-"
        ("22", 22162050, 1, "-"),
        ("22", 22165000, 0, "-"),
        ("22", 22155000, 2, "-"),
        ("22", 22156000, 2, "-"),
        ("22", 22220000, 0, "-"),
        ("22", 22125000, 6, "-")
    )

    @TestUtils.requiresDefaultDB()
    @data_provider_decorator(segment_start_data_negative_strand)
    def test_determine_exons_affected_by_start_negative_strand(self, chrom, start, gt_exon_id, gt_exon_direction):

        chosen_tx, transcript_ds = self._get_chosen_tx_and_transcript_ds(chrom, start)

        result_tuple = transcript_ds._determine_exons_affected_by_start(start, chosen_tx)

        self.assertTrue(result_tuple[0] == gt_exon_id, "GT did not match guess exon ID... GT exon ID: %d    Seen: %d " % (gt_exon_id, result_tuple[0]))
        self.assertTrue(result_tuple[1] == gt_exon_direction)

    segment_end_data_negative_strand = lambda: (
        ("22", 22123000, 8, "+"),
        ("22", 22165000, 1, "+"),
        ("22", 22155000, 3, "+"),
        ("22", 22156000, 3, "+"),
        ("22", 22220000, 1, "+"),
        ("22", 22125000, 7, "+"),
        ("22", 22162050, 1, "+"), # in exon
    )

    @TestUtils.requiresDefaultDB()
    @data_provider_decorator(segment_end_data_negative_strand)
    def test_determine_exons_affected_by_end_negative_strand(self, chrom, end, gt_exon_id, gt_exon_direction):

        chosen_tx, transcript_ds = self._get_chosen_tx_and_transcript_ds(chrom, end)

        result_tuple = transcript_ds._determine_exons_affected_by_end(end, chosen_tx)

        self.assertTrue(result_tuple[0] == gt_exon_id, "GT did not match guess exon ID... GT exon ID: %d    Seen: %d " % (gt_exon_id, result_tuple[0]))
        self.assertTrue(result_tuple[1] == gt_exon_direction)


    segment_start_data_positive_strand = lambda: (
        ("3", 178920000, 4, "+"),
        ("3", 178921000, 4, "+"),
        ("3", 178919500, 4, "+"),
        ("3", 178917500, 2, "+"), # in exon
    )

    @TestUtils.requiresDefaultDB()
    @data_provider_decorator(segment_start_data_positive_strand)
    def test_determine_exons_affected_by_start_positive_strand(self, chrom, start, gt_exon_id, gt_exon_direction):

        chosen_tx, transcript_ds = self._get_chosen_tx_and_transcript_ds(chrom, start)

        result_tuple = transcript_ds._determine_exons_affected_by_start(start, chosen_tx)

        self.assertTrue(result_tuple[0] == gt_exon_id, "GT did not match guess exon ID... GT exon ID: %d    Seen: %d " % (gt_exon_id, result_tuple[0]))
        self.assertTrue(result_tuple[1] == gt_exon_direction)

    segment_end_data_positive_strand = lambda: (
        ("3", 178920000, 3, "-"),
        ("3", 178921000, 3, "-"),
        ("3", 178919500, 3, "-"),
        ("3", 178917500, 2, "-"), # in exon
    )

    def _get_chosen_tx_and_transcript_ds(self, chrom, loc):
        config = TestUtils.createUnitTestConfig()
        transcript_ds = TestUtils.createTranscriptProviderDatasource(config)
        transcript_ds.set_tx_mode(TranscriptProvider.TX_MODE_CANONICAL)
        start_txs = transcript_ds.get_transcripts_by_pos(chr=chrom, start=str(loc), end=str(loc))
        chosen_tx = transcript_ds._choose_transcript(start_txs, transcript_ds.get_tx_mode(),
                                                     VariantClassification.VT_SNP, "", "", str(loc), str(loc))
        return chosen_tx, transcript_ds

    @TestUtils.requiresDefaultDB()
    @data_provider_decorator(segment_end_data_positive_strand)
    def test_determine_exons_affected_by_end_positive_strand(self, chrom, end, gt_exon_id, gt_exon_direction):

        chosen_tx, transcript_ds = self._get_chosen_tx_and_transcript_ds(chrom, end)

        result_tuple = transcript_ds._determine_exons_affected_by_end(end, chosen_tx)

        self.assertTrue(str(result_tuple[0])+result_tuple[1] == str(gt_exon_id)+gt_exon_direction, "GT did not match guess exon ID... GT exon ID: %s    Seen: %s " % (str(gt_exon_id) + gt_exon_direction, str(result_tuple[0]) + result_tuple[1]))

    segment_igr_overlaps = lambda: (
        ("3", 178990000, -1, ""), # IGR
        ("22", 22100000, -1, "")
    )

    @TestUtils.requiresDefaultDB()
    @data_provider_decorator(segment_igr_overlaps)
    def test_determine_exons_affected_for_start_for_IGR_segment(self, chrom, start, gt_exon_id, gt_exon_direction):
        """Test exon inclusion for a segment that has a start position in IGR"""
        chosen_tx, transcript_ds = self._get_chosen_tx_and_transcript_ds(chrom, start)

        result_tuple = transcript_ds._determine_exons_affected_by_start(start, chosen_tx)

        self.assertTrue(result_tuple is None, "Result should have been None for IGR overlap, but saw: %s " % str(result_tuple))

    @TestUtils.requiresDefaultDB()
    @data_provider_decorator(segment_igr_overlaps)
    def test_determine_exons_affected_for_end_for_IGR_segment(self, chrom, start, gt_exon_id, gt_exon_direction):
        """Test exon inclusion for a segment that has a start position in IGR"""
        chosen_tx, transcript_ds = self._get_chosen_tx_and_transcript_ds(chrom, start)

        result_tuple = transcript_ds._determine_exons_affected_by_end(start, chosen_tx)

        self.assertTrue(result_tuple is None, "Result should have been None for IGR overlap, but saw: %s " % str(result_tuple))

    @TestUtils.requiresDefaultDB()
    def test_continuous_exons_in_segments(self):
        """Test that all exons are accounted when annotating adjacent segments that skip an exon. """
        # SPECC1L 10+	    22	24734447	SPECC1L	10+	41783674	TEF	1-	1215.0	-0.04975556624325125		hg19	CESC.TCGA.BI.A0VR.Tumor.SM.1RACM
        # SPECC1L 8-	    22	16282318	POTEH	2-	24730543	SPECC1L	8-	433.0	-0.00781166374668759		hg19	CESC.TCGA.BI.A0VR.Tumor.SM.1RACM
        # SPECC1L-ADORA2A	22	24734447	SPECC1L	10+	41783674	TEF	1-	1215.0	-0.04975556624325125		hg19	CESC.TCGA.BI.A0VR.Tumor.SM.1RACM

        seg1 = MutationData()
        seg1.chr = "22"
        seg1.start = "24734447" # Just passed the exon 9 (0-based)
        seg1.end = "41783674"

        seg2 = MutationData()
        seg2.chr = "22"
        seg2.start = "16282318"
        seg2.end = "24730543" # Just passed the exon 8 (0-based)

        segs = [seg1, seg2]

        # 'ENST00000314328.9' for GENCODE v19
        chosen_tx, transcript_ds = self._get_chosen_tx_and_transcript_ds(seg1.chr, seg1.start)
        result_tuple = transcript_ds._determine_exons_affected_by_start(seg1.start, chosen_tx)

        self.assertTrue(result_tuple == (10, '+'))

        result_tuple = transcript_ds._determine_exons_affected_by_end(seg2.end, chosen_tx)
        self.assertTrue(result_tuple == (8, '-'))

    def test_retrieve_transcript_by_gene(self):
        """Simple test of retrieve_transcript_by_gene """
        gene = "MAPK1"

        ds = TestUtils._create_test_gencode_v19_ds("out/test_retrieve_transcript_by_gene_")
        txs = ds.retrieve_transcripts_by_gene(gene)

        self.assertTrue(len(txs) > 2)

        tx_ids = [tx.get_transcript_id() for tx in txs]

        self.assertTrue("ENST00000398822.3" in tx_ids, "ENST00000398822.3 not in gene %s -- is the version number correct?" % gene)
        self.assertTrue("ENST00000215832.6" in tx_ids, "ENST00000215832.6 not in gene %s -- is the version number correct?" % gene)

        for tx in txs:
            self.assertTrue(tx.get_gene() == gene)

    def test_arbitrary_rankings(self):
        """test that _select_best_with_multiple_criteria can sort with mutliple criteria and get the right answer"""
        a = (0,1)
        b = (1,1)
        c = (1,2)
        d = (2,1)
        e = (2,2)
        f = (0,4)
        g = (-1,5)
        input = [a,b,c,d,e,f,g]
        #sort by left minimum, right minimum
        result = EnsemblTranscriptDatasource._select_best_with_multiple_criteria(input, [(lambda x: x[0], min),(lambda x: x[1], min)])
        self.assertEqual(result[0], g)

        #sort by right minimum, left minimum
        result = EnsemblTranscriptDatasource._select_best_with_multiple_criteria(input,[(lambda x: x[1], min),(lambda x: x[0], min)])
        self.assertEqual(result[0], a)

        #sort by left maximum, right minimum
        result = EnsemblTranscriptDatasource._select_best_with_multiple_criteria(input, [(lambda x: x[0], max),(lambda x: x[0], min)])
        self.assertEqual(result[0], d)

        #sort by sum, then right maximum
        result = EnsemblTranscriptDatasource._select_best_with_multiple_criteria(input,[(sum, max), (lambda x: x[1],max)])
        self.assertEqual(result[0], g)

    def test_tie_breaking_rankings(self):
        """test that _select_best_with_multiple_criteria works with ties"""
        a = (0,0,1)
        b = (0,0,2)
        c = (0,0,3)
        input =[a,b,c]
        result = EnsemblTranscriptDatasource._select_best_with_multiple_criteria(input, [(lambda x: x[0], max),
                                                                                (lambda x: 3,min),
                                                                                (lambda x: x[1], max),
                                                                                (lambda x: x[2],max)])
        self.assertEqual(result[0],c)

    def test_canonical_tx_list(self):
        """Test that specifying the canonical list will actually change the transcript selected. """
        ds = TestUtils._create_test_gencode_v19_ds("out/test_canonical_tx_list_")
        m = MutationData()
        m.chr = "22"
        m.start = "22142650"
        m.end = "22142650"
        m.ref_allele = "T"
        m.alt_allele = "A"
        ds.set_custom_canonical_txs(["ENST00000544786"])
        ds.set_tx_mode(TranscriptProvider.TX_MODE_BEST_EFFECT)

        # NOTE: tx list overrides best effect
        m2 = ds.annotate_mutation(m)
        self.assertTrue(m2['annotation_transcript'].startswith("ENST00000544786"))
        self.assertTrue(m2['variant_classification'] == VariantClassification.INTRON)

        ds.set_custom_canonical_txs([])
        m2 = ds.annotate_mutation(m)
        self.assertTrue(m2['variant_classification'] == VariantClassification.MISSENSE)
        self.assertFalse(m2['annotation_transcript'].startswith("ENST00000544786"))


    def test_canonical_tx_list_miss(self):
        """Test that specifying the canonical list will do nothing otherwise."""
        ds = TestUtils._create_test_gencode_v19_ds("out/test_canonical_tx_list_")
        m = MutationData()
        m.chr = "22"
        m.start = "22142650"
        m.end = "22142650"
        m.ref_allele = "T"
        m.alt_allele = "A"
        ds.set_custom_canonical_txs(["ENST00000123456"])

        m2 = ds.annotate_mutation(m)
        self.assertFalse(m2['annotation_transcript'].startswith("ENST00000544786"))
        self.assertFalse(m2['variant_classification'] == VariantClassification.INTRON)

        ds.set_custom_canonical_txs([])
        m2 = ds.annotate_mutation(m)
        self.assertTrue(m2['variant_classification'] == VariantClassification.MISSENSE)
        self.assertFalse(m2['annotation_transcript'].startswith("ENST00000544786"))

    def test_canonical_tx_list_empty(self):
        """Test that not specifying the canonical list will do nothing."""
        ds = TestUtils._create_test_gencode_v19_ds("out/test_canonical_tx_list_")
        m = MutationData()
        m.chr = "22"
        m.start = "22142650"
        m.end = "22142650"
        m.ref_allele = "T"
        m.alt_allele = "A"

        m2 = ds.annotate_mutation(m)
        self.assertFalse(m2['annotation_transcript'].startswith("ENST00000544786"))
        self.assertFalse(m2['variant_classification'] == VariantClassification.INTRON)

        ds.set_custom_canonical_txs([])
        m2 = ds.annotate_mutation(m)
        self.assertTrue(m2['variant_classification'] == VariantClassification.MISSENSE)
        self.assertFalse(m2['annotation_transcript'].startswith("ENST00000544786"))

    def test_hashcode_changes_when_tx_mode_changes(self):
        """Test that a call to set_tx_mode will change the md5 hash for the datasource"""
        ds = TestUtils._create_test_gencode_v19_ds("out/test_hashcode_changes_when_tx_mode_changes_")
        ds.set_tx_mode(TranscriptProvider.TX_MODE_CANONICAL)
        dummy_seed = "dummy"
        ds.set_hashcode(dummy_seed)

        initial_hash = ds.get_hashcode()
        self.assertTrue(initial_hash != dummy_seed)

        ds.set_tx_mode(TranscriptProvider.TX_MODE_BEST_EFFECT)
        be_hash = ds.get_hashcode()
        self.assertTrue(initial_hash != be_hash)

        ds.set_tx_mode(TranscriptProvider.TX_MODE_CANONICAL)
        test_hash = ds.get_hashcode()
        self.assertTrue(test_hash == initial_hash)

        new_dummy_seed = "new_dummy"
        ds.set_hashcode(new_dummy_seed)

        # MAke sure new_dummy changes the hash.
        initial_hash2 = ds.get_hashcode()
        self.assertTrue(initial_hash2 != initial_hash)

        ds.set_tx_mode(TranscriptProvider.TX_MODE_BEST_EFFECT)
        be_hash2 = ds.get_hashcode()
        self.assertTrue(initial_hash2 != be_hash2)
        self.assertTrue(be_hash != be_hash2)

        ds.set_tx_mode(TranscriptProvider.TX_MODE_CANONICAL)
        test_hash = ds.get_hashcode()
        self.assertTrue(test_hash == initial_hash2)


if __name__ == '__main__':
    unittest.main()
