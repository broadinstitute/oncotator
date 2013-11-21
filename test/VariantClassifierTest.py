import shutil
import Bio
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources import EnsemblTranscriptDatasource
from oncotator.index.gaf import region2bin, region2bins
from oncotator.utils.VariantClassification import VariantClassification
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from test.MUC16ChangeTestdata import muc16_change_testdata
from test.TestUtils import TestUtils
from MUC16Testdata import muc16testdata
import unittest
from unittest_data_provider import data_provider
from oncotator.utils.VariantClassifier import VariantClassifier

TestUtils.setupLogging(__file__, __name__)
class VariantClassifierTest(unittest.TestCase):
    _multiprocess_can_split_ = True
    # TODO: Get recently downloaded test data and use that.

    def _create_ensembl_ds_from_testdata(self, gene, fixed_id):
        gencode_input_gtf = "testdata/gencode/" + gene + ".gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/" + gene + ".gencode.v18.pc_transcripts.fa"
        base_output_filename = "out/test_variant_classification_" + fixed_id
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([gencode_input_gtf], [gencode_input_fasta], base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, title="GENCODE", version="v18")
        return ensembl_ds

    def _test_variant_classification(self, alt, chr, end, gt_vc, ref, start, vt, fixed_id, gene="MAPK1"):
        ensembl_ds = self._create_ensembl_ds_from_testdata(gene, fixed_id)
        recs = ensembl_ds.get_overlapping_transcripts(chr, start, end)
        txs = ensembl_ds._filter_transcripts(recs)
        tx = ensembl_ds._choose_transcript(txs, EnsemblTranscriptDatasource.TX_MODE_CANONICAL, vt, ref, alt, start, end)
        self.assertTrue(len(recs) != 0, "Issue with test...No transcripts found for: " + str([chr, start, end]))

        vcer = VariantClassifier()
        variant_classification = vcer.variant_classify(tx, ref, alt, start, end, vt, dist=2)
        vc = variant_classification.get_vc()
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))
        return variant_classification, tx

    variants_indels_MAPK1 = lambda: (
        ("22", "22221645", "22221647", "In_Frame_Del", "DEL", "GAG", "-"),
        ("22", "22221645", "22221645", "Frame_Shift_Del", "DEL", "G", "-"),
        ("22",	"22221645", "22221645", "Frame_Shift_Ins", "INS", "-", "A")
    )
    @data_provider(variants_indels_MAPK1)
    def test_variant_classification_indels_simple(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test a small set of indels, from MAPK1.  These are all in frame or frame shifts cleanly within exons."""
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, fixed_id="variant_classification_indels_simple")

    variants_indels_MUC16 = lambda: (
        ("19", "8979051", "8979051", "Intron", "DEL", "A", "-"),
        ("19", "9058690", "9058690", "Frame_Shift_Del", "DEL", "T", "-"),
        ("19", "9064536", "9064536", "Frame_Shift_Del", "DEL", "T", "-"),
        ("19", "9006152", "9006153", "Intron", "INS", "-", "AA"),
        ("19", "9006225", "9006226", "Intron", "INS", "-", "CT"),
        ("19", "9073231", "9073232", "Frame_Shift_Ins", "INS", "-", "A")
    )
    @data_provider(variants_indels_MUC16)
    def test_variant_classification_indels_muc16(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test small set of indels, from MUC16.  These are all introns or frame shifts cleanly within exon /intron."""
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, fixed_id="variant_classification_indels_muc16", gene="MUC16")

    variants_snps_missense = lambda: (
        ("22", "22127164", "22127164", "Missense_Mutation", "SNP", "C", "G"),
        ("22", "22143049", "22143049", "Nonsense_Mutation", "SNP", "C", "A"),
        ("22", "22162014", "22162014", "Missense_Mutation", "SNP", "C", "T")
    )
    @data_provider(variants_snps_missense)
    def test_variant_classification_missense(self, chr, start, end, gt_vc, vt, ref, alt):
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, fixed_id="variant_classification_missense")

    def test_region_queries(self):
        b = region2bin(22221612, 22221919)
        bins = region2bins(22221645, 22221645)
        self.assertTrue(b in bins)

    frameshift_indels = lambda: (
        ("INS", 10, 11,  "A", True),
        ("INS", 10, 12,  "ATC", False),
        ("DEL", 10, 10,  "-", True),
        ("DEL", 10, 12,  "-", False)
    )

    @data_provider(frameshift_indels)
    def test_is_framshift_indel(self, vt, s, e, alt, gt):
        vcer = VariantClassifier()
        self.assertTrue(vcer.is_frameshift_indel(vt, s, e, alt) == gt)

    # ("MUC16", "19", "9057555", "9057555", "Missense_Mutation", "SNP", "C", "A", "g.chr19:9057555C>A", "-", "c.29891G>T", "c.(29890-29892)gGg>gTg", "p.G9964V"),
    @data_provider(muc16_change_testdata)
    def test_muc16_change_genome(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, strand, transcript_change_gt, codon_change_gt, protein_change_gt):
        """ Test all of the MUC16 changes (protein, genome, codon, and transcript)."""
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, fixed_id="muc16_change_genome", gene="MUC16")
        vcer = VariantClassifier()
        genome_change = TranscriptProviderUtils.determine_genome_change(chr, start, end, ref, alt, vc.get_vt())
        self.assertTrue(genome_change == genome_change_gt, "Genome change did not match gt (%s): %s" %(genome_change_gt, genome_change))

    @data_provider(muc16_change_testdata)
    def test_muc16_change_transcript(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, strand, transcript_change_gt, codon_change_gt, protein_change_gt):
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, fixed_id="muc16_change_transcript",gene="MUC16")
        vcer = VariantClassifier()
        transcript_change = vcer.generate_transcript_change_from_tx(tx, vt, vc, int(start), int(end), ref, alt)
        self.assertTrue(transcript_change == transcript_change_gt, "Transcript change did not match gt (%s): %s" % (transcript_change_gt, transcript_change))

    @data_provider(muc16_change_testdata)
    def test_muc16_change_codon(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, strand, transcript_change_gt, codon_change_gt, protein_change_gt):
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, fixed_id="muc16_change_codon",gene="MUC16")
        vcer = VariantClassifier()
        codon_change = vcer.generate_codon_change_from_vc(tx, int(start), int(end), vc)
        self.assertTrue(codon_change == codon_change_gt, "Codon change did not match gt (%s): (%s)" % (codon_change_gt, codon_change))

    @data_provider(muc16_change_testdata)
    def test_muc16_change_protein(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, strand, transcript_change_gt, codon_change_gt, protein_change_gt):
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, fixed_id="muc16_change_protein", gene="MUC16")
        vcer = VariantClassifier()
        protein_change = vcer.generate_protein_change_from_vc(vc)
        self.assertTrue(protein_change == protein_change_gt, "Protein change did not match gt (%s): (%s)" % (protein_change_gt, protein_change))

    @data_provider(muc16testdata)
    def test_muc16_snps(self, chr, start, end, gt_vc, vt, ref, alt):
        """ Test all of the MUC16 SNPs."""
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt, fixed_id="muc16_snps", gene="MUC16")

    def test_snp_vc_on_one_transcript_5UTR(self):
        """Take the test transcript (ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand) and test the entire 5'UTR"""
        tx = self.retrieve_test_transcript_MAPK1(fixed_id="snp_vc_on_one_transcript_5UTR")
        vcer = VariantClassifier()

        chr = 22
        for i in range(0, 200):
            start = end = (22221919 - i)
            ref = tx.get_seq()[i]
            alt = 'G'
            vc = vcer.variant_classify(tx, ref, alt, start, end, "SNP").get_vc()
            if i < 189:
                self.assertTrue(vc == "5'UTR", "Should be 5'UTR, but saw " + vc + ".  For " + str([ref, alt, start, end]))
            if i < 200 and i >= 189:
                self.assertTrue(vc != "5'UTR", "Should not be 5'UTR, but saw " + vc + ".  For " + str([ref, alt, start, end]))

    def retrieve_test_transcript_MAPK1(self, fixed_id):
        ensembl_ds = self._create_ensembl_ds_from_testdata("MAPK1", fixed_id=fixed_id)
        tx = ensembl_ds.transcript_db['ENST00000215832.6']
        self.assertTrue(tx is not None, "Unit test appears to be misconfigured or a bug exists in the ensembl datasource code.")
        return tx

    variants_snps_splice_sites = lambda: (
        ("22", "22162138", "22162138", "Intron", "SNP", "T", "G"),
        ("22", "22162137", "22162137", "Splice_Site", "SNP", "T", "T"),
        ("22", "22162136", "22162136", "Splice_Site", "SNP", "C", "T"),
        ("22", "22162135", "22162135", "Splice_Site", "SNP", "G", "T"),
        ("22", "22162134", "22162134", "Splice_Site", "SNP", "A", "T"),
        ("22", "22162133", "22162133", "Missense_Mutation", "SNP", "G", "T"),
        ("22", "22162132", "22162132", "Silent", "SNP", "A", "T"),
        ("22", "22127164", "22127164", "Silent", "SNP", "C", "C"),
        ("22", "22127163", "22127163", "Splice_Site", "SNP", "T", "G"),
        ("22", "22127162", "22127162", "Splice_Site", "SNP", "C", "A"),
        ("22", "22127161", "22127161", "Splice_Site", "SNP", "C", "A"),
        ("22", "22127160", "22127160", "Splice_Site", "SNP", "A", "T"),
        ("22", "22127159", "22127159", "Intron", "SNP", "T", "A")
    )
    @data_provider(variants_snps_splice_sites)
    def test_snp_vc_on_one_transcript_splice_site(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test several positions on one MAPK1 transcript: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self.retrieve_test_transcript_MAPK1(fixed_id="snp_vc_on_one_transcript_splice_site")
        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref, alt, start, end, vt).get_vc()
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))

    # Stop codons: TAA, TAG, TGA
    # Reminder: this is a negative transcript and the stop codons listed on previous line are not reverse complemented
    variants_nonsense_nonstop = lambda: (
        ("22", "22123493", "22123493", VariantClassification.SILENT, "SNP", "T", "C"),
        ("22", "22123494", "22123494", VariantClassification.NONSTOP, "SNP", "T", "A"),
        ("22", "22123495", "22123495", VariantClassification.NONSTOP, "SNP", "A", "G"),
        ("22", "22123563", "22123563", VariantClassification.NONSENSE, "SNP", "A", "T")
    )
    @data_provider(variants_nonsense_nonstop)
    def test_vc_on_one_transcript_nonsense_nonstop(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test several positions on one MAPK1 transcript: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self.retrieve_test_transcript_MAPK1(fixed_id="vc_on_one_transcript_nonsense_nonstop")
        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref, alt, start, end, vt).get_vc()
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))


    variants_indels_splice_sites = lambda: (
        ("22", "22127155", "22127158", "Intron", "DEL", "CTCT", "-"),
        ("22", "22127155", "22127160", "Splice_Site", "DEL", "CTCTTA", "-"),
        ("22", "22127155", "22127163", "Splice_Site", "DEL", "CTCTTACCT", "-"),
        ("22", "22127155", "22127166", "Splice_Site", "DEL", "CTCTTACCTCGT", "-"),
        ("22", "22127163", "22127166", "Splice_Site", "DEL", "TCGT", "-"),
        ("22", "22127163", "22127166", "Splice_Site", "INS", "-", "AAAA"),
        ("22", "22127162", "22127166", "Splice_Site", "INS", "-", "AAAAA"),
        ("22", "22127167", "22127171", "Frame_Shift_Ins", "INS", "-", "AAAAA"),
        ("22", "22127155", "22127166", "Intron", "INS", "-", "AAAAAAAAAAAA"),
        ("22", "22221735", "22221735", "5'UTR", "INS", "-", "A"),
        ("22", "22123486", "22123486", "3'UTR", "INS", "-", "A"),
        ("22", "22123494", "22123494", "Stop_Codon_Del", "DEL", "A", "-"),
        ("22", "22221729", "22221729", "Start_Codon_Del", "DEL", "A", "-"),
        ("22", "22221735", "22221735", "5'UTR", "DEL", "G", "-"),
        # ref is wrong here, but test should still pass
        ("22", "22123486", "22123486", "3'UTR", "DEL", "A", "-")

        #TODO: Test RNA VCs
        #TODO: Test "+" strand transcript
        #TODO: Fix in frame calculation when only partially overlapping an exon.  (No need until secondary vc is implemented)
    )
    @data_provider(variants_indels_splice_sites)
    def test_indels_vc_on_one_transcript_splice_site(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test several positions on one MAPK1 transcript for indels: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self.retrieve_test_transcript_MAPK1(fixed_id="indels_vc_on_one_transcript_splice_site")
        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref, alt, start, end, vt).get_vc()
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))

    def test_determine_cds_in_exon_space(self):
        tx = self.retrieve_test_transcript_MAPK1(fixed_id="determine_cds_in_exon_space")
        vcer = VariantClassifier()
        cds_start, cds_stop = TranscriptProviderUtils.determine_cds_in_exon_space(tx)
        s,e = TranscriptProviderUtils.convert_genomic_space_to_exon_space(tx.get_start(), tx.get_end(), tx)
        self.assertTrue(s == 0, "Incorrect exon start: %d, gt: %d" % (s, 0))
        self.assertTrue(e == 11022, "Incorrect exon end: %d, gt: %d" % (e, 11022))
        self.assertTrue(cds_start == 189, "incorrect cds_start: %d, gt: %d" % (cds_start, 189))
        self.assertTrue(cds_stop == 1269, "incorrect cds_stop: %d, gt: %d" % (cds_stop, 1269))

    denovo_and_start_codon_test_data = lambda: (
        ("22", "22221729", "22221729", VariantClassification.START_CODON_SNP, "SNP", "A", "T"),
        ("22", "22221735", "22221737", VariantClassification.DE_NOVO_START_OUT_FRAME, "INS", "-", "CAT"),
        ("22", "22221735", "22221740", VariantClassification.DE_NOVO_START_OUT_FRAME, "INS", "-", "ACATAA"),
        ("22", "22221735", "22221741", VariantClassification.DE_NOVO_START_IN_FRAME, "INS", "-", "AACATAA"),
        ("22", "22221754", "22221754", VariantClassification.DE_NOVO_START_OUT_FRAME, "SNP", "C", "A"),
        ("22", "22221737", "22221737", VariantClassification.DE_NOVO_START_OUT_FRAME, "INS", "-", "A")
    )
    @data_provider(denovo_and_start_codon_test_data)
    def test_denovo_and_start_codon_on_one_transcript(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test several positions on one MAPK1 transcript for de novo and start codon: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self.retrieve_test_transcript_MAPK1(fixed_id="denovo_and_start_codon_on_one_transcript")
        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref, alt, start, end, vt).get_vc()
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))


    test_mutating_sequences = lambda: (
        ("AGGC", 0, 1, 1, "T", "SNP", "-", "AAGC"),
        ("AGGC", 0, 1, 1, "-", "DEL", "-", "AGC"),
        ("AGGC", 0, 1, 1, "-", "DEL", "+", "AGC"),
        ("AGGC", 0, 1, 3, "-", "DEL", "+", "A"),
        ("TAGGC", 0, 1, 3, "AAA", "INS", "+", "TAAAAGGC"),
        ("TTAGGC", 0, 1, 3, "AAA", "INS", "-", "TTTTTAGGC")

    )
    @data_provider(test_mutating_sequences)
    def test_mutating_sequences(self, seq_stranded, seq_index, exon_position_start, exon_position_end, alt_allele, variant_type, strand, gt):
        """Test that we can take a sequence, apply various mutations, and retrieve a correct sequence back."""
        observed_allele_stranded = alt_allele
        if strand == "-":
            observed_allele_stranded = Bio.Seq.reverse_complement(alt_allele)
        mutated_codon_seq = TranscriptProviderUtils.mutate_reference_sequence(seq_stranded, seq_index, exon_position_start, exon_position_end, observed_allele_stranded, variant_type)
        self.assertTrue(gt == mutated_codon_seq, "GT: " + gt + " Guess: " + mutated_codon_seq + "  " + str([seq_stranded, seq_index, exon_position_start, exon_position_end, alt_allele, variant_type, strand, gt]))

    #TODO: Test secondary VC
        #TODO: Test Flank (if not already done in MUC16 test)

if __name__ == '__main__':
    unittest.main()
