from functools import wraps
import shutil
import Bio

from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.index.gaf import region2bin, region2bins
from oncotator.utils.VariantClassification import VariantClassification
from test.MUC16ChangeTestdata import change_testdata
from test.PIK3CATestdata import pik3ca_testdata
from test.TestUtils import TestUtils,data_provider_decorator
from test.MUC16Testdata import muc16testdata
import unittest
from oncotator.utils.VariantClassifier import VariantClassifier

TestUtils.setupLogging(__file__, __name__)
class VariantClassifierTest(unittest.TestCase):

    # Do not uncomment this, unless you know what you are doing
    # _multiprocess_can_split_ = True

    GLOBAL_DS_BASENAME = "out/test_variant_classification_"

    def setUp(self):
        ensembl_ds = TestUtils._create_test_gencode_ds(VariantClassifierTest.GLOBAL_DS_BASENAME)
        self.ds = ensembl_ds

    def _create_ensembl_ds_trimmed(self):
        return self.ds

    def _determine_test_transcript(self, chr, start, end, alt, ref, vt):
        ensembl_ds = self._create_ensembl_ds_trimmed()
        recs = ensembl_ds.get_overlapping_transcripts(chr, start, end)
        txs = ensembl_ds._filter_transcripts(recs)
        tx = ensembl_ds._choose_transcript(txs, EnsemblTranscriptDatasource.TX_MODE_CANONICAL, vt, ref, alt, start, end)
        self.assertTrue(len(recs) != 0, "Issue with test...No transcripts found for: " + str([chr, start, end]))
        return tx

    def _test_variant_classification(self, alt, chr, end, gt_vc, ref, start, vt):
        tx = self._determine_test_transcript(chr, start, end, alt, ref, vt)

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
    @data_provider_decorator(variants_indels_MAPK1)
    def test_variant_classification_indels_simple(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test a small set of indels, from MAPK1.  These are all in frame or frame shifts cleanly within exons."""
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)

    variants_indels_MUC16 = lambda: (
        ("19", "8979051", "8979051", "Intron", "DEL", "A", "-"),
        ("19", "9058690", "9058690", "Frame_Shift_Del", "DEL", "T", "-"),
        ("19", "9064536", "9064536", "Frame_Shift_Del", "DEL", "T", "-"),
        ("19", "9006152", "9006153", "Intron", "INS", "-", "AA"),
        ("19", "9006225", "9006226", "Intron", "INS", "-", "CT"),
        ("19", "9073231", "9073232", "Frame_Shift_Ins", "INS", "-", "A")
    )
    @data_provider_decorator(variants_indels_MUC16)
    def test_variant_classification_indels_muc16(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test small set of indels, from MUC16.  These are all introns or frame shifts cleanly within exon /intron."""
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)

    variants_snps_missense = lambda: (
        ("22", "22127164", "22127164", "Missense_Mutation", "SNP", "C", "G"),
        ("22", "22143049", "22143049", "Nonsense_Mutation", "SNP", "C", "A"),
        ("22", "22162014", "22162014", "Missense_Mutation", "SNP", "C", "T")
    )
    @data_provider_decorator(variants_snps_missense)
    def test_variant_classification_missense(self, chr, start, end, gt_vc, vt, ref, alt):
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)

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

    @data_provider_decorator(frameshift_indels)
    def test_is_framshift_indel(self, vt, s, e, alt, gt):
        vcer = VariantClassifier()
        self.assertTrue(vcer.is_frameshift_indel(vt, s, e, alt) == gt)

    # ("MUC16", "19", "9057555", "9057555", "Missense_Mutation", "SNP", "C", "A", "g.chr19:9057555C>A", "-", "c.29891G>T", "c.(29890-29892)gGg>gTg", "p.G9964V"),
    @data_provider_decorator(change_testdata)
    def test_muc16_change_genome(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, strand, transcript_change_gt, codon_change_gt, protein_change_gt):
        """ Test all of the MUC16 changes (protein, genome, codon, and transcript)."""
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        genome_change = TranscriptProviderUtils.determine_genome_change(chr, start, end, ref, alt, vc.get_vt())
        self.assertTrue(genome_change == genome_change_gt, "Genome change did not match gt (%s): %s" %(genome_change_gt, genome_change))

    @data_provider_decorator(change_testdata)
    def test_muc16_change_transcript(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, strand, transcript_change_gt, codon_change_gt, protein_change_gt):
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        transcript_change = vcer.generate_transcript_change_from_tx(tx, vt, vc, int(start), int(end), ref, alt)
        self.assertTrue(transcript_change == transcript_change_gt, "Transcript change did not match gt (%s): %s  for %s" % (transcript_change_gt, transcript_change, str([chr, start, end, gt_vc, vt, ref, alt, vc.get_secondary_vc()])))

    @data_provider_decorator(change_testdata)
    def test_muc16_change_codon(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, strand, transcript_change_gt, codon_change_gt, protein_change_gt):
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        codon_change = vcer.generate_codon_change_from_vc(tx, int(start), int(end), vc)
        self.assertTrue(codon_change == codon_change_gt, "Codon change did not match gt (%s): (%s)" % (codon_change_gt, codon_change))

    @data_provider_decorator(change_testdata)
    def test_muc16_change_protein(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, strand, transcript_change_gt, codon_change_gt, protein_change_gt):
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        protein_change = vcer.generate_protein_change_from_vc(vc)
        self.assertTrue(protein_change == protein_change_gt, "Protein change did not match gt (%s): (%s) for %s" % (protein_change_gt, protein_change, str([chr, start, end, gt_vc, vt, ref, alt, vc.get_secondary_vc()])))

    @data_provider_decorator(muc16testdata)
    def test_muc16_snps(self, chr, start, end, gt_vc, vt, ref, alt):
        """ Test all of the MUC16 SNPs."""
        self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)

    @data_provider_decorator(pik3ca_testdata)
    def test_pik3ca_change_protein(self, gene, chr, start, end, gt_vc, ref, alt, vt, genome_change_gt, transcript_change_gt, protein_change_gt):
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        protein_change = vcer.generate_protein_change_from_vc(vc)
        self.assertTrue(protein_change == protein_change_gt, "Protein change did not match gt (%s): (%s) for %s" % (protein_change_gt, protein_change, str([chr, start, end, gt_vc, vt, ref, alt, vc.get_secondary_vc()])))

    @data_provider_decorator(pik3ca_testdata)
    def test_pik3ca_change_genome(self, gene, chr, start, end, gt_vc, ref, alt, vt, genome_change_gt, transcript_change_gt, protein_change_gt):
        """ Test all of the MUC16 changes (protein, genome, codon, and transcript)."""
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        genome_change = TranscriptProviderUtils.determine_genome_change(chr, start, end, ref, alt, vc.get_vt())
        self.assertTrue(genome_change == genome_change_gt, "Genome change did not match gt (%s): %s" %(genome_change_gt, genome_change))

    @data_provider_decorator(pik3ca_testdata)
    def test_pik3ca_change_transcript(self, gene, chr, start, end, gt_vc, ref, alt, vt, genome_change_gt, transcript_change_gt, protein_change_gt):
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        transcript_change = vcer.generate_transcript_change_from_tx(tx, vt, vc, int(start), int(end), ref, alt)
        self.assertTrue(transcript_change == transcript_change_gt, "Transcript change did not match gt (%s): %s  for %s" % (transcript_change_gt, transcript_change, str([chr, start, end, gt_vc, vt, ref, alt, vc.get_secondary_vc()])))

    variant_codons_to_check = lambda: (
        ("PIK3CA", "3", "178916619", "178916620", "Frame_Shift_Ins", "INS", "-", "CG", "g.chr3:178916619_178916620insCG", "c.6_7insCG", "c.(7-9)ccafs", "p.P3fs"),
        ("PIK3CA", "3", "178916619", "178916620", "In_Frame_Ins", "INS", "-", "CGA", "g.chr3:178916619_178916620insCGA", "c.6_7insCGA", "c.(7-9)cca>CGAcca", "p.2_3insR"),
        ("PIK3CA", "3", "178948154", "178948155", "Frame_Shift_Ins", "INS", "-", "GAATT", "g.chr3:178948154_178948155insGAATT", "c.2926_2930insGAATT", "c.(2926-2928)gaafs", "p.E976fs"),
        ("PIK3CA", "3", "178948155", "178948156", "Frame_Shift_Ins", "INS", "-", "GAATT", "g.chr3:178948155_178948156insGAATT", "c.2927_2931insGAATT", "c.(2926-2931)gaatttfs", "p.F977fs"),  # issue 109 is around this entry
        ("PIK3CA", "3", "178948159", "178948160", "Frame_Shift_Ins", "INS", "-", "GA",  "g.chr3:178948159_178948160insGA", "c.2931_2932insGA",  "c.(2932-2934)gagfs","p.E978fs"),
        ("PIK3CA", "3", "178916938", "178916940", "In_Frame_Del", "DEL", "GAA", "-", "g.chr3:178916938_178916940delGAA", "c.325_327delGAA", "c.(325-327)gaadel", "p.E110del"),
        ("PIK3CA", "3", "178948159", "178948160", "In_Frame_Ins", "INS", "-", "GAG",  "g.chr3:178948159_178948160insGAG", "c.2931_2932insGAG",  "c.(2932-2934)gag>GAGgag","p.978_978E>EE"),
        ("PIK3CA", "3", "178948160", "178948162", "In_Frame_Del", "DEL", "GAG", "-",  "g.chr3:178948160_178948162delGAG", "c.2932_2934delGAG",  "c.(2932-2934)gagdel", "p.E978del"),
        ("PIK3CA", "3", "178948160", "178948161", "Frame_Shift_Del", "DEL", "GA", "-",  "g.chr3:178948160_178948161delGA", "c.2932_2933delGA",  "c.(2932-2934)gagfs", "p.E978fs"),
        ("PIK3CA", "3", "178948160", "178948164", "Splice_Site", "DEL", "GAGAG", "-", "g.chr3:178948160_178948164delGAGAG", "c.2936_splice",  "c.e20+1", "p.ER978_splice"),
        ("PIK3CA", "3", "178948154", "178948158", "Frame_Shift_Del", "DEL", "GAATT", "-", "g.chr3:178948154_178948158delGAATT", "c.2926_2930delGAATT", "c.(2926-2931)gaatttfs", "p.EF976fs"),
        ("PIK3CA", "3", "178948154", "178948157", "Frame_Shift_Del", "DEL", "GAAT", "-", "g.chr3:178948154_178948158delGAAT", "c.2926_2929delGAAT", "c.(2926-2931)gaatttfs", "p.EF976fs"),
    )
    # chr3:178,916,611-178,916,632
    @data_provider_decorator(variant_codons_to_check)
    def test_pik3ca_change_codons_indels(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, transcript_change_gt, codon_change_gt, protein_change_gt):
        """Verify the codon change on a positive transcript for indels."""
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        codon_change = vcer.generate_codon_change_from_vc(tx, int(start), int(end), vc)
        self.assertTrue(codon_change == codon_change_gt, "Codon change did not match gt (%s): %s" %(codon_change_gt, codon_change))

    @data_provider_decorator(variant_codons_to_check)
    def test_pik3ca_change_proteins_indels(self, gene, chr, start, end, gt_vc, vt, ref, alt, genome_change_gt, transcript_change_gt, codon_change_gt, protein_change_gt):
        """Verify the protein change on a positive transcript for indels (issue 107)."""
        vc, tx = self._test_variant_classification(alt, chr, end, gt_vc, ref, start, vt)
        vcer = VariantClassifier()
        protein_change = vcer.generate_protein_change_from_vc(vc)
        self.assertTrue(protein_change == protein_change_gt, "Protein change did not match gt (%s): (%s) for %s" % (protein_change_gt, protein_change, str([chr, start, end, gt_vc, vt, ref, alt, vc.get_secondary_vc()])))


    def test_snp_vc_on_one_transcript_5UTR(self):
        """Take test transcript (ENST00000215832.6 (chr 22: 22108789:22221919) "-" strand) and test the entire 5'UTR"""
        tx = self._retrieve_test_transcript_MAPK1()
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

    def _retrieve_test_transcript(self, tx_id):
        ensembl_ds = self._create_ensembl_ds_trimmed()
        tx = ensembl_ds.transcript_db[tx_id]
        self.assertTrue(tx is not None,
                        "Unit test appears to be misconfigured or a bug exists in the ensembl datasource code.")
        return tx

    def _retrieve_test_transcript_MAPK1(self):
        tx_id = 'ENST00000215832.6'
        return self._retrieve_test_transcript(tx_id)

    def _retrieve_test_transcript_PIK3CA(self):
        tx_id = 'ENST00000263967.3'
        return self._retrieve_test_transcript(tx_id)

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
    @data_provider_decorator(variants_snps_splice_sites)
    def test_snp_vc_on_one_transcript_splice_site(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test several positions on one MAPK1 transcript: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self._retrieve_test_transcript_MAPK1()
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
    @data_provider_decorator(variants_nonsense_nonstop)
    def test_vc_on_one_transcript_nonsense_nonstop(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test several positions on one MAPK1 transcript: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self._retrieve_test_transcript_MAPK1()
        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref, alt, start, end, vt).get_vc()
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))


    variants_indels_splice_sites = lambda: (
        ("22", "22221735", "22221735", "5'UTR", "INS", "-", "A"),
        ("22", "22127155", "22127158", "Intron", "DEL", "CTCT", "-"),
        ("22", "22127155", "22127160", "Splice_Site", "DEL", "CTCTTA", "-"),
        ("22", "22127155", "22127163", "Splice_Site", "DEL", "CTCTTACCT", "-"),
        ("22", "22127155", "22127166", "Splice_Site", "DEL", "CTCTTACCTCGT", "-"),
        ("22", "22127163", "22127166", "Splice_Site", "DEL", "TCGT", "-"),
        ("22", "22127163", "22127166", "Splice_Site", "INS", "-", "AAAA"),
        ("22", "22127162", "22127166", "Splice_Site", "INS", "-", "AAAAA"),
        ("22", "22127167", "22127171", "Frame_Shift_Ins", "INS", "-", "AAAAA"),
        ("22", "22127155", "22127166", "Intron", "INS", "-", "AAAAAAAAAAAA"),
        ("22", "22123486", "22123486", "3'UTR", "INS", "-", "A"),
        ("22", "22123494", "22123494", "Stop_Codon_Del", "DEL", "A", "-"),
        ("22", "22221729", "22221729", "Start_Codon_Del", "DEL", "A", "-"),
        ("22", "22221735", "22221735", "5'UTR", "DEL", "G", "-"),
        # ref is wrong here, but test should still pass
        ("22", "22123486", "22123486", "3'UTR", "DEL", "A", "-")

        #TODO: Fix in frame calculation when only partially overlapping an exon.  (No need until secondary vc is implemented)
    )
    @data_provider_decorator(variants_indels_splice_sites)
    def test_indels_vc_on_one_transcript_splice_site(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test several positions on one MAPK1 transcript for indels: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self._retrieve_test_transcript_MAPK1()
        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref, alt, start, end, vt).get_vc()
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))

    def test_determine_cds_in_exon_space(self):
        tx = self._retrieve_test_transcript_MAPK1()
        vcer = VariantClassifier()
        cds_start, cds_stop = TranscriptProviderUtils.determine_cds_in_exon_space(tx)
        s,e = TranscriptProviderUtils.convert_genomic_space_to_exon_space(tx.get_start(), tx.get_end(), tx)
        self.assertTrue(s == 0, "Incorrect exon start: %d, gt: %d" % (s, 0))
        self.assertTrue(e == 11022, "Incorrect exon end: %d, gt: %d" % (e, 11022))
        self.assertTrue(cds_start == 189, "incorrect cds_start: %d, gt: %d" % (cds_start, 189))
        self.assertTrue(cds_stop == 1269, "incorrect cds_stop: %d, gt: %d" % (cds_stop, 1269))

    denovo_and_start_codon_test_data = lambda: (
        ("22", "22221904", "22221904", VariantClassification.DE_NOVO_START_OUT_FRAME, "DEL", "C", "-"),
        ("22", "22221794", "22221794", VariantClassification.DE_NOVO_START_OUT_FRAME, "DEL", "GC", "-"),
        ("22", "22221735", "22221740", VariantClassification.DE_NOVO_START_OUT_FRAME, "INS", "-", "ACATAA"),
        ("22", "22221735", "22221741", VariantClassification.DE_NOVO_START_IN_FRAME, "INS", "-", "AACATAA"),
        ("22", "22221729", "22221729", VariantClassification.START_CODON_SNP, "SNP", "A", "T"),
        ("22", "22221735", "22221737", VariantClassification.DE_NOVO_START_OUT_FRAME, "INS", "-", "CAT"),
        ("22", "22221754", "22221754", VariantClassification.DE_NOVO_START_OUT_FRAME, "SNP", "C", "A"),
        ("22", "22221737", "22221737", VariantClassification.DE_NOVO_START_OUT_FRAME, "INS", "-", "A")
    )
    @data_provider_decorator(denovo_and_start_codon_test_data)
    def test_denovo_and_start_codon_on_one_transcript(self, chr, start, end, gt_vc, vt, ref, alt):
        """Test several positions on one MAPK1 transcript for de novo and start codon: ENST00000215832.6 (chr 22: 22108789:22221919) Negative strand
        >ENST00000215832.6|ENSG00000100030.10|OTTHUMG00000030508.2|OTTHUMT00000075396.2|MAPK1-001|MAPK1|11022|UTR5:1-189|CDS:190-1272|UTR3:1273-11022|
        """
        tx = self._retrieve_test_transcript_MAPK1()
        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref, alt, start, end, vt).get_vc()
        self.assertTrue(gt_vc == vc, "Should have been " + gt_vc + ", but saw " + vc + "  with transcript " + tx.get_transcript_id() + " at " + str([chr, start, end, ref, alt]))


    transcript_ids = lambda: (
        # PIK3CA "+" chr 3
        ('ENST00000263967.3', ""),
        # MAPK1 "-"
        ('ENST00000215832.6', "")
    )
    @data_provider_decorator(transcript_ids)
    def test_basic_protein_position(self, tx_id, dummy):
        """ Test a 1-indexed protein position and splice sites are rendered properly for the entire test tx. """
        tx = self._retrieve_test_transcript(tx_id)
        strand = tx.get_strand()

        vcer = VariantClassifier()
        ctr = 0
        for j, cd in enumerate(tx.get_cds()):

            idx = 0
            if strand == "-":
                # First index since this is a negative strand.
                idx = 1

            start_0 = tx.get_cds()[j][idx]
            exon_0_range = tx.get_cds()[j][1] - tx.get_cds()[j][0]

            # Test against the codons in the cds
            for i in range(0, exon_0_range):
                if strand == "-":
                    # minus i, since this is a negative strand.
                    start = start_0 - i
                else:
                    start = start_0 + i + 1
                end = start

                protein_position_gt = 1 + ctr/3
                ctr += 1
                vt = VariantClassification.VT_SNP
                ref = "C"
                alt = "A"
                vc = vcer.variant_classify(tx, ref, alt, start, end, vt)
                prot_position_start = vc.get_ref_protein_start()
                self.assertTrue(vc.get_vc() != VariantClassification.INTRON and vc.get_secondary_vc() != VariantClassification.INTRON, "Intron should not have been seen here.  Failed on the %d position on the %d exon... %d .... %s" % (i, j, start, vc.get_vc()))
                self.assertTrue(prot_position_start == protein_position_gt, "Protein position failed on the %d position on the %d exon... (gt/guess) %d/%d  --- %d .... %s" % (i, j, protein_position_gt, prot_position_start, start, vc.get_vc()))

                # These need to be adjusted for a positive strand test
                is_at_last_exon = (j == (len(tx.get_cds())-1))
                is_at_first_exon = (j == 0)
                if strand == "-":
                    if start == (tx.get_cds()[j][0] + 1) or start == (tx.get_cds()[j][0] + 2):
                        self.assertTrue((vc.get_vc() == VariantClassification.SPLICE_SITE) or is_at_last_exon, "Not a splice site on the %d position on the %d far end of the exon... %d .... %s" % (i, j, start, vc.get_vc()))
                    if start == (tx.get_cds()[j][1] - 1) or start == (tx.get_cds()[j][1]):
                        self.assertTrue((vc.get_vc() == VariantClassification.SPLICE_SITE) or is_at_first_exon, "Not a splice site on the %d position on the %d close end of the exon... %d .... %s" % (i, j, start, vc.get_vc()))
                else:
                    if start == (tx.get_cds()[j][0] + 1) or start == (tx.get_cds()[j][0] + 2):
                        self.assertTrue((vc.get_vc() == VariantClassification.SPLICE_SITE) or is_at_first_exon, "Not a splice site on the %d position on the %d close end of the exon... %d .... %s" % (i, j, start, vc.get_vc()))
                    if start == (tx.get_cds()[j][1] - 1) or start == (tx.get_cds()[j][1]):
                        self.assertTrue((vc.get_vc() == VariantClassification.SPLICE_SITE) or is_at_last_exon, "Not a splice site on the %d position on the %d far end of the exon... %d .... %s" % (i, j, start, vc.get_vc()))

    mutating_sequences = lambda: (
        ("AGGC", 0, 1, 1, "T", "SNP", "-", "AAGC"),
        ("AGGC", 0, 1, 1, "-", "DEL", "-", "AGC"),
        ("AGGC", 0, 1, 1, "-", "DEL", "+", "AGC"),
        ("AGGC", 0, 1, 3, "-", "DEL", "+", "A"),
        ("TAGGC", 0, 1, 3, "AAA", "INS", "+", "TAAAAGGC"),
        ("TTAGGC", 0, 1, 3, "AAA", "INS", "-", "TTTTTAGGC")

    )
    @data_provider_decorator(mutating_sequences)
    def test_mutating_sequences(self, seq_stranded, seq_index, exon_position_start, exon_position_end, alt_allele, variant_type, strand, gt):
        """Test that we can take a sequence, apply various mutations, and retrieve a correct sequence back."""
        observed_allele_stranded = alt_allele
        if strand == "-":
            observed_allele_stranded = Bio.Seq.reverse_complement(alt_allele)
        mutated_codon_seq = TranscriptProviderUtils.mutate_reference_sequence(seq_stranded, seq_index, exon_position_start, exon_position_end, observed_allele_stranded, variant_type)
        self.assertTrue(gt == mutated_codon_seq, "GT: " + gt + " Guess: " + mutated_codon_seq + "  " + str([seq_stranded, seq_index, exon_position_start, exon_position_end, alt_allele, variant_type, strand, gt]))

    #TODO: Test secondary VC
    #TODO: Test Flank (if not already done in MUC16 test)

    mutating_exons = lambda: (
        ("DEL", "GG", "-", 22221734, 22221734, "AGAA"),
        ("DEL", "GGCT", "-", 22221734, 22221734, "GCAA"),
        ("SNP", "G", "T", 22221734, 22221734, "GCAAA"),
        ("SNP", "G", "C", 22221734, 22221734, "GCGAA"),
        ("DEL", "G", "-", 22221734, 22221734, "GCAA"),
        ("DEL", "GTTGGCT", "-", 22221731, 22221731, "GCAT"),
        ("INS", "-", "A", 22221733, 22221734, "CCTAA"),
        ("INS", "-", "GAG", 22221733, 22221734, "CCCTCAA"),
        ("INS", "-", "GAGA", 22221733, 22221734, "CCTCTCAA"),
        ("INS", "-", "GAGAAA", 22221733, 22221734, "CCTTTCTCAA")
    )

    @data_provider_decorator(mutating_exons)
    def test_mutate_exon_negative_stranded(self, vt, ref, alt, start, end, mutated_exon_gt):
        """Test that we can get the proper obs allele when mutating (negative strand). """
        tx = self._retrieve_test_transcript_MAPK1()
        vcer = VariantClassifier()
        exon_start, exon_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(start, end, tx)
        mutated_exon = vcer._mutate_exon(tx, ref, alt, vt, exon_start, buffer=2)
        self.assertTrue(mutated_exon == mutated_exon_gt, "GT/Guess: %s/%s" % (mutated_exon_gt, mutated_exon))

    # Essentially, the rendering of the reference sequence for inserts vs. non-insert
    reference_seq_testdata = lambda: (
        (3, 178948149, 178948150, "ACAAGA", [2], VariantClassification.VT_INS),
        (3, 178948150, 178948151, "AGA", [3], VariantClassification.VT_INS),
        (3, 178948148, 178948149, "ACA", [1], VariantClassification.VT_INS),
        (3, 178948147, 178948148, "ACA", [0], VariantClassification.VT_INS),
        (3, 178948146, 178948147, "AAGACA", [2], VariantClassification.VT_INS),
        (3, 178948150, 178948150, "ACA", [2], VariantClassification.VT_SNP),
        (3, 178948149, 178948150, "ACA", [1,2], VariantClassification.VT_SNP),
        (3, 178948149, 178948151, "ACAAGA", [1,2,0], VariantClassification.VT_SNP),
        (3, 178948149, 178948149, "ACA", [1], VariantClassification.VT_SNP),
        (3, 178948145, 178948145, "AAG", [0], VariantClassification.VT_SNP),
    )
    @data_provider_decorator(reference_seq_testdata)
    def test_reference_sequence_codon_construction_positive_strand(self, chr, start, end, ref_codon_sequence_gt, ref_codon_positions_gt, vt):
        tx = self._retrieve_test_transcript_PIK3CA()
        transcript_position_start, transcript_position_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(
            start, end, tx)
        if tx.get_strand() == "+" and not vt == VariantClassification.VT_INS:
            transcript_position_start -= 1
            transcript_position_end -= 1
        cds_start, cds_stop = TranscriptProviderUtils.determine_cds_in_exon_space(tx)
        protein_position_start, protein_position_end = TranscriptProviderUtils.get_protein_positions(
            transcript_position_start,
            transcript_position_end, cds_start)

        cds_codon_start, cds_codon_end = TranscriptProviderUtils.get_cds_codon_positions(protein_position_start, protein_position_end, cds_start)
        reference_codon_seq = tx.get_seq()[cds_codon_start:cds_codon_end+1]
        self.assertTrue(reference_codon_seq == ref_codon_sequence_gt, "Codon sequence did not match for reference (GT/Guess): %s/%s " % (ref_codon_sequence_gt, reference_codon_seq))
        # reference_codon_seq = TranscriptProviderUtils.mutate_reference_sequence(tx.get_seq()[cds_codon_start:cds_codon_end+1].lower(), cds_codon_start, transcript_position_start, transcript_position_end, reference_allele_stranded, variant_type)

    indel_testdata_for_change_pik3ca = lambda : (
        # In Frame
        ("3", "178916621","178916622", "In_Frame_Ins", "INS", "-", "TAT", "c.(7-12)ccacga>ccTATacga", "p.3_4PR>PIR"),
        ("3", "178916622","178916623", "In_Frame_Ins", "INS", "-", "TAT", "c.(10-12)cga>TATcga", "p.3_4insY"),
        ("3", "178916619","178916620", "In_Frame_Ins", "INS", "-", "TAT", "c.(7-9)cca>TATcca", "p.2_3insY"),
        ("3", "178916620","178916621", "In_Frame_Ins", "INS", "-", "TAT", "c.(7-9)cca>cTATca", "p.3_3P>LS"),
        ("3", "178916622","178916623", "In_Frame_Ins", "INS", "-", "CTTGAAGAA", "c.(10-12)cga>CTTGAAGAAcga", "p.3_4insLEE"),
        #fs
        ("3", "178916621","178916622", "Frame_Shift_Ins", "INS", "-", "TA", "c.(7-12)ccacgafs", "p.R4fs"),
        ("3", "178916622","178916623", "Frame_Shift_Ins", "INS", "-", "TA", "c.(10-12)cgafs", "p.R4fs"),
        ("3", "178916619","178916620", "Frame_Shift_Ins", "INS", "-", "TA", "c.(7-9)ccafs", "p.P3fs"),
        ("3", "178916620","178916621", "Frame_Shift_Ins", "INS", "-", "TA", "c.(7-9)ccafs", "p.P3fs"),
        ("3", "178916621","178916622", "Frame_Shift_Ins", "INS", "-", "T", "c.(7-12)ccacgafs", "p.R4fs"),
        ("3", "178916622","178916623", "Frame_Shift_Ins", "INS", "-", "T", "c.(10-12)cgafs", "p.R4fs"),
        ("3", "178916619","178916620", "Frame_Shift_Ins", "INS", "-", "T", "c.(7-9)ccafs", "p.P3fs"),
        ("3", "178916620","178916621", "Frame_Shift_Ins", "INS", "-", "T", "c.(7-9)ccafs", "p.P3fs"),
        ("3", "178916621","178916622", "Frame_Shift_Ins", "INS", "-", "TATT", "c.(7-12)ccacgafs", "p.R4fs"),
        ("3", "178916622","178916623", "Frame_Shift_Ins", "INS", "-", "TATT", "c.(10-12)cgafs", "p.R4fs"),
        ("3", "178916619","178916620", "Frame_Shift_Ins", "INS", "-", "TATT", "c.(7-9)ccafs", "p.P3fs"),
        ("3", "178916620","178916621", "Frame_Shift_Ins", "INS", "-", "TATT", "c.(7-9)ccafs", "p.P3fs"), # 17

        # DEL
        # In Frame
        ("3", "178916619", "178916621", "In_Frame_Del", "DEL", "TCC", "-", "c.(4-9)cctcca>cca", "p.2_3PP>P"),
        ("3", "178916620", "178916622", "In_Frame_Del", "DEL", "CCA", "-", "c.(7-9)ccadel", "p.P3del"),
        ("3", "178916621", "178916623", "In_Frame_Del", "DEL", "CAC", "-", "c.(7-12)ccacga>cga", "p.P3del"),
        ("3", "178916622", "178916624", "In_Frame_Del", "DEL", "ACG", "-", "c.(7-12)ccacga>cca", "p.R4del"), #21

        #fs
        ("3", "178916619","178916620", "Frame_Shift_Del", "DEL", "TC","-", "c.(4-9)cctccafs", "p.PP2fs"),
        ("3", "178916620","178916621", "Frame_Shift_Del", "DEL", "CC","-", "c.(7-9)ccafs", "p.P3fs"),
        ("3", "178916621","178916622", "Frame_Shift_Del", "DEL", "CA","-", "c.(7-9)ccafs", "p.P3fs"),
        ("3", "178916622","178916623", "Frame_Shift_Del", "DEL", "AC","-", "c.(7-12)ccacgafs", "p.R4fs"), # This is correct, since the first amino acid remains the same
        ("3", "178916619","178916619", "Frame_Shift_Del", "DEL", "T","-", "c.(4-6)cctfs", "p.P3fs"), #26 -- correct as written, since there are two P's in a row
        ("3", "178916620","178916620", "Frame_Shift_Del", "DEL", "C","-", "c.(7-9)ccafs", "p.P3fs"),
        ("3", "178916621","178916621", "Frame_Shift_Del", "DEL", "C","-", "c.(7-9)ccafs", "p.P3fs"), #28
        ("3", "178916622","178916622", "Frame_Shift_Del", "DEL", "A","-", "c.(7-9)ccafs", "p.P3fs"), #29

        ("3", "178916619","178916622", "Frame_Shift_Del", "DEL", "TCCA","-", "c.(4-9)cctccafs", "p.PP2fs"),
        ("3", "178916620","178916623", "Frame_Shift_Del", "DEL", "CCAC","-", "c.(7-12)ccacgafs", "p.PR3fs"),
        ("3", "178916621","178916624", "Frame_Shift_Del", "DEL", "CACG","-", "c.(7-12)ccacgafs", "p.PR3fs"),
        ("3", "178916622","178916625", "Frame_Shift_Del", "DEL", "ACGA","-", "c.(7-12)ccacgafs", "p.PR3fs"),

        ("3", "178916619","178916625", "Frame_Shift_Del", "DEL", "TCCACGA","-", "c.(4-12)cctccacgafs", "p.PPR2fs"),
        ("3", "178916620","178916626", "Frame_Shift_Del", "DEL", "CCACGAC","-", "c.(7-15)ccacgaccafs", "p.PRP3fs"),
        ("3", "178916621","178916627", "Frame_Shift_Del", "DEL", "CACGACC","-", "c.(7-15)ccacgaccafs", "p.PRP3fs"), #36
        ("3", "178916622","178916628", "Frame_Shift_Del", "DEL", "ACGACCA","-", "c.(7-15)ccacgaccafs", "p.PRP3fs"),
    )
    @data_provider_decorator(indel_testdata_for_change_pik3ca)
    def test_reference_change_construction_positive_strand(self, chr, start, end, vc_gt, vt, ref_allele, alt_allele, codon_change_gt, protein_change_gt):
        """Test that different indel configurations are rendered correctly on a positive stranded transcript."""
        tx = self._retrieve_test_transcript_PIK3CA()

        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref_allele, alt_allele, start, end, vt, dist=2)
        protein_change = vcer.generate_protein_change_from_vc(vc)
        codon_change = vcer.generate_codon_change_from_vc(tx, int(start), int(end), vc)
        self.assertTrue(codon_change == codon_change_gt, "GT/Guess: %s/%s" % (codon_change_gt, codon_change))
        self.assertTrue(protein_change == protein_change_gt, "GT/Guess: %s/%s" % (protein_change_gt, protein_change))

    indel_testdata_for_change_mapk1 = lambda : (
        # In Frame
        ("22", "22221703","22221704", "In_Frame_Ins", "INS", "-", "TAT", "c.(25-30)gcgggc>gcgATAggc", "p.9_10AG>AIG"),  # gcccgc
        ("22", "22221702","22221703", "In_Frame_Ins", "INS", "-", "TAT", "c.(28-30)ggc>gATAgc", "p.10_10G>DS"),
        ("22", "22221701","22221702", "In_Frame_Ins", "INS", "-", "TAT", "c.(28-30)ggc>ggATAc", "p.10_11insY"),  # gcccgc -- this is correct since G>GY, we are just inserting a Y
        ("22", "22221700","22221701", "In_Frame_Ins", "INS", "-", "TAT", "c.(28-33)ggcccg>ggcATAccg", "p.10_11GP>GIP"),
        ("22", "22221700","22221701", "In_Frame_Ins", "INS", "-", "TTCTTCAAG", "c.(28-33)ggcccg>ggcCTTGAAGAAccg", "p.10_11GP>GLEEP"),

        #fs
        ("22", "22221724","22221725", "Frame_Shift_Ins", "INS", "-", "TATT", "c.(4-9)gcggcgfs", "p.A3fs"),
        ("22", "22221723","22221724", "Frame_Shift_Ins", "INS", "-", "TATT", "c.(7-9)gcgfs", "p.A3fs"),
        ("22", "22221722","22221723", "Frame_Shift_Ins", "INS", "-", "TATT", "c.(7-9)gcgfs", "p.-3fs"), # This is correct, since the protein does not change, but is fs # CGC > CTATTGC  gcg>gATAAcg
        ("22", "22221721","22221722", "Frame_Shift_Ins", "INS", "-", "TATT", "c.(7-12)gcggcgfs", "p.A4fs"),
        ("22", "22221724","22221725", "Frame_Shift_Ins", "INS", "-", "TA", "c.(4-9)gcggcgfs", "p.A3fs"),
        ("22", "22221723","22221724", "Frame_Shift_Ins", "INS", "-", "TA", "c.(7-9)gcgfs", "p.A3fs"),
        ("22", "22221722","22221723", "Frame_Shift_Ins", "INS", "-", "TA", "c.(7-9)gcgfs", "p.A3fs"),
        ("22", "22221721","22221722", "Frame_Shift_Ins", "INS", "-", "TA", "c.(7-12)gcggcgfs", "p.A4fs"),
        ("22", "22221724","22221725", "Frame_Shift_Ins", "INS", "-", "T", "c.(4-9)gcggcgfs", "p.A3fs"), #14
        ("22", "22221723","22221724", "Frame_Shift_Ins", "INS", "-", "T", "c.(7-9)gcgfs", "p.A3fs"),
        ("22", "22221722","22221723", "Frame_Shift_Ins", "INS", "-", "T", "c.(7-9)gcgfs", "p.A3fs"),
        ("22", "22221721","22221722", "Frame_Shift_Ins", "INS", "-", "T", "c.(7-12)gcggcgfs", "p.A4fs"),

        # DEL
        # In Frame
        ("22", "22221720", "22221722", "In_Frame_Del", "DEL", "GCC", "-", "c.(7-12)gcggcg>gcg", "p.3_4AA>A"), #18
        ("22", "22221719", "22221721", "In_Frame_Del", "DEL", "CGC", "-", "c.(10-12)gcgdel", "p.A7del"),  # This is technically correct, since there are A's upstream.
        ("22", "22221718", "22221720", "In_Frame_Del", "DEL", "CCG", "-", "c.(10-15)gcggcg>gcg", "p.4_5AA>A"),
        ("22", "22221717", "22221719", "In_Frame_Del", "DEL", "GCC", "-", "c.(10-15)gcggcg>gcg", "p.4_5AA>A"), #21

        #fs
        ("22", "22221700","22221701", "Frame_Shift_Del", "DEL", "GG","-", "c.(28-33)ggcccgfs", "p.P11fs"),
        ("22", "22221699","22221700", "Frame_Shift_Del", "DEL", "GG","-", "c.(31-33)ccgfs", "p.P11fs"),
        ("22", "22221698","22221699", "Frame_Shift_Del", "DEL", "CG","-", "c.(31-33)ccgfs", "p.P11fs"),
        ("22", "22221697","22221698", "Frame_Shift_Del", "DEL", "CC","-", "c.(31-36)ccggagfs", "p.E12fs"),
        ("22", "22221700","22221700", "Frame_Shift_Del", "DEL", "G","-", "c.(31-33)ccgfs", "p.P11fs"), #26
        ("22", "22221699","22221699", "Frame_Shift_Del", "DEL", "G","-", "c.(31-33)ccgfs", "p.P11fs"),
        ("22", "22221698","22221698", "Frame_Shift_Del", "DEL", "C","-", "c.(31-33)ccgfs", "p.P11fs"), #28
        ("22", "22221697","22221697", "Frame_Shift_Del", "DEL", "C","-", "c.(34-36)gagfs", "p.E12fs"), #29

        ("22", "22221694","22221697", "Frame_Shift_Del", "DEL", "TCTC","-", "c.(34-39)gagatgfs", "p.EM12fs"),
        ("22", "22221693","22221696", "Frame_Shift_Del", "DEL", "ATCT","-", "c.(34-39)gagatgfs", "p.EM12fs"),
        ("22", "22221692","22221695", "Frame_Shift_Del", "DEL", "CATC","-", "c.(34-39)gagatgfs", "p.EM12fs"),
        ("22", "22221691","22221694", "Frame_Shift_Del", "DEL", "CCAT","-", "c.(37-42)atggtcfs", "p.MV13fs"),

        ("22", "22221690","22221696", "Frame_Shift_Del", "DEL", "ACCATCT","-", "c.(34-42)gagatggtcfs", "p.EMV12fs"),
        ("22", "22221689","22221695", "Frame_Shift_Del", "DEL", "GACCATC","-", "c.(34-42)gagatggtcfs", "p.EMV12fs"),
        ("22", "22221688","22221694", "Frame_Shift_Del", "DEL", "GGACCAT","-", "c.(37-45)atggtccgcfs", "p.MVR13fs"), #36
        ("22", "22221687","22221693", "Frame_Shift_Del", "DEL", "CGGACCA","-", "c.(37-45)atggtccgcfs", "p.MVR13fs"),
    )
    @data_provider_decorator(indel_testdata_for_change_mapk1)
    def test_reference_change_construction_negative_strand(self, chr, start, end, vc_gt, vt, ref_allele, alt_allele, codon_change_gt, protein_change_gt):
        """Test that different indel configurations are rendered correctly on a negative stranded transcript."""
        tx = self._retrieve_test_transcript_MAPK1()

        vcer = VariantClassifier()
        vc = vcer.variant_classify(tx, ref_allele, alt_allele, start, end, vt, dist=2)
        protein_change = vcer.generate_protein_change_from_vc(vc)
        codon_change = vcer.generate_codon_change_from_vc(tx, int(start), int(end), vc)
        self.assertTrue(codon_change == codon_change_gt, "GT/Guess: %s/%s" % (codon_change_gt, codon_change))
        self.assertTrue(protein_change == protein_change_gt, "GT/Guess: %s/%s" % (protein_change_gt, protein_change))

if __name__ == '__main__':
    unittest.main()
