import unittest
import logging
import os.path as op

from oncotator.MutationData import MutationData
from oncotator.datasources.HgvsChangeTransformingDatasource import HgvsChangeTransformingDatasource
from TestUtils import TestUtils



TestUtils.setupLogging(__file__, __name__)
class HgvsChangeTransformingDatasourceTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()
        genecode_ds_path = op.join(self.config.get('DEFAULT', 'dbDir'), 'gencode_out2/hg19/gencode.v18.annotation.gtf')
        self.hgvs_datasource = HgvsChangeTransformingDatasource(genecode_ds_path)

    def test_annotate_SNP_missense(self):
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', 'Missense_Mutation')
        m.createAnnotation('genome_change', 'g.chr13:32914782C>T')
        m.createAnnotation('transcript_change', 'c.6290C>T')
        m.createAnnotation('protein_change', 'p.T2097M')
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'g.32914782C>T')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'c.6290C>T')
        self.assertEqual(m.annotations['HGVS_coding_protein_change'].getValue(), 'p.Thr2097Met')

    def test_annotate_SNP_nonsense(self):
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', 'Nonsense_Mutation')
        m.createAnnotation('genome_change', 'g.chr5:45303809G>A')
        m.createAnnotation('transcript_change', 'c.1510C>T')
        m.createAnnotation('protein_change', 'p.R504*')
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'g.45303809G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'c.1510C>T')
        self.assertEqual(m.annotations['HGVS_coding_protein_change'].getValue(), 'p.Arg504Ter')

    def test_annotate_SNP_intron(self):
        #+ strand transcript
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', 'Intron')
        m.createAnnotation('annotation_transcript', 'ENST00000402739.4')
        m.createAnnotation('transcript_strand', '+')
        m.createAnnotation('genome_change', 'g.chr2:80529551A>C')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'g.80529551A>C')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'c.1057-90785A>C')
        self.assertEqual(m.annotations['HGVS_coding_protein_change'].getValue(), '')

        #- strand transcript
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', 'Intron')
        m.createAnnotation('annotation_transcript', 'ENST00000356200.3')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr1:1636274C>T')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'g.1636274C>T')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'c.1428+69G>A')
        self.assertEqual(m.annotations['HGVS_coding_protein_change'].getValue(), '')

    def test_annotate_SNP_5_utr(self):
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', "5'UTR")
        m.createAnnotation('annotation_transcript', 'ENST00000402739.4')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr7:6865862G>C')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'g:6865862G>C')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'c.-34C>G')
        self.assertEqual(m.annotations['HGVS_coding_protein_change'].getValue(), '')

    def test_annotate_SNP_3_utr(self):
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', "3'UTR")
        m.createAnnotation('annotation_transcript', 'ENST00000521253.1')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('genome_change', 'g.chr8:27145409G>A')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'g:27145409G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'c.*220C>T')
        self.assertEqual(m.annotations['HGVS_coding_protein_change'].getValue(), '')

    def test_annotate_SNP_igr(self):
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', 'IGR')
        m.createAnnotation('chr', 15)
        m.createAnnotation('start', 30938316)
        m.createAnnotation('end', 30938316)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'A')
        m.createAnnotation('genome_change', '')
        m.createAnnotation('transcript_change', '')
        m.createAnnotation('protein_change', '')
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'g.30938316G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), '')
        self.assertEqual(m.annotations['HGVS_coding_protein_change'].getValue(), '')

    def test_annotate_SNP_silent(self):
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('variant_classification', 'Silent')
        m.createAnnotation('genome_change', 'g.chr1:19549914G>A')
        m.createAnnotation('transcript_change', 'c.2352C>T')
        m.createAnnotation('protein_change', 'p.I784I')
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'g.19549914G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'c.2352C>T')
        self.assertEqual(m.annotations['HGVS_coding_protein_change'].getValue(), '')

    def test_annotate_SNP_5_flank(self):
        pass

    def test_annotate_SNP_3_flank(self):
        pass

    def test_annotate_SNP_nonstop(self):
        pass

    def test_annotate_SNP_splice_site(self):
        pass

    def test_annotate_ONP_missense(self):
        #tests for genomic changes
        #tests for coding dna changes
        #tests for protein changes
        pass

    def test_annotate_ONP_nonsense(self):
        pass

    def test_annotate_ONP_intron(self):
        pass

    def test_annotate_ONP_5_utr(self):
        pass

    def test_annotate_ONP_3_utr(self):
        pass

    def test_annotate_ONP_igr(self):
        pass

    def test_annotate_ONP_5_flank(self):
        pass

    def test_annotate_ONP_3_flank(self):
        pass

    def test_annotate_ONP_nonstop(self):
        pass

    def test_annotate_ONP_silent(self):
        pass

    def test_annotate_ONP_splice_site(self):
        pass

    def test_annotate_INS_inframe(self):
        pass

    def test_annotate_INS_frameshift(self):
        pass

    def test_annotate_INS_start_codon(self):
        pass

    def test_annotate_INS_stop_codon(self):
        pass

    def test_annotate_DEL_inframe(self):
        pass

    def test_annotate_DEL_frameshift(self):
        pass

    def test_annotate_DEL_start_codon(self):
        pass

    def test_annotate_DEL_stop_codon(self):
        pass

