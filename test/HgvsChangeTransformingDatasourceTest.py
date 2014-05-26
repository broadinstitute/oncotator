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

    ### TODO need test to assert that all necessary fields are present

    def test_annotate_SNP_missense(self):
        #rs80358866
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr13.hg19:g.32914782C>T')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000380152.3:c.6290C>T')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000369497:p.Thr2097Met')

    def test_annotate_SNP_nonsense(self):
        #rs35229491
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr5.hg19:g.45303809G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000303230.4:c.1510C>T')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000307342:p.Arg504*')

    def test_annotate_renders_with_no_build(self):
        #rs148119501
        """If mutation instance being annotated does not have a build value or is '', annotate should
        return a genome_change value with just chr. i.e. chr2:g.80529551A>C vs. chr2.hg19:g.80529551A>C"""
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr2:g.80529551A>C')

    def test_annotate_SNP_intron(self):
        #rs148119501
        #+ strand transcript
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr2.hg19:g.80529551A>C')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000402739.4:c.1057-90785A>C')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

        #- strand transcript
        #rs78420771
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr10.hg19:g.118891993A>G')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000277905.2:c.430-5T>C')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_SNP_5_utr(self):
        #rs141173433
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr7.hg19:g.6865862G>C')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000316731.8:c.-34C>G')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_SNP_3_utr(self):
        #rs143436239
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr8.hg19:g.27145409G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000521253.1:c.*220C>T')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_SNP_igr(self):
        #rs112615235
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr15.hg19:g.30938316G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), '')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_SNP_silent(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.19549914G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000477853.1:c.2352C>T')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_SNP_splice_site(self):
        #splice site mutation occuring in intron prior to coding start position
        #rs61191258
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr19.hg19:g.52994576G>A')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000421239.2:c.-121-1G>A')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

        #splice site mutation occuring in intron after coding start position
        #rs144524702
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr5.hg19:g.484634C>T')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000264938.3:c.932+1G>A')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_SNP_de_novo_start_OutOfFrame(self):
        #rs114472931
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.45140082G>T')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000372237.3:c.-19C>A')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_ONP_missense(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr22.hg19:g.27003913_27003914delinsAT')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000215939.2:c.371_372delinsAT')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000215939:p.Trp124Tyr')

    def test_annotate_INS_inframe_1(self):
        m = MutationData()
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
        
        m = self.hgvs_datasource.annotate_mutation(m)

        #this ins of GCC occurs in a GCC-repeat region and thus need to 3' adjust position for HGVS compliance
        # it is technically a duplication
        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr5.hg19:g.113698641_113698643dupGCC')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000512097.3:c.169_171dupGCC')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000427120:p.Ala58dup')

    def test_annotate_INS_inframe_2(self):
        m = MutationData()
        m.createAnnotation('variant_type', 'INS')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('chr', '7')
        m.createAnnotation('start', 11871469)
        m.createAnnotation('end', 11871470)
        m.createAnnotation('ref_allele', '-')
        m.createAnnotation('alt_allele', 'GCAGCG')
        m.createAnnotation('transcript_strand', '-')
        m.createAnnotation('variant_classification', 'In_Frame_Ins')
        m.createAnnotation('annotation_transcript', 'ENST00000423059.3')
        m.createAnnotation('genome_change', 'g.chr7:11871469_11871470insGCAGCG')
        m.createAnnotation('transcript_change', 'c.103_104insCGCTGC')
        m.createAnnotation('protein_change', 'p.34_35insPL')
        #m.createAnnotation('ref_context', 'cagcagcaggagcagcggcagc')
        m.createAnnotation('ref_context', 'CGCAGCCCTGCCGGCGCCCGGGCGTAGCAGCAGCAGCAGGAGCAGCGGCAGCGGCAGCGGCAGCGGCAGCAGCTGCAGGACG') # need a larger ref_context size to get the correct mapping
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr7.hg19:g.11871488_11871493dupGCAGCG')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000423059.3:c.98_103dupCGCTGC')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000406482:p.Pro33_Leu34dup')

    def test_annotate_INS_inframe_3(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr8.hg19:g.10467629_10467630insTTC')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000382483.3:c.3978_3979insGAA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000371923:p.Lys1326_Thr1327insGlu')

    def test_annotate_INS_inframe_4(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr8.hg19:g.10467628_10467629insCCC')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000382483.3:c.3979_3980insGGG')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000371923:p.Thr1327delinsArgAla')

    def test_annotate_INS_inframe_5(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr2.hg19:g.3197915_3197917dupCAT')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000398659.4:c.755_757dupATG')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000381652:p.Asn252_Ala253insAsp')

    def test_annotate_INS_inframe_6(self):
        #This is an insertion of a STOP in between two amino acids
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.248637602_248637603insTGA')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000359594.2:c.951_952insTGA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000352604:p.Lys318*')

    def test_annotate_INS_inframe_7(self):
        #This is an insertion of a STOP right before a stop
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.248637605_248637606insTGA')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000359594.2:c.954_955insTGA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_INS_frameshift(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        #this ins of CG does NOT occurs next to a CG and does not need to be position adjusted
        # it is technically an insertion
        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr4.hg19:g.1388441_1388442insCG')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000324803.4:c.142_143insCG')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000323978:p.Met48fs')

    def test_annotate_INS_frameshift_2(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr9.hg19:g.135977872_135977875dupCGCT')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000393160.3:c.1832_1835dupAGCG')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000376867:p.Arg612fs')

    def test_annotate_DEL_inframe(self):
        #rs141326765
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        #this deletion is straightforward, no position adjustments necessary
        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr14.hg19:g.70924869_70924871delATG')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000603540.1:c.653_655delATG')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000474385:p.Asp219del')

    def test_annotate_DEL_inframe_2(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        #this deletion is straightforward, no position adjustments necessary
        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr12.hg19:g.50156659_50156667delAAGAAGAAA')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000552699.1:c.868_876delAAGAAGAAA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000446734:p.Lys290_Lys292del')

    def test_annotate_DEL_inframe_3(self):
        #rs141326765
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        #this deletion is straightforward, no position adjustments necessary
        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr19.hg19:g.40900189_40900191delTCC')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000324001.7:c.4077_4079delGGA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000326018:p.Glu1361del')

    def test_annotate_DEL_frameshift(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        #Here only the genomic change needs to get '3 shifted because the transcript is negative strand 
        #and the coding postion is already the most 3'
        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr19.hg19:g.11348960delG')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000294618.7:c.1664delC')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000294618:p.Pro555fs')

    def test_annotate_SNP_nonstop(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr7.hg19:g.55273310A>G')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000275493.2:c.3633A>G')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000275493:p.*1211Trpext*6') #6 new amino acids added until another stop codon is encountered
        # "p.*1211Trpext?" would describe a variant in the stop codon at position 1211 changing it to a codon for Tryptophan (Trp, W) and adding a tail of new amino acids of unknown length since the shifted frame does not contain a new stop codon.

    def test_annotate_stop_codon_DEL_1(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        #Here only the genomic change needs to get '3 shifted because the transcript is negative strand 
        #and the coding postion is already the most 3'
        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.248637607delA')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000359594.2:c.956delA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000352604:p.*319Cysext*?')

    def test_annotate_stop_codon_DEL_2(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        #Here only the genomic change needs to get '3 shifted because the transcript is negative strand 
        #and the coding postion is already the most 3'
        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.248637605_248637606delGT')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000359594.2:c.954_955delGT')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000352604:p.*319Valext*?')

    def test_annotate_stop_codon_DEL_3(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.248637608_248637610delGTA')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000359594.2:c.957_*2delGTA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), '')

    def test_annotate_stop_codon_DEL_4(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.248637602_248637610delAAAGTAGTA')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000359594.2:c.951_*2delAAAGTAGTA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000352604:p.Glu317_*319delinsGluext*?')

    def test_annotate_stop_codon_DEL_5(self):
        #negative strand transcript
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr2.hg19:g.29416090_29416091delTC')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000389048.3:c.4862_4863delGA')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000373700:p.*1621Cysext*53')

    def test_annotate_stop_codon_DEL_6(self):
        #negative strand transcript
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr2.hg19:g.29416092_29416094delAGG')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000389048.3:c.4859_4861delCCT')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000373700:p.Pro1620_*1621delinsArgext*41')

    def test_annotate_stop_codon_DEL_7(self):
        #negative strand transcript
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr2.hg19:g.29416085_29416090delCGAGCT')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000389048.3:c.4863_*5delAGCTCG')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000373700:p.*1621Trpext*39')

    def test_annotate_stop_codon_INS(self):
        m = MutationData()
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
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.248637606_248637607insCAT')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000359594.2:c.955_956insCAT')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000352604:p.*319Serext*1')

    def test_annotate_stop_codon_ONP(self):
        m = MutationData()
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
        
        m = self.hgvs_datasource.annotate_mutation(m)

        self.assertEqual(m.annotations['HGVS_genomic_change'].getValue(), 'chr1.hg19:g.248637605_248637606delinsCC')
        self.assertEqual(m.annotations['HGVS_coding_DNA_change'].getValue(), 'ENST00000359594.2:c.954_955delinsCC')
        self.assertEqual(m.annotations['HGVS_protein_change'].getValue(), 'ENSP00000352604:p.Lys318_*319delinsAsnGlnext*1')

