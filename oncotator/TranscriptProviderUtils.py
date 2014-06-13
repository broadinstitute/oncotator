import math
import Bio
from oncotator.utils.VariantClassification import VariantClassification


class TranscriptProviderUtils(object):

    @staticmethod
    def infer_variant_type(reference_allele, observed_allele):
        """ To go completely in annotate method.  Returns variant type string."""
        if reference_allele == '-': #is insertion
            return 'INS'
        elif observed_allele == '-': #is deletion
            return 'DEL'
        else:
            if len(observed_allele) > len(reference_allele):
                return 'INS'
            elif len(observed_allele) < len(reference_allele):
                return 'DEL'
            elif len(reference_allele) == 1:
                return 'SNP'
            elif len(reference_allele) == 2:
                return 'DNP'
            elif len(reference_allele) == 3:
                return 'TNP'
            elif len(reference_allele) > 3:
                return 'ONP'

        raise Exception('Variant Type cannot be inferred from reference and observed allele: (%s, %s)' % (reference_allele, observed_allele))

    @staticmethod
    def retrieve_effect_dict():
        """Used for ordering the effect given the variant_classification"""
        return {
        'De_novo_Start_OutOfFrame':0,
        'Nonsense_Mutation':0,
        'Nonstop_Mutation':0,
        'Missense_Mutation':1,
        'De_novo_Start_InFrame':1,
        'In_Frame_Del':1,
        'In_Frame_Ins':1,
        'Frame_Shift_Del':2,
        'Frame_Shift_Ins':2,
        'Frame_Shift_Sub':2,
        'Start_Codon_SNP':3,
        'Start_Codon_Del':3,
        'Start_Codon_Ins':3,
        'Start_Codon_DNP':3,
        'Start_Codon_TNP':3,
        'Start_Codon_ONP':3,
        'Stop_Codon_SNP':3,
        'Stop_Codon_Del':3,
        'Stop_Codon_Ins':3,
        'Stop_Codon_DNP':3,
        'Stop_Codon_TNP':3,
        'Stop_Codon_ONP':3,
        'Splice_Site':4,
        'Splice_Site_SNP':4,
        'Splice_Site_Del':4,
        'Splice_Site_Ins':4,
        'Splice_Site_DNP':4,
        'Splice_Site_TNP':4,
        'Splice_Site_ONP':4,
        'Splice_Site':4,
        'miRNA':4,
        'Silent':5,
        "3'UTR":6,
        "5'UTR":6,
        'Intron':7,
        "5'Flank":8,
        "3'Flank":8,
        'Non-coding_Transcript':9,
        'IGR':20,
        'TX-REF-MISMATCH':100
        }

    @staticmethod
    def test_overlap(a_st, a_en, b_st, b_en):
        """Test if two genomic start, end tuples overlap.  """
        if (a_st >= b_st and a_st <= b_en) or (a_en >= b_st and a_en <= b_en) or \
            (a_st <= b_st and a_en >= b_en):
            return True
        else:
            return False

    @staticmethod
    def test_feature_overlap(a_st, a_en, tuple_list):
        result = -1
        for i,f in enumerate(tuple_list):
            # Must add one to keep convention in line with genome browser
            b_st = f[0]+1
            b_en = f[1]
            is_overlap = TranscriptProviderUtils.test_overlap(a_st, a_en, b_st, b_en)
            if is_overlap:
                return i
        return result

    @staticmethod
    def test_overlap_with_strand(a_st, a_en, b_st, b_en, strand):
        if strand == '+':
            if a_st <= b_st and a_en >= b_en: return 'a_encompasses_b'
            elif a_st >= b_st and a_en <= b_en: return 'a_within_b'
            elif a_st < b_st and a_en >= b_st: return 'a_overlaps_b_left_border'
            elif a_st <= b_en and a_en >= b_en: return 'a_overlaps_b_right_border'
            else: return None
        elif strand == '-':
            b_st, b_en = b_en, b_st
            if a_st <= b_st and a_en >= b_en: return 'a_encompasses_b'
            elif a_st >= b_st and a_en <= b_en: return 'a_within_b'
            elif a_st < b_st and a_en >= b_st: return 'a_overlaps_b_right_border'
            elif a_st <= b_en and a_en >= b_en: return 'a_overlaps_b_left_border'
            else: return None

    @staticmethod
    def determine_genome_change(chr, start, end, ref_allele, alt_allele, variant_type):
        genome_change = ''
        start = int(start)
        end = int(end)
        if variant_type == 'SNP':
            genome_change = 'g.chr%s:%d%s>%s' % (chr, start, ref_allele,
                                                 alt_allele)
        elif variant_type.endswith('NP'):
            genome_change = 'g.chr%s:%d_%d%s>%s' % (chr, start, end,
                                                    ref_allele, alt_allele)
        elif variant_type == 'DEL':
            if start == end:
                genome_change = 'g.chr%s:%ddel%s' % (chr, start, ref_allele)
            else:
                genome_change = 'g.chr%s:%d_%ddel%s' % (chr, start, end, ref_allele)
        elif variant_type == 'INS':
            genome_change = 'g.chr%s:%d_%dins%s' % (chr, start, end,
                                                    alt_allele)
        return genome_change

    @staticmethod
    def render_transcript_change(variant_type, variant_classification, exon_position_start, exon_position_end, ref_allele_stranded, alt_allele_stranded):
        """

        :param variant_type:
        :param variant_classification:
        :param exon_position_start: Coordinates in transcript/exon space
        :param exon_position_end: Coordinates in transcript/exon space
        :param ref_allele_stranded: ref_allele with strand already accounted
        :param alt_allele_stranded: alt_allele with strand already accounted
        """
        transcript_change = ""
        if variant_classification.startswith('Splice_Site'):
            return 'c.%d_splice' % (exon_position_start)

        if variant_type == VariantClassification.VT_SNP:
            transcript_change = 'c.%d%s>%s' % (exon_position_start, ref_allele_stranded, alt_allele_stranded)
        elif variant_type.endswith('NP'):
            transcript_change = 'c.%d_%d%s>%s' % (exon_position_start,
                exon_position_end, ref_allele_stranded, alt_allele_stranded)
        elif variant_type == VariantClassification.VT_DEL:
            if exon_position_start == exon_position_end:
                transcript_change = 'c.%ddel%s' % (exon_position_start, ref_allele_stranded)
            else:
                transcript_change = 'c.%d_%ddel%s' % (exon_position_start,
                    exon_position_end, ref_allele_stranded)
        elif variant_type == VariantClassification.VT_INS:
            transcript_change = 'c.%d_%dins%s' % (exon_position_start,
                exon_position_end, alt_allele_stranded)

        return transcript_change

    @staticmethod
    def render_transcript_position(genomic_start, genomic_end, tx):
        """
        Create a transcript position that looks like [exon_start]_[exon_end] or [exon_position] (when start == end)
        :param genomic_start: int
        :param genomic_end: int
        :return: a string representation
        """
        exon_start, exon_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(genomic_start, genomic_end, tx)
        if exon_start == exon_end:
            return "%d" % exon_start
        else:
            return "%d_%d" % (exon_start, exon_end)


    @staticmethod
    def _render_basic_protein_change(variant_classification, ref_prot_allele, alt_prot_allele, prot_position_end,
                                     prot_position_start):
        """
        Helper method that handles the main variant classifications: INDELS and SNPs... no splice sites.
        """
        if variant_classification.startswith('Frame_Shift'):
            protein_change = 'p.%s%dfs' % (ref_prot_allele, prot_position_start)
        elif alt_prot_allele == '-' or alt_prot_allele == '':
            protein_change = 'p.%s%ddel' % (ref_prot_allele, prot_position_start)
        elif ref_prot_allele == '-' or ref_prot_allele == '':
            protein_change = 'p.%d_%dins%s' % (prot_position_start,
                                               prot_position_end, alt_prot_allele)
        elif len(ref_prot_allele) == 1 and len(alt_prot_allele) == 1:
            protein_change = 'p.%s%d%s' % (ref_prot_allele, prot_position_start, alt_prot_allele)
        else:
            protein_change = 'p.%d_%d%s>%s' % (prot_position_start, prot_position_end, ref_prot_allele, alt_prot_allele)
        return protein_change

    @staticmethod
    def render_protein_change(variant_type, variant_classification, prot_position_start, prot_position_end, ref_prot_allele, alt_prot_allele, secondary_variant_classification=""):
        if variant_classification.startswith('Splice_Site'):
            # True when exon side splice site.
            if prot_position_start > 0:
                # protein_change = 'p.%s%d_splice' % (ref_prot_allele, prot_position_start)
                protein_change = TranscriptProviderUtils._render_basic_protein_change(secondary_variant_classification, ref_prot_allele,
                                                                              alt_prot_allele, prot_position_end,
                                                                              prot_position_start)
            else:
                protein_change = ""
        else:
            protein_change = TranscriptProviderUtils._render_basic_protein_change(variant_classification, ref_prot_allele,
                                                                              alt_prot_allele, prot_position_end,
                                                                              prot_position_start)
        return protein_change

    @staticmethod
    def render_splice_site_codon_change(dist_from_exon, exon_i):
        if dist_from_exon < 0:
            dist_from_exon = str(dist_from_exon)
        else:
            dist_from_exon = ''.join(['+', str(dist_from_exon)])
        codon_change = 'c.e%d%s' % (exon_i+1, dist_from_exon)
        return codon_change

    @staticmethod
    def render_codon_change(variant_type, variant_classification, codon_position_start, codon_position_end, ref_codon_seq, alt_codon_seq, dist_from_exon, exon_i):
        """
        :param tx:
        :param variant_type:
        :param variant_classification: (str)
        :param codon_position_start:
        :param codon_position_end:
        :param ref_codon_seq:
        :param alt_codon_seq:
        :return:
        """
        updated_ref_seq = ref_codon_seq
        updated_alt_seq = alt_codon_seq
        codon_change = ''
        if variant_classification.startswith('Frame_Shift'):
            codon_change = 'c.(%d-%d)%sfs' % (codon_position_start, codon_position_end, updated_ref_seq)

        elif variant_classification == VariantClassification.SPLICE_SITE:
            codon_change = TranscriptProviderUtils.render_splice_site_codon_change(dist_from_exon, exon_i)

        elif variant_type.endswith('NP'):
            codon_change = 'c.(%d-%d)%s>%s' % (codon_position_start, codon_position_end, updated_ref_seq, updated_alt_seq)
        elif variant_type == 'DEL':
            if alt_codon_seq == '': #full codon deleted
                codon_change = 'c.(%d-%d)%sdel' % (codon_position_start, codon_position_end, updated_ref_seq)
            else:
                codon_change = 'c.(%d-%d)%s>%s' % (codon_position_start, codon_position_end, updated_ref_seq, updated_alt_seq)
        elif variant_type == 'INS':
            if ref_codon_seq == '': #insertion between codons
                codon_change = 'c.(%d-%d)ins%s' % (codon_position_start, codon_position_end, alt_codon_seq)
            else:
                codon_change = 'c.(%d-%d)%s>%s' % (codon_position_start, codon_position_end, updated_ref_seq, updated_alt_seq)
        return codon_change

    # TODO: These transforms should be in a separate class.
    @staticmethod
    def _transform_to_feature_space(exons, s, strand):
        """
        Assumes exons are in order given strand, though their start > end
        :param exons:
        :param s:
        :param strand:
        :return:
        """
        d = 0
        for exon in exons:

            if strand =="+":
                if s > exon[0] and s <= exon[1]:
                    d += (s - exon[0])
                    break
            else:
                if s >= exon[0] and s < exon[1]:
                    d += (exon[1] - s)
                    break

            if (s > exon[1] and strand == "+") or (s < exon[0] and strand == "-"):
                d += (exon[1] - exon[0])
            else:
                break
        return d

    @staticmethod
    def _convert_genomic_space_to_feature_space(s, e, features, strand):
        d_start = TranscriptProviderUtils._transform_to_feature_space(features, s, strand)
        d_end = TranscriptProviderUtils._transform_to_feature_space(features, e, strand)
        if strand == "-":
            tmp = d_start
            d_start = d_end
            d_end = tmp
        return d_start, d_end

    @staticmethod
    def convert_genomic_space_to_exon_space(start, end, tx):
        s = int(start)
        e = int(end)
        strand = tx.get_strand()
        exons = tx.get_exons()
        return TranscriptProviderUtils._convert_genomic_space_to_feature_space(s, e, exons, strand)

    @staticmethod
    def convert_genomic_space_to_cds_space(start, end, tx):
        s = int(start)
        e = int(end)
        strand = tx.get_strand()
        cds = tx.get_cds()
        return TranscriptProviderUtils._convert_genomic_space_to_feature_space(s, e, cds, strand)

    @staticmethod
    def convert_genomic_space_to_transcript_space(start, end, tx):
        """ start <= end, regardless of strand for this method
        This includes all exons, UTR, but not padding.

        :param start: str position in genome space  start <= end
        :param end: str position in genome space    start <= end
        :return (list of str and length 2) ["-1", "-1"] if position cannot be mapped.  Note that here start cannot be greater or less than
            end (start, end).  The strand of transcript is taken into account.
        """
        s = int(start)
        e = int(end)

        tx_start = int(tx.determine_transcript_start())
        tx_end = int(tx.determine_transcript_stop())

        result = ["-1", "-1"]

        if tx.get_strand() == "-":
            if e <= tx_start:
                result[0] = tx_start - e
            if s >= tx_end:
                result[1] = tx_start - s
        else:
            if e <= tx_end:
                result[1] = e - tx_start
            if s >= tx_start:
                result[0] = s - tx_start
        return result[0], result[1]

    @staticmethod
    def get_protein_positions(transcript_position_start, transcript_position_end, cds_start):
        """Parameters are all in transcript space """
        protein_position_start = int(math.ceil((float(transcript_position_start) - float(cds_start) + 1) /3))
        if transcript_position_end == transcript_position_start:
            protein_position_end = protein_position_start
        else:
            protein_position_end = int(math.ceil((float(transcript_position_end) - float(cds_start) + 1) /3))
        return protein_position_start, protein_position_end

    @staticmethod
    def get_cds_codon_positions(protein_position_start, protein_position_end, cds_start):
        """cds_Start is in transcript space """
        cds_codon_start = (protein_position_start * 3 - 2) + cds_start - 1
        cds_codon_end = (protein_position_end * 3) + cds_start - 1
        return cds_codon_start, cds_codon_end

    @staticmethod
    def mutate_reference_sequence(seq, seq_start_pos, mutation_start_pos, mutation_end_pos, observed_allele, variant_type):
        """
        :param seq: (str) the sequence (stranded).
        :param seq_start_pos: (int) Start position of the sequence.  Would normally be zero, but does not have to be.
        :param mutation_start_pos: (int)  This is not relative to the seq_start_pos.  This is on the same scale.
            In other words, this is in exon space.
        :param mutation_end_pos: (int) This is not relative to the seq_start_pos.  This is on the same scale.
            In other words, this is in exon space.
        :param observed_allele: (str) alt allele, stranded.
        :param variant_type:
        :return: Updated sequence using the observed allele
        """
        mutated_seq = list(seq)
        if variant_type == 'DEL':
            del(mutated_seq[mutation_start_pos - seq_start_pos: \
                mutation_end_pos - seq_start_pos + 1])
        elif variant_type == 'INS':
            mutated_seq.insert(mutation_end_pos - seq_start_pos - 1, observed_allele)
        else: #SNP or ONP
            mutated_seq[mutation_start_pos - seq_start_pos: \
                mutation_end_pos - seq_start_pos + 1] = str(observed_allele)
        mutated_seq = ''.join(mutated_seq)
        return mutated_seq

    @staticmethod
    def determine_cds_in_exon_space(tx):
        cds_start_genomic_space, cds_stop_genomic_space = tx.determine_cds_footprint()
        cds_start, cds_stop = TranscriptProviderUtils.convert_genomic_space_to_exon_space(cds_start_genomic_space, cds_stop_genomic_space, tx)
        return cds_start, cds_stop

    @staticmethod
    def render_splice_site_transcript_change(tx, dist_from_exon, exon_i, is_intron):
        if tx.get_strand() == "-":
            dist_from_exon *= -1
        if dist_from_exon < 1:
            cDNA_position = tx.get_exons()[exon_i][0]
        else:
            cDNA_position = tx.get_exons()[exon_i][1]

        cDNA_position_in_exon, dummy = TranscriptProviderUtils.convert_genomic_space_to_cds_space(cDNA_position, cDNA_position, tx)
        if is_intron:
            cDNA_position_in_exon += 1

        return "c.%d_splice" % (cDNA_position_in_exon)
