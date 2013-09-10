class TranscriptProviderUtils(object):

    @staticmethod
    def infer_variant_type(reference_allele, observed_allele):
        """ To go completely in annotate method.  Returns variant type string."""
        if reference_allele == '-': #is insertion
            return 'INS'
        elif observed_allele == '-': #is deletion
            return 'DEL'
        else:
            if len(observed_allele) != len(reference_allele):
                return 'ONP'
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
