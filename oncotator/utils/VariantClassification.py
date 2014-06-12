
class VariantClassification(object):
    """ Stores the names and values of a variant classification.
    """
    INTRON = "Intron"
    FIVE_PRIME_UTR = "5'UTR"
    THREE_PRIME_UTR = "3'UTR"
    IGR = "IGR"
    FIVE_PRIME_PRIME_FLANK = "5'Flank"
    THREE_PRIME_PRIME_FLANK = "3'Flank"
    MISSENSE = "Missense_Mutation"
    NONSENSE = "Nonsense_Mutation"
    NONSTOP = "Nonstop_Mutation"
    SILENT = "Silent"
    SPLICE_SITE = "Splice_Site"
    IN_FRAME_DEL = "In_Frame_Del"
    IN_FRAME_INS = "In_Frame_Ins"
    FRAME_SHIFT_INS = "Frame_Shift_Ins"
    FRAME_SHIFT_DEL = "Frame_Shift_Del"
    START_CODON_SNP = "Start_Codon_SNP"
    START_CODON_INS = "Start_Codon_Ins"
    START_CODON_DEL = "Start_Codon_Del"
    STOP_CODON_INS = "Stop_Codon_Ins"
    STOP_CODON_DEL = "Stop_Codon_Del"
    # Note: A STOP_CODON_SNP is a nonstop mutation (or Silent)
    DE_NOVO_START_IN_FRAME = "De_novo_Start_InFrame"
    DE_NOVO_START_OUT_FRAME = "De_novo_Start_OutOfFrame"
    RNA = "RNA"
    LINCRNA = "lincRNA"

    VT_INS = "INS"
    VT_DEL = "DEL"
    VT_SNP = "SNP"

    def __init__(self, vc_primary, vt, transcript_id="", vc_secondary="", alt_codon="", alt_codon_start_in_exon="",
        alt_codon_end_in_exon="", ref_codon="", ref_codon_start_in_exon="", ref_codon_end_in_exon="", ref_aa="",
        ref_protein_start="", ref_protein_end="", alt_aa="", alt_protein_start="", alt_protein_end="",
        cds_start_in_exon_space="", ref_allele_stranded="", alt_allele_stranded="", exon_i=-1):
    
        self._vc_primary = vc_primary
        self._vc_secondary = vc_secondary
        self._transcript_id = transcript_id
        self._alt_codon = alt_codon
        self._ref_codon = ref_codon
        self._ref_aa = ref_aa
        self._ref_protein_start = ref_protein_start
        self._ref_protein_end = ref_protein_end
        self._alt_aa = alt_aa
        self._alt_protein_start = alt_protein_start
        self._alt_protein_end = alt_protein_end
        self._vt = vt
        self._alt_codon_start_in_exon = alt_codon_start_in_exon
        self._alt_codon_end_in_exon = alt_codon_end_in_exon
        self._ref_codon_start_in_exon = ref_codon_start_in_exon
        self._ref_codon_end_in_exon = ref_codon_end_in_exon
        self._cds_start_in_exon_space = cds_start_in_exon_space
        self._ref_allele_stranded = ref_allele_stranded
        self._alt_allele_stranded = alt_allele_stranded
        self._exon_i = exon_i

    def get_vc(self):
        """Returns primary variant classification as string """
        return self._vc_primary

    def get_secondary_vc(self):
        return self._vc_secondary

    def get_transcript_id(self):
        return self._transcript_id

    def get_alt_codon(self):
        return self._alt_codon

    def get_vt(self):
        return self._vt

    def get_ref_codon(self):
        return self._ref_codon

    def get_ref_aa(self):
        return self._ref_aa

    def get_ref_protein_start(self):
        return self._ref_protein_start

    def get_ref_protein_end(self):
        return self._ref_protein_end

    def get_alt_aa(self):
        return self._alt_aa

    def get_alt_protein_start(self):
        return self._alt_protein_start

    def get_alt_protein_end(self):
        return self._alt_protein_end

    def get_alt_codon_start_in_exon(self):
        return self._alt_codon_start_in_exon

    def get_alt_codon_end_in_exon(self):
        return self._alt_codon_end_in_exon

    def get_ref_codon_start_in_exon(self):
        return self._ref_codon_start_in_exon

    def get_ref_codon_end_in_exon(self):
        return self._ref_codon_end_in_exon

    def get_cds_start_in_exon_space(self):
        return self._cds_start_in_exon_space

    def get_ref_allele_stranded(self):
        return self._ref_allele_stranded

    def get_alt_allele_stranded(self):
        return self._alt_allele_stranded

    def get_exon_i(self):
        return self._exon_i