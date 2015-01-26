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

import collections
import math
from oncotator.utils.VariantClassification import VariantClassification


class TranscriptProviderUtils(object):

    @staticmethod
    def is_valid_xNP(variant_type, ref, alt):
        """Call this method when this may be an invalid variant type, but it is xNP.
        :param variant_type: One of the VariantClassificaton.VT_SNP, VariantClassificaton.VT_DNP,
            VariantClassificaton.VT_TNP, or VariantClassificaton.VT_ONP
        :param ref: (str) ref allele
        :param alt: (str) alt allele

        Adapted from http://code.activestate.com/recipes/499304-hamming-distance/

        Returns a string with a message if there is a difference.  Empty string means no issue seen.

        """
        diffs = 0
        for ch1, ch2 in zip(ref, alt):
                if ch1 != ch2:
                        diffs += 1

        if ref != "-" and alt != "-" and len(ref) != len(alt):
            return "ref and alt are of different lengths.  This is an indel, but variant type is: %s" % (variant_type)

        message = ""

        if diffs == 0:
            message = "No difference between ref and alt."
        if diffs == 1 and variant_type != VariantClassification.VT_SNP:
            message = "%s>%s is a SNP, but is found to be: %s." %(ref, alt, variant_type)
        if diffs == 2 and variant_type != VariantClassification.VT_DNP:
            message = "%s>%s is a DNP, but is found to be: %s." %(ref, alt, variant_type)
        if diffs == 3 and variant_type != VariantClassification.VT_TNP:
            message = "%s>%s is a TNP, but is found to be: %s." %(ref, alt, variant_type)
        if diffs > 3 and variant_type != VariantClassification.VT_ONP:
            message = "%s>%s is a ONP, but is found to be: %s." %(ref, alt, variant_type)

        return message

    @staticmethod
    def is_xnp(variant_type):
        return variant_type.endswith(VariantClassification.VT_xNP)

    @staticmethod
    def infer_variant_type(reference_allele, observed_allele):
        """ To go completely in annotate method.  Returns variant type string."""
        if (reference_allele == '-' or reference_allele == '') and (observed_allele == '-' or observed_allele == ''):
            raise Exception('Variant Type cannot be inferred from reference and observed allele because both are empty: (%s, %s)' % (reference_allele, observed_allele))

        if reference_allele == '-' or reference_allele == "": #is insertion
            return VariantClassification.VT_INS
        elif observed_allele == '-' or observed_allele == "": #is deletion
            return VariantClassification.VT_DEL
        else:
            if len(observed_allele) > len(reference_allele):
                return VariantClassification.VT_INS
            elif len(observed_allele) < len(reference_allele):
                return VariantClassification.VT_DEL
            elif len(reference_allele) == 1:
                return VariantClassification.VT_SNP
            elif len(reference_allele) == 2:
                return VariantClassification.VT_DNP
            elif len(reference_allele) == 3:
                return VariantClassification.VT_TNP
            elif len(reference_allele) > 3:
                return VariantClassification.VT_ONP

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



    NO_APPRIS_VALUE = 100
    #Used for ordering transcripts based on appris values
    #See http://www.gencodegenes.org/gencode_tags.html
    APPRIS_TAGS =[ 'appris_principal',
            'appris_candidate_highest_score',
            'appris_candidate_longest_ccds',
            'appris_candidate_ccds',
            'appris_candidate_longest_seq',
            'appris_candidate_longest',
            'appris_candidate',
            None]
    APPRIS_RANKING_DICT = collections.OrderedDict()
    for (i, tag) in enumerate(APPRIS_TAGS):
        APPRIS_RANKING_DICT[tag] = i



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
    def render_transcript_change(variant_type, variant_classification, exon_position_start, exon_position_end, ref_allele_stranded, alt_allele_stranded, secondary_vc):
        """



        :param variant_type:
        :param variant_classification:
        :param exon_position_start: Coordinates in transcript/exon space
        :param exon_position_end: Coordinates in transcript/exon space
        :param ref_allele_stranded: ref_allele with strand already accounted
        :param alt_allele_stranded: alt_allele with strand already accounted
        :param secondary_vc:
        """
        transcript_change = ""
        if variant_classification.startswith(VariantClassification.SPLICE_SITE) and secondary_vc == VariantClassification.INTRON:
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
    def render_intronic_splice_site_codon_change(dist_from_exon, exon_i):
        if dist_from_exon < 0:
            dist_from_exon = str(dist_from_exon)
        else:
            dist_from_exon = ''.join(['+', str(dist_from_exon)])
        codon_change = 'c.e%d%s' % (exon_i+1, dist_from_exon)
        return codon_change

    @staticmethod
    def render_codon_change(variant_type, variant_classification, codon_position_start, codon_position_end, ref_codon_seq, alt_codon_seq, dist_from_exon, exon_i, secondary_vc):
        """
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

        elif variant_classification == VariantClassification.SPLICE_SITE and secondary_vc == VariantClassification.INTRON:
            codon_change = TranscriptProviderUtils.render_intronic_splice_site_codon_change(dist_from_exon, exon_i)

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
            del(mutated_seq[mutation_start_pos - seq_start_pos:
                mutation_end_pos - seq_start_pos + 1])
        elif variant_type == 'INS':
            mutated_seq.insert(mutation_end_pos - seq_start_pos - 1, observed_allele)
        else: #SNP or ONP
            mutated_seq[mutation_start_pos - seq_start_pos:
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

    @staticmethod
    def determine_closest_distance_from_exon(start_genomic, end_genomic, exon_i,  t):
        """
        Return the distance of the given position to the exon in the transcript.

        :param int start_genomic: start position in genomic space
        :param int end_genomic: end position in genomic space (inclusive)
        :param int exon_i: exon index (0-based) of the transcript t
        :param Transcript t: transcript
        :return tuple: left_diff, right_diff -- distance from specified position from the left- (and right-) most position of the exon in genomic space.
        """
        left_start_diff = t.get_exons()[exon_i][0] - start_genomic
        left_end_diff = t.get_exons()[exon_i][0] - end_genomic
        right_start_diff = t.get_exons()[exon_i][1] - start_genomic
        right_end_diff = t.get_exons()[exon_i][1] - end_genomic
        left_diff = min(left_start_diff, left_end_diff)
        right_diff = max(right_start_diff, right_end_diff)
        return left_diff, right_diff

    @staticmethod
    def determine_closest_exon(t, start, end):
        """
        Return the closest exon index (0-based) in the given transcript for the given start and end position.

        :param Transcript t: transcript
        :param int start: genomic coordinate for start position
        :param int end: genomic coordinate for end position
        :return int: exon index.  None if t is None
        """
        if t is None:
            return None
        tmp_distances = []
        for i,exon in enumerate(t.get_exons()):
            left_diff, right_diff = TranscriptProviderUtils.determine_closest_distance_from_exon(start, end, i,  t)
            tmp_distances.append((abs(left_diff), abs(right_diff)))

        min_dist = float("Inf")
        min_index = -1
        for i,ds in enumerate(tmp_distances):
            if ds[0] < min_dist:
                min_dist = ds[0]
                min_index = i
            if ds[1] < min_dist:
                min_dist = ds[1]
                min_index = i
        return min_index


    @staticmethod
    def determine_if_exon_overlap(s, e, tx, variant_type):
        return TranscriptProviderUtils.determine_exon_index(s, e, tx, variant_type) != -1

    @staticmethod
    def determine_exon_index(s, e, tx, variant_type):
        if variant_type != VariantClassification.VT_INS:
            return TranscriptProviderUtils.test_feature_overlap(s, e, tx.get_exons())
        else:
            return TranscriptProviderUtils.test_feature_overlap(s, s, tx.get_exons())