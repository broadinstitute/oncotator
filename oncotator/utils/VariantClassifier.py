import Bio
from Bio import Seq
import itertools
import math
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.utils.gaf_annotation import chop


class VariantClassifier(object):

    def __init__(self):
        pass

    def _adjust_protein_position_and_alleles(self, protein_seq, protein_position_start, protein_position_end, reference_aa, observed_aa):
        ## Adjust positions and alleles for indels and ONPs when leading or trailing amino acids
        ## are the same in both the refence and observed, thus not changed in the final mutant protein.
        ## e.g. L > LK at position 813 should be - > L at position 814
        if reference_aa == observed_aa: #don't adjust on silent DNPs
            return reference_aa, observed_aa, protein_position_start, protein_position_end
        if reference_aa == '': reference_aa = '-'
        if observed_aa == '': observed_aa = '-'


        adjust_alleles_and_positions, is_reverse = False, False
        if reference_aa[0] == observed_aa[0] and reference_aa[-1] == observed_aa[-1]:
            left_adjust_count, right_adjust_count = 0, 0
            for i, alleles in enumerate(itertools.izip(reference_aa, observed_aa)):
                if alleles[0] == alleles[1]:
                    left_adjust_count += 1
                else:
                    break
            for i, alleles in enumerate(itertools.izip(reference_aa[::-1], observed_aa[::-1])):
                if alleles[0] == alleles[1]:
                    right_adjust_count += 1
                else:
                    break

            if left_adjust_count > right_adjust_count:
                adjust_alleles_and_positions = True
                adjusted_ref_aa, adjusted_obs_aa = list(reference_aa), list(observed_aa)
            elif right_adjust_count > left_adjust_count:
                adjust_alleles_and_positions = True
                is_reverse = True
                adjusted_ref_aa, adjusted_obs_aa = list(reference_aa[::-1]), list(observed_aa[::-1])

        elif reference_aa[0] == observed_aa[0]:
            adjust_alleles_and_positions = True
            adjusted_ref_aa, adjusted_obs_aa = list(reference_aa), list(observed_aa)
        elif reference_aa[-1] == observed_aa[-1]:
            adjust_alleles_and_positions = True
            is_reverse = True
            adjusted_ref_aa, adjusted_obs_aa = list(reference_aa[::-1]), list(observed_aa[::-1])

        if adjust_alleles_and_positions:
            allele_idxs_to_delete = list()
            for i, alleles in enumerate(itertools.izip(adjusted_ref_aa, adjusted_obs_aa)):
                if alleles[0] == alleles[1]:
                    allele_idxs_to_delete.append(i)
                else:
                    break

            for idx in allele_idxs_to_delete[::-1]:
                del(adjusted_ref_aa[idx])
                del(adjusted_obs_aa[idx])

            if is_reverse:
                adjusted_ref_aa, adjusted_obs_aa = adjusted_ref_aa[::-1], adjusted_obs_aa[::-1]
                adjusted_start = protein_position_start
                adjusted_end = protein_position_end - len(allele_idxs_to_delete)
            else:
                adjusted_start = protein_position_start + len(allele_idxs_to_delete)
                adjusted_end = protein_position_end

            adjusted_ref_aa, adjusted_obs_aa = ''.join(adjusted_ref_aa), ''.join(adjusted_obs_aa)
            if adjusted_ref_aa == '': adjusted_ref_aa = '-'
            if adjusted_obs_aa == '': adjusted_obs_aa = '-'
            if adjusted_ref_aa == '-': #insertion between codons:
                adjusted_start, adjusted_end = min(adjusted_start, adjusted_end), max(adjusted_start, adjusted_end)
                if adjusted_end-adjusted_start != 1: raise Exception

            reference_aa, observed_aa = adjusted_ref_aa, adjusted_obs_aa
            protein_position_start, protein_position_end = adjusted_start, adjusted_end

        ## Shift protein postion upstream if adjacent amino acids are the same
        look_upstream_and_adjust_position = False
        if reference_aa == '-': #insertion
            slice = len(observed_aa)
            pos = protein_position_end
            q_aa = observed_aa
            look_upstream_and_adjust_position = True
        elif observed_aa == '-': #deletion
            slice = len(reference_aa)
            pos = protein_position_end + 1
            q_aa = reference_aa
            look_upstream_and_adjust_position = True
        if look_upstream_and_adjust_position:
            for aa in chop(protein_seq[pos-1:], slice):
                if q_aa == ''.join(aa):
                    protein_position_start += slice
                    protein_position_end += slice
                else:
                    break
        return reference_aa, observed_aa, protein_position_start, protein_position_end


    def infer_variant_classification(self, variant_type, reference_aa, observed_aa, reference_allele, observed_allele, is_frameshift_indel=False, is_splice_site=False):
        """ In a nutshell:
        if indel, then return in frame or frame shift ins/del
        if SNP, then return splice_site (if applicable), otherwise the proper vc.

        :param variant_type:
        :param reference_aa:
        :param observed_aa:
        :param reference_allele:
        :param observed_allele:
        :param is_frameshift_indel:
        :param is_splice_site:
        :return:
        """
        if variant_type == 'INS' or (variant_type == 'ONP' and len(reference_allele) < len(observed_allele)):
            if not is_frameshift_indel:
                vc = 'In_Frame_Ins'
            else:
                vc = "Frame_Shift_Ins"
        elif variant_type == 'DEL' or (variant_type == 'ONP' and len(reference_allele) > len(observed_allele)):
            if not is_frameshift_indel:
                vc = 'In_Frame_Del'
            else:
                vc = "Frame_Shift_Del"
        else:
            if is_splice_site:
                vc = "Splice_Site"
            else:
                if reference_aa == observed_aa:
                    vc = 'Silent'
                elif reference_aa.find('*') > -1 and observed_aa.find('*') == -1:
                    vc = 'Nonstop_Mutation'
                elif reference_aa.find('*') == -1 and observed_aa.find('*') > -1:
                    vc = 'Nonsense_Mutation'
                elif reference_aa != observed_aa:
                    vc = 'Missense_Mutation'
        return vc

    def _determine_cds_in_exon_space(self, tx):
        cds_start_genomic_space, cds_stop_genomic_space = tx.determine_cds_footprint()
        cds_start, cds_stop = TranscriptProviderUtils.convert_genomic_space_to_exon_space(cds_start_genomic_space, cds_stop_genomic_space, tx)
        return cds_start, cds_stop

    def variant_classify_old(self, tx, variant_type, ref_allele, alt_allele, start, end):

        reference_allele, observed_allele = str(ref_allele), str(alt_allele)
        if tx.get_gene_type() == 'miRNA':
            is_mirna = True
        else:
            is_mirna = False

        transcript_position_start, transcript_position_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(start, end, tx)

        if tx.get_strand() == '-':
            reference_allele, observed_allele = Bio.Seq.reverse_complement(reference_allele), Bio.Seq.reverse_complement(observed_allele)

        transcript_seq = tx.get_seq()

        # Fix to correct reference transcript sequence when it differs from genomic reference
        # -1 here to account for zero-base python lists
        if variant_type == 'SNP' and transcript_seq[transcript_position_start-1:transcript_position_end] != reference_allele:
            new_transcript_seq = list(transcript_seq)
            new_transcript_seq[transcript_position_start-1:transcript_position_end] = reference_allele
            transcript_seq = ''.join(new_transcript_seq)
            ref_tx_seq_has_been_changed = True
        else:
            ref_tx_seq_has_been_changed = False

        protein_seq = tx.get_protein_seq()
        cds_start, cds_stop = self._determine_cds_in_exon_space(tx)

        #always use '+' here because strand doesn't matter, since inputs are in transcript space
        cds_overlap_type = TranscriptProviderUtils.test_overlap_with_strand(transcript_position_start, transcript_position_end,
            cds_start, cds_stop, '+')

        if cds_overlap_type == 'a_within_b':
            protein_position_start, protein_position_end = TranscriptProviderUtils.get_protein_positions(transcript_position_start,
                transcript_position_end, cds_start)

            cds_codon_start, cds_codon_end = TranscriptProviderUtils.get_cds_codon_positions(protein_position_start,
                protein_position_end, cds_start)

            reference_codon_seq = transcript_seq[cds_codon_start-1:cds_codon_end]
            reference_aa = protein_seq[protein_position_start-1:protein_position_end]

            is_mut_a_frameshift_indel = self.is_framshift_indel(variant_type, int(start), int(end),  observed_allele)

            if variant_type == 'INS' and protein_position_start != protein_position_end:
                #treat differently if insertion falls between codons
                reference_codon_seq = ''
                mutated_codon_seq = observed_allele
                reference_aa = ''
            else:
                mutated_codon_seq = TranscriptProviderUtils.mutate_reference_sequence(reference_codon_seq,
                    cds_codon_start, transcript_position_start, transcript_position_end, observed_allele, variant_type)

            if ref_tx_seq_has_been_changed:
                reference_aa = Bio.Seq.translate(reference_codon_seq)

            observed_aa = Bio.Seq.translate(mutated_codon_seq)

            is_splice_site = self._determine_if_splice_site(int(start), int(end), tx, dist=2)

            variant_classification = self.infer_variant_classification(variant_type,
                reference_aa, observed_aa, reference_allele, observed_allele, is_frameshift_indel=is_mut_a_frameshift_indel, is_splice_site=is_splice_site)

            # TODO: Bring this code back, so that we can handle alternate values for silent mutations
            # # If silent mutation w/in 2 bp of a splice junction, then change to splice site
            # if variant_classification.lower() == "silent":
            #     self._determine_if_splice_site(int(start), int(end), tx, dist=2)

            if variant_type != 'SNP':
                reference_aa, observed_aa, protein_position_start, protein_position_end = \
                    self._adjust_protein_position_and_alleles(protein_seq, protein_position_start,
                        protein_position_end, reference_aa, observed_aa)
        else:
            variant_classification = self.annotate_mutations_not_fully_within_cds(transcript_seq, observed_allele, cds_overlap_type, transcript_position_start,
                transcript_position_end, variant_type, tx)

        return variant_classification

    def variant_classify(self, tx, variant_type, ref_allele, alt_allele, start, end):
        """
        Important note:  Does not support ONP.
        :param tx:
        :param variant_type:
        :param ref_allele:
        :param alt_allele:
        :param start:
        :param end:
        :return:
        """
        if variant_type == "SNP":
            return self._variant_classify_snp(tx, ref_allele, alt_allele, start, end)
        else:
            return self.variant_classify_old(tx, variant_type, ref_allele, alt_allele, start, end)

    def is_framshift_indel(self, variant_type, start, end,  observed_allele):
        if variant_type in ['DEL','INS']:
            if variant_type == 'DEL':
                dl = (end-start+1) % 3
            elif variant_type == 'INS':
                dl = len(observed_allele) % 3
            if dl != 0:
                return True
        return False

    def _determine_if_splice_site(self, start_genomic_space, end_genomic_space, tx, dist=2):

        """
         If overlap is detected, but the start or end is within dist bp, then this is a splice site.
        :param end:
        :param start:
        :param t:

        """
        exon_i, ldist, rdist = self._determine_closest_distances_from_exon(tx.get_exons(), start_genomic_space, end_genomic_space)

        # If we are on the exon side, we need to correct by one base.
        if start_genomic_space > tx.get_exons()[exon_i][0] and end_genomic_space < tx.get_exons()[exon_i][1]:
            if abs(start_genomic_space - tx.get_exons()[exon_i][0]) > abs(end_genomic_space - tx.get_exons()[exon_i][1]):
                rdist += 1
            if abs(start_genomic_space - tx.get_exons()[exon_i][0]) < abs(end_genomic_space - tx.get_exons()[exon_i][1]):
                ldist += 1

        if abs(ldist) <= dist or abs(rdist) <= dist:
            return True
        return False

    def _determine_closest_distances_from_exon(self, exons, start_genomic_space, end_genomic_space):
        """
        start_genomic_space <= end_genomic_space

        Everything in genomic space.

        :param exons:
        :param start_genomic_space:
        :param end_genomic_space:
        :return: exon_i, ldist, rdist
        ldist is the distance from the closest edge (start/end) to the left side of the exon (lower genomic position)
        rdist is the distance from the closest edge (start/end) to the right side of the exon (lower genomic position)
        """
        exon_i = ldist = rdist = -1
        min_val = 10000000000

        if start_genomic_space > end_genomic_space:
            raise ValueError("start value must be less than (or equal) to end value.  " + str([start_genomic_space,end_genomic_space]))

        for i in range(0, len(exons)):
            exon = exons[i]

            # Distance from the larger genomic position
            left_diff = abs(start_genomic_space - exon[0])
            right_diff = abs(end_genomic_space - exon[1])

            if min_val > left_diff or min_val > right_diff:
                ldist = left_diff
                rdist = right_diff
                exon_i = i
                min_val = min(left_diff, right_diff)

        return exon_i, ldist, rdist

    def annotate_mutations_not_fully_within_cds(self, transcript_seq, observed_allele, cds_overlap_type, transcript_position_start, transcript_position_end, variant_type, t):
        if variant_type in ['DEL','INS']:
            vt = ''.join([variant_type[0], variant_type[1:].lower()])
        else:
            vt = 'variant_type'

        result = ""
        if cds_overlap_type == 'a_overlaps_b_left_border':
            result = 'Start_Codon_' + vt
        elif cds_overlap_type == 'a_overlaps_b_right_border':
            result = 'Stop_Codon_' + vt
        else:
            cds_start_in_exon_space,cds_end_in_exon_space = self._determine_cds_in_exon_space(t)

            if transcript_position_start < cds_start_in_exon_space and transcript_position_start < cds_end_in_exon_space: p = 5
            elif transcript_position_start > cds_start_in_exon_space and transcript_position_start > cds_end_in_exon_space: p = 3
            vc = "%d'UTR" % (p)
            result = vc
            if vc == "5'UTR":
                utr_region_start, utr_region_end = transcript_position_start-2, transcript_position_end+2

                utr_region_seq = transcript_seq[utr_region_start-1:utr_region_end]
                # TranscriptProviderUtils.mutate_reference_sequence(reference_codon_seq,
                #     cds_codon_start, transcript_position_start, transcript_position_end, observed_allele, variant_type)
                mutated_utr_region_seq = TranscriptProviderUtils.mutate_reference_sequence(utr_region_seq, utr_region_start,
                    transcript_position_start, transcript_position_end, observed_allele, variant_type)

                # Check for Denovo
                ATG_position = mutated_utr_region_seq.find('ATG')
                if ATG_position > -1:
                    ATG_position = utr_region_start + ATG_position
                    if (t['cds_start'] - ATG_position) % 3 == 0:
                        frameness = 'InFrame'
                    else:
                        frameness = 'OutOfFrame'
                    result = 'De_novo_Start_' + frameness
        return result

    def _determine_snp_vc_for_cds_overlap(self, start, end, ref_allele, alt_allele, is_splice_site, tx):

        reference_allele_stranded, observed_allele_stranded = ref_allele, alt_allele
        if tx.get_strand() == '-':
            reference_allele_stranded, observed_allele_stranded = Bio.Seq.reverse_complement(ref_allele), Bio.Seq.reverse_complement(alt_allele)

        transcript_position_start, transcript_position_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(
            start, end, tx)
        transcript_seq = tx.get_seq()
        protein_seq = tx.get_protein_seq()
        cds_start, cds_stop = self._determine_cds_in_exon_space(tx)
        protein_position_start, protein_position_end = TranscriptProviderUtils.get_protein_positions(
            transcript_position_start,
            transcript_position_end, cds_start)
        if transcript_seq[transcript_position_start-1:transcript_position_end] != reference_allele_stranded:
            new_transcript_seq = list(transcript_seq)
            new_transcript_seq[transcript_position_start-1:transcript_position_end] = reference_allele_stranded
            transcript_seq = ''.join(new_transcript_seq)
            ref_tx_seq_has_been_changed = True
        else:
            ref_tx_seq_has_been_changed = False
        cds_codon_start, cds_codon_end = TranscriptProviderUtils.get_cds_codon_positions(protein_position_start,
                                                                                         protein_position_end,
                                                                                         cds_start)
        reference_codon_seq = transcript_seq[cds_codon_start-1:cds_codon_end]
        mutated_codon_seq = TranscriptProviderUtils.mutate_reference_sequence(reference_codon_seq,
                                                                              cds_codon_start,
                                                                              transcript_position_start,
                                                                              transcript_position_end, observed_allele_stranded,
                                                                              "SNP")
        # if tx.get_strand() == "-":
        #     # Get the AA for the reversed mutated codon.
        #     observed_aa = Bio.Seq.translate(mutated_codon_seq[::-1])
        # else:
        observed_aa = Bio.Seq.translate(mutated_codon_seq)
        # if ref_tx_seq_has_been_changed:
        #     reference_aa = Bio.Seq.translate(reference_codon_seq)
        # else:
        reference_aa = protein_seq[protein_position_start+1:protein_position_end+2]
        vc_tmp = self.infer_variant_classification("SNP", reference_aa, observed_aa, ref_allele, alt_allele,
                                                   is_frameshift_indel=False, is_splice_site=is_splice_site)
        return vc_tmp

    def _variant_classify_snp(self, tx, ref_allele, alt_allele, start, end, dist=2):
        """Perform classifications for SNPs only.
        Everything handled in genomic space

        *RNA*
        x'UTR
        Splice_Site (Intron)
        Intron
        Splice_Site (Exon)
        {Missense, Silent}
        {Nonsense, Silent}
        {Nonstop, Silent}
        IGR
        x'Flank
        De_novo_Start

        """
        gene_type = tx.get_gene_type()
        if gene_type != "protein_coding":
            return gene_type

        s = int(start)
        e = int(end)
        is_exon_overlap = TranscriptProviderUtils.test_feature_overlap(s, e, tx.get_exons())
        is_splice_site = self._determine_if_splice_site(s, e, tx, dist)
        is_beyond_exons, side, is_flank = self._determine_beyond_exon_info(int(start), int(end), tx)
        if not is_exon_overlap and not is_beyond_exons:
            if is_splice_site:
                # Intron Splice Site
                return "Splice_Site"
            else:
                return "Intron"

        if not is_exon_overlap and is_beyond_exons:
            if is_flank:
                return side + "Flank"
            else:
                return "IGR"

        is_cds_overlap = TranscriptProviderUtils.test_feature_overlap(s, e, tx.get_cds())
        # is_start_codon_overlap = TranscriptProviderUtils.test_overlap(s, e, tx.get_start_codon()[0], tx.get_start_codon()[1])
        # is_stop_codon_overlap = TranscriptProviderUtils.test_overlap(s, e, tx.get_stop_codon()[0], tx.get_stop_codon()[1])

        if is_exon_overlap and not is_cds_overlap:
            # UTR
            vc_tmp = side + "UTR"
            transcript_position_exon_space_start, transcript_position_exon_space_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(start, end, tx)
            vc = self._determine_de_novo(vc_tmp, transcript_position_exon_space_start, transcript_position_exon_space_end, ref_allele, alt_allele, tx, "SNP")
            return vc

        # We have a clean SNP in the CDS.  No start codon or stop codon.
        if is_cds_overlap:
            return self._determine_snp_vc_for_cds_overlap(start, end, ref_allele, alt_allele, is_splice_site, tx)

        raise ValueError("Could not determine variant classification.")


    def _determine_beyond_exon_info(self, start, end, tx, flank_padding_5prime=3000, flank_padding_3prime=0):
        """
        start and end must be completely non-overlapping of the tx exons.  Reminder that tx exons include UTRs.
        Also, this will return False, closest side (3' or 5'), False if Intron.

        :param start: int
        :param end: int
        :param tx: Transcript
        :return: triplet of is_beyond_exons, side (always a str of 3' or 5' or "" if not is_beyond_exons), is_flank
        """
        is_beyond_exons_left = (start < tx.get_start() and end < tx.get_start())
        is_beyond_exons_right = (start > tx.get_end() and end > tx.get_end())
        is_beyond_exons = is_beyond_exons_left or is_beyond_exons_right
        side = self._determine_strand_side(start, end, tx)
        if not is_beyond_exons:
            return False, side, False
        if is_beyond_exons_left:
            d = min(abs(start - tx.get_start()), abs(end - tx.get_start()))
        else:
            d = min(abs(start - tx.get_end()), abs(end - tx.get_end()))

        if side == "5'":
            is_flank = (d < flank_padding_5prime)
        else:
            is_flank = (d < flank_padding_3prime)
        return is_beyond_exons, side, is_flank

    def _determine_strand_side(self, start, end, tx):
        """ Determine the closest side (3' or 5') given start and end (in genomic space) and a transcript.
        This method takes into account transcript strand

        Checks the transcript beginning and end for distance to start and end.

        :param start: int
        :param end: int
        :param tx: Transcript
        :return: str ("3'" or "5'")
        """
        strand = tx.get_strand()
        is_beyond_exons_left = (start < tx.get_start() and end < tx.get_start())
        is_beyond_exons_right = (start > tx.get_end() and end > tx.get_end())
        if is_beyond_exons_left and strand == "+":
            return "5'"
        if is_beyond_exons_left and strand == "-":
            return "3'"
        if is_beyond_exons_right and strand == "+":
            return "3'"
        if is_beyond_exons_right and strand == "-":
            return "5'"

        d_left = min(abs(start - tx.get_start()), abs(end - tx.get_start()))
        d_right = min(abs(start - tx.get_end()), abs(end - tx.get_end()))

        if d_left <= d_right:
            if strand == "+":
                return "5'"
            else:
                return "3'"
        else:
            if strand == "+":
                return "3'"
            else:
                return "5'"

    def _determine_de_novo(self, vc, transcript_position_start, transcript_position_end, ref, alt, tx, variant_type):
        """Returns input vc if not de Novo.  Otherwise, returns updated variant classification.

         Will always return original vc if the vc is not None."""
        result = vc

        if vc == "5'UTR" and ref != alt:
            tx_seq = tx.get_seq()
            utr_region_start, utr_region_end = transcript_position_start-2, transcript_position_end+2

            utr_region_seq = tx_seq[utr_region_start-1:utr_region_end]

            mutated_utr_region_seq = TranscriptProviderUtils.mutate_reference_sequence(utr_region_seq, utr_region_start,
                transcript_position_start, transcript_position_end, alt, variant_type)

            # Check for Denovo
            ATG_position = mutated_utr_region_seq.find('ATG')
            if ATG_position > -1:
                ATG_position = utr_region_start + ATG_position
                if (t['cds_start'] - ATG_position) % 3 == 0:
                    frameness = 'InFrame'
                else:
                    frameness = 'OutOfFrame'
                result = 'De_novo_Start_' + frameness
        return result