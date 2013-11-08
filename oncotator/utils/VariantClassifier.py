import Bio
from Bio import Seq
import itertools
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


    def infer_variant_classification(self, variant_type, reference_aa, observed_aa, reference_allele, observed_allele):
        if variant_type == 'INS' or (variant_type == 'ONP' and len(reference_allele) < len(observed_allele)):
            vc = 'In_Frame_Ins'
        elif variant_type == 'DEL' or (variant_type == 'ONP' and len(reference_allele) > len(observed_allele)):
            vc = 'In_Frame_Del'
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
        cds_start_genomic_space, cds_stop_genomic_space = self._determine_cds_footprint(tx)
        cds_start, cds_stop = TranscriptProviderUtils.convert_genomic_space_to_exon_space(cds_start_genomic_space, cds_stop_genomic_space, tx)
        return cds_start, cds_stop

    def _determine_cds_footprint(self, tx):
        """ Returns the cds in genomic space.  Note that strand is ignored, so the first return value is always lower.
        :param tx:
        :return: cds_start, cds_stop in genomic coordinates.
        """
        s = cds_start = tx.determine_cds_start()
        e = cds_stop = tx.determine_cds_stop()
        if cds_stop < cds_start:
            s = cds_stop
            e = cds_start
        return s, e

    def _determine_protein_seq(self, tx):
        cds_start, cds_stop = self._determine_cds_footprint(tx)
        protein_seq = self.get_protein_sequence(tx, cds_start, cds_stop)
        protein_seq = ''.join([protein_seq, '*'])
        return protein_seq

    def variant_classify(self, tx, variant_type, ref_allele, alt_allele, start, end):

        reference_allele, observed_allele = str(ref_allele), str(alt_allele)
        if tx.get_gene_type() == 'miRNA':
            is_mirna = True
        else:
            is_mirna = False


        transcript_position_start, transcript_position_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(start, end, tx)

        # TODO: Fix exon_affected
        exon_affected = "-1"

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

        protein_seq = self._determine_protein_seq(tx)
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
            # TODO: Handle a frameshift here, since the rest of the code will not handle it.

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

            variant_classification = self.infer_variant_classification(variant_type,
                reference_aa, observed_aa, reference_allele, observed_allele)

            # If silent mutation w/in 2 bp of a splice junction, then change to splice site
            if variant_classification.lower() == "silent":
                determineIfSpliceSiteThrown(end, exon_affected,gaf,start,t, dist=2)

            if variant_type != 'SNP':
                reference_aa, observed_aa, protein_position_start, protein_position_end = \
                    self._adjust_protein_position_and_alleles(protein_seq, protein_position_start,
                        protein_position_end, reference_aa, observed_aa)
        else:
            annotate_mutations_not_fully_within_cds(transcript_seq, observed_allele, cds_overlap_type, transcript_position_start,
                transcript_position_end, m, t, exon_affected)

        return variant_classification

    def get_protein_sequence(self, tx, cds_start_genomic_space, cds_stop_genomic_space):
        tx_seq = tx.get_seq()
        cds_start_exon_space, cds_stop_exon_space = TranscriptProviderUtils.convert_genomic_space_to_exon_space(cds_start_genomic_space, cds_stop_genomic_space, tx)

        prot_seq = Seq.translate(tx_seq[int(cds_start_exon_space):int(cds_stop_exon_space)])
        if prot_seq[-1] == '*':
            prot_seq = prot_seq[:-1]

        return prot_seq

    def is_framshift_indel(self, variant_type, start, end,  observed_allele):
        if variant_type in ['DEL','INS']:
            if variant_type == 'DEL':
                dl = (end-start+1) % 3
            elif variant_type == 'INS':
                dl = len(observed_allele) % 3
            if dl != 0:
                return True
        return False