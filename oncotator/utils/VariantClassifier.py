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
import logging

import Bio
import itertools
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.utils.InvalidVariantException import InvalidVariantException
from oncotator.utils.VariantClassification import VariantClassification
from oncotator.utils.gaf_annotation import chop
from oncotator.utils.MutUtils import MutUtils


class VariantClassifier(object):

    # TODO: Make these configurable
    FLANK_PADDING_5PRIME = 3000
    FLANK_PADDING_3PRIME = 0

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
                if adjusted_end-adjusted_start != 1:
                    raise InvalidVariantException("Insertion of protein between codons cannot be rendered properly.")

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


    def infer_variant_classification(self, variant_type, reference_aa, observed_aa, reference_allele, observed_allele, is_frameshift_indel=False, is_splice_site=False, is_start_codon=False):
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
        :return: primary classification, secondary classification
        """
        #TODO: Replace the rest with constants in the VariantClassification class
        #TODO: Cleanup since the flag to start codon is a bit of a hack.
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
            if reference_aa == observed_aa:
                vc = 'Silent'
            elif reference_aa.find('*') > -1 and observed_aa.find('*') == -1:
                vc = 'Nonstop_Mutation'
            elif reference_aa.find('*') == -1 and observed_aa.find('*') > -1:
                vc = 'Nonsense_Mutation'
            elif reference_aa != observed_aa:
                vc = VariantClassification.MISSENSE
        if is_splice_site:
            return VariantClassification.SPLICE_SITE, vc
        if is_start_codon:
            return VariantClassification.START_CODON_SNP, vc
        return vc, ""


    def is_frameshift_indel(self, variant_type, start, end,  observed_allele):
        if variant_type in ['DEL','INS']:
            if variant_type == 'DEL':
                dl = (end-start+1) % 3
            elif variant_type == 'INS':
                dl = len(observed_allele) % 3
            if dl != 0:
                return True
        return False

    def _determine_if_splice_site_overlap(self, start_genomic_space, end_genomic_space, tx, variant_type, dist=2):

        """

        Overlap of start and stop codon (i.e. start of first exon and end of last exon -- stranded) will not be a
            Splice_Site.  This method will return is_splice_site_overlap of False

         If overlap is detected, but the start or end is within dist bp, then this is a splice site.
         start <= end
        INS events only call splice site when they start in the splice site

        :param start_genomic_space: int in genomic space
        :param end_genomic_space: int in genomic space
        :param tx: Transcript
        :param variant_type:
        :param dist:
        :return is_splice_site_overlap, exon_i, is_right_overlap (Higher genomic position --> True)

        """
        exons = tx.get_exons()
        strand = tx.get_strand()

        # If this is an insertion, we only want to count a splice site if it starts in the splice site regions
        if variant_type == VariantClassification.VT_INS:
            end_genomic_space = start_genomic_space

        for i,exon in enumerate(exons):
            is_internal_exon = (i > 0) and (i < (len(exons)-1))
            is_check_left = is_internal_exon or (strand == "-" and i == 0) or (strand == "+" and i == (len(exons)-1))
            is_check_right = is_internal_exon or (strand == "+" and i == 0) or (strand == "-" and i == (len(exons)-1))
            if is_check_left:
                splice_site_left = (exon[0]-dist+1, exon[0]+(dist-1)+1)
                overlap_type_left = TranscriptProviderUtils.test_overlap(start_genomic_space, end_genomic_space, splice_site_left[0], splice_site_left[1])
                if overlap_type_left:
                    return True, i, False
            if is_check_right:
                splice_site_right = (exon[1]-(dist-1), exon[1] + dist)
                overlap_type_right = TranscriptProviderUtils.test_overlap(start_genomic_space, end_genomic_space, splice_site_right[0], splice_site_right[1])
                if overlap_type_right:
                    return True, i, True

        return False, -1, None, False

    def _get_stranded_alleles(self, ref_allele, alt_allele, tx):
        reference_allele_stranded, observed_allele_stranded = ref_allele, alt_allele
        if tx.get_strand() == '-':
            reference_allele_stranded, observed_allele_stranded = Bio.Seq.reverse_complement(
                ref_allele), Bio.Seq.reverse_complement(alt_allele)
        return observed_allele_stranded, reference_allele_stranded

    def _determine_vc_for_cds_overlap(self, start, end, ref_allele, alt_allele, is_frameshift_indel, is_splice_site, tx, variant_type, is_start_codon):
        """
        Note: This method can also handle start and stop codons.

        :param start:
        :param end:
        :param ref_allele:
        :param alt_allele:
        :param is_frameshift_indel:
        :param is_splice_site:
        :param tx:
        :param variant_type:
        :return:
        """
        observed_allele_stranded, reference_allele_stranded = self._get_stranded_alleles(ref_allele, alt_allele, tx)
        transcript_position_start, transcript_position_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(
            start, end, tx)

        if tx.get_strand() == "+" and not variant_type == VariantClassification.VT_INS:
            transcript_position_start -= 1
            transcript_position_end -= 1

        transcript_seq = tx.get_seq()
        protein_seq = tx.get_protein_seq()
        cds_start, cds_stop = TranscriptProviderUtils.determine_cds_in_exon_space(tx)
        protein_position_start, protein_position_end = TranscriptProviderUtils.get_protein_positions(
            transcript_position_start,
            transcript_position_end, cds_start)
        new_ref_transcript_seq = transcript_seq
        if (transcript_seq[transcript_position_start:transcript_position_end+1] != reference_allele_stranded) and variant_type != VariantClassification.VT_INS:
            new_ref_transcript_seq = list(transcript_seq)
            new_ref_transcript_seq[transcript_position_start:transcript_position_end+1] = reference_allele_stranded
            new_ref_transcript_seq = ''.join(new_ref_transcript_seq)
            ref_tx_seq_has_been_changed = True
        else:
            ref_tx_seq_has_been_changed = False
        cds_codon_start, cds_codon_end = TranscriptProviderUtils.get_cds_codon_positions(protein_position_start, protein_position_end, cds_start)

        if variant_type == "DEL":
            reference_codon_seq = new_ref_transcript_seq[cds_codon_start:cds_codon_end+1].lower()
        else:
            reference_codon_seq = TranscriptProviderUtils.mutate_reference_sequence(new_ref_transcript_seq[cds_codon_start:cds_codon_end+1].lower(), cds_codon_start, transcript_position_start, transcript_position_end, reference_allele_stranded, variant_type)

        if variant_type == "INS" and tx.get_strand() == "-":
            mutated_codon_seq = TranscriptProviderUtils.mutate_reference_sequence(reference_codon_seq.lower(), cds_codon_start - 1, transcript_position_start, transcript_position_end, observed_allele_stranded, variant_type)
        else:
            mutated_codon_seq = TranscriptProviderUtils.mutate_reference_sequence(reference_codon_seq.lower(), cds_codon_start, transcript_position_start, transcript_position_end, observed_allele_stranded, variant_type)


        observed_aa = MutUtils.translate_sequence(mutated_codon_seq)
        if ref_tx_seq_has_been_changed:
            reference_aa = MutUtils.translate_sequence(reference_codon_seq)
        else:
            reference_aa = protein_seq[protein_position_start-1:protein_position_end]

        if variant_type != VariantClassification.VT_SNP:

            try:
                reference_aa, observed_aa, protein_position_start, protein_position_end = \
                    self._adjust_protein_position_and_alleles(protein_seq, protein_position_start,
                        protein_position_end, reference_aa, observed_aa)
            except InvalidVariantException as ive:
                logging.getLogger(__name__).error("Could not properly adjust protein position for variant: %s, %s, %s, %s, %s VT: %s" % (tx.get_contig(), start, end, ref_allele, alt_allele, variant_type))
                logging.getLogger(__name__).error(str(ive))
                logging.getLogger(__name__).warn("Above error may not have exact start and end positions if this is a VCF input.")
                logging.getLogger(__name__).warn("Variant type is likely incorrect.  This can happen with some GATK VCFs")
                logging.getLogger(__name__).warn(TranscriptProviderUtils.is_valid_xNP(variant_type, ref_allele, alt_allele))
                logging.getLogger(__name__).warn("The protein_change annotation may not be properly rendered.")

        vc_tmp, vc_tmp_secondary = self.infer_variant_classification(variant_type, reference_aa, observed_aa, ref_allele, alt_allele,
                                                   is_frameshift_indel=is_frameshift_indel, is_splice_site=is_splice_site, is_start_codon=is_start_codon)

        cds_start_exon_space, cds_end_exon_space = TranscriptProviderUtils.determine_cds_in_exon_space(tx)
        exon_i = TranscriptProviderUtils.determine_exon_index(int(start), int(end), tx, variant_type)
        final_vc = VariantClassification(vc_tmp, variant_type, transcript_id=tx.get_transcript_id(), alt_codon=mutated_codon_seq, ref_codon=reference_codon_seq, ref_aa=reference_aa, ref_protein_start=protein_position_start, ref_protein_end=protein_position_end, alt_aa=observed_aa, alt_codon_start_in_exon=cds_codon_start, alt_codon_end_in_exon=cds_codon_end, ref_codon_start_in_exon=cds_codon_start, ref_codon_end_in_exon=cds_codon_end, cds_start_in_exon_space=cds_start_exon_space, ref_allele_stranded=reference_allele_stranded, alt_allele_stranded=observed_allele_stranded, exon_i=exon_i, vc_secondary=vc_tmp_secondary)
        return final_vc

    def _determine_beyond_exon_info_vt(self, start, end, tx, variant_type):
        if variant_type != VariantClassification.VT_INS:
            is_beyond_exons, side, is_flank = self._determine_beyond_exon_info(int(start), int(end), tx)
        else:
            is_beyond_exons, side, is_flank = self._determine_beyond_exon_info(int(start), int(start), tx)
        return is_beyond_exons, side, is_flank

    def _determine_if_cds_overlap(self, s, e, tx, variant_type):
        if variant_type == VariantClassification.VT_INS:
            is_cds_overlap = TranscriptProviderUtils.test_feature_overlap(s, s, tx.get_cds()) != -1
        else:
            is_cds_overlap = TranscriptProviderUtils.test_feature_overlap(s, e, tx.get_cds()) != -1
        return is_cds_overlap

    def _determine_codon_overlap(self, s, e, codon_tuple, variant_type):
        if codon_tuple is None:
            return False
        if variant_type == VariantClassification.VT_INS:
            is_codon_overlap = TranscriptProviderUtils.test_overlap(s, s, codon_tuple[0]+1, codon_tuple[1])
        else:
            is_codon_overlap = TranscriptProviderUtils.test_overlap(s, e, codon_tuple[0]+1, codon_tuple[1])
        return is_codon_overlap

    def variant_classify(self, tx, ref_allele, alt_allele, start, end, variant_type, dist=2):
        """Perform classifications.

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
            if gene_type == VariantClassification.LINCRNA:
                return VariantClassification(VariantClassification.LINCRNA, variant_type, tx.get_transcript_id())
            else:
                return VariantClassification(VariantClassification.RNA, variant_type, tx.get_transcript_id())

        if ref_allele == "-":
            ref_allele = ""
        if alt_allele == "-":
            alt_allele = ""

        s = int(start)
        e = int(end)
        is_exon_overlap = TranscriptProviderUtils.determine_if_exon_overlap(s, e, tx, variant_type)

        is_splice_site_tuple = self._determine_if_splice_site_overlap(s, e, tx, variant_type, dist)
        is_splice_site = is_splice_site_tuple[0]

        is_beyond_exons, side, is_flank = self._determine_beyond_exon_info_vt(start, end, tx, variant_type)

        if not is_exon_overlap and not is_beyond_exons:
            exon_i = TranscriptProviderUtils.determine_closest_exon(tx, int(start), int(end))
            if is_splice_site:
                # Intron Splice Site
                return VariantClassification(VariantClassification.SPLICE_SITE, variant_type, tx.get_transcript_id(), vc_secondary=VariantClassification.INTRON, exon_i=exon_i)
            else:
                return VariantClassification(VariantClassification.INTRON, variant_type, tx.get_transcript_id(), exon_i=exon_i)

        if not is_exon_overlap and is_beyond_exons:
            if is_flank:
                # Flanks
                if side.startswith("3"):
                    return VariantClassification(VariantClassification.THREE_PRIME_PRIME_FLANK, variant_type, transcript_id=tx.get_transcript_id())
                else:
                    return VariantClassification(VariantClassification.FIVE_PRIME_PRIME_FLANK, variant_type, transcript_id=tx.get_transcript_id())

            else:
                # IGR
                return VariantClassification(VariantClassification.IGR, variant_type)

        is_start_codon_overlap = self._determine_codon_overlap(s, e, tx.get_start_codon(), variant_type)
        is_stop_codon_overlap = self._determine_codon_overlap(s, e, tx.get_stop_codon(), variant_type)

        if is_start_codon_overlap and not variant_type.endswith("NP"):
            return VariantClassification('Start_Codon_' + variant_type.capitalize(), variant_type, transcript_id=tx.get_transcript_id())
        if is_stop_codon_overlap and not variant_type.endswith("NP"):
            return VariantClassification('Stop_Codon_' + variant_type.capitalize(), variant_type, transcript_id=tx.get_transcript_id())

        is_cds_overlap = self._determine_if_cds_overlap(s, e, tx, variant_type)
        if is_exon_overlap and not is_cds_overlap and not is_start_codon_overlap and not is_stop_codon_overlap:
            # UTR
            if side.startswith("3"):
                vc_tmp = VariantClassification.THREE_PRIME_UTR
            else:
                vc_tmp = VariantClassification.FIVE_PRIME_UTR
            transcript_position_exon_space_start, transcript_position_exon_space_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(start, end, tx)
            vc = self._determine_de_novo(vc_tmp, transcript_position_exon_space_start, ref_allele, alt_allele, tx, variant_type)
            return VariantClassification(vc, variant_type, transcript_id=tx.get_transcript_id(), )

        # We have a clean overlap in the CDS.  Includes start codon or stop codon.
        if is_cds_overlap or is_stop_codon_overlap or is_start_codon_overlap:
            is_frameshift_indel = self.is_frameshift_indel(variant_type, int(start), int(end), alt_allele)
            return self._determine_vc_for_cds_overlap(start, end, ref_allele, alt_allele, is_frameshift_indel, is_splice_site, tx, variant_type, is_start_codon_overlap)

        raise ValueError("Could not determine variant classification:  " + tx.get_trancript_id() + " " + str([ref_allele, alt_allele, start, end]))

    def _determine_beyond_exon_info(self, start, end, tx, flank_padding_5prime=None, flank_padding_3prime=None):
        """
        start and end must be completely non-overlapping of the tx exons.  Reminder that tx exons include UTRs.
        Also, this will return False, closest side (3' or 5'), False if Intron.

        :param start: int
        :param end: int
        :param tx: Transcript
        :return: triplet of is_beyond_exons, side (always a str of 3' or 5' or "" if not is_beyond_exons), is_flank
        """
        if flank_padding_5prime is None:
            flank_padding_5prime = VariantClassifier.FLANK_PADDING_5PRIME
        if flank_padding_3prime is None:
            flank_padding_3prime = VariantClassifier.FLANK_PADDING_3PRIME

        tx_start = tx.get_start()
        tx_end = tx.get_end()

        # TODO: This correction is needed here, but not sure if this is an issue with bcbio gff parsing (issue #140)
        if tx.get_strand() == "+":
            tx_start += 1

        is_beyond_exons_left = (start < tx_start and end < tx_start)
        is_beyond_exons_right = (start > tx_end and end > tx_end)
        is_beyond_exons = is_beyond_exons_left or is_beyond_exons_right
        side = self._determine_strand_side(start, end, tx)
        if not is_beyond_exons:
            return False, side, False
        if is_beyond_exons_left:
            d = min(abs(start - tx_start), abs(end - tx_start))
        else:
            d = min(abs(start - tx_end), abs(end - tx_end))

        if side == "5'":
            is_flank = ((d <= flank_padding_5prime) and (flank_padding_5prime > 0))
        else:
            is_flank = ((d <= flank_padding_3prime) and (flank_padding_3prime > 0))
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
        cds_start, cds_end = tx.determine_cds_footprint()
        strand = tx.get_strand()

        is_beyond_cds_left = (start < cds_start and end < cds_start)
        is_beyond_cds_right = (start > cds_end and end > cds_end)
        if is_beyond_cds_left and strand == "+":
            return "5'"
        if is_beyond_cds_left and strand == "-":
            return "3'"
        if is_beyond_cds_right and strand == "+":
            return "3'"
        if is_beyond_cds_right and strand == "-":
            return "5'"

        d_left = min(abs(start - cds_start), abs(end - cds_start))
        d_right = min(abs(start - cds_end), abs(end - cds_end))

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

    def _mutate_exon(self, tx, ref, alt, vt, exon_start, buffer):
        """
        :param tx:
        :param vc: an instance of VariantClassification
        :param ref: (str) not stranded
        :param alt:  (str) not stranded
        :param vt: variant type
        :param exon_start:  an exon start position in exon space
        :return:
        """
        if alt == "-":
            alt = ""

        observed_alt = self._determine_stranded_allele(alt, tx.get_strand())

        # Inject the alternate allele
        if vt == VariantClassification.VT_INS:
            exon_end = exon_start
            relevant_seq_with_buffer = tx.get_seq()[exon_start-buffer+1:exon_end+buffer+1]
            mutated_seq = relevant_seq_with_buffer[:buffer] + observed_alt + relevant_seq_with_buffer[buffer:]
        elif vt == VariantClassification.VT_DEL:
            buffered_seq = tx.get_seq()[exon_start - len(ref) - buffer + 1: exon_start + 1 + buffer]
            mutated_seq = buffered_seq[0:buffer] + buffered_seq[-buffer:]
        else:
            exon_end = exon_start + len(alt) - 1
            relevant_seq_with_buffer = tx.get_seq()[exon_start-buffer:exon_end+buffer+1]
            mutated_seq = relevant_seq_with_buffer[0:buffer] + observed_alt + relevant_seq_with_buffer[buffer + len(observed_alt):len(relevant_seq_with_buffer)]

        return mutated_seq

    def _determine_de_novo(self, vc_str, exon_start, ref, alt, tx, variant_type, buffer=2 ):
        """Returns input vc if not de Novo.  Otherwise, returns updated variant classification.

        :param exon_start:
        :param buffer:
        :param vc_str: Current variant classification.  Note that if this is not 5'UTR, this method will just return this input.
        :param ref: (str) Does not take into account strandedness (e.g. m.ref_allele)
        :param alt: (str) Does not take into account strandedness (e.g. m.alt_allele)
        :param tx: transcript
        :param variant_type:
         Will always return original vc if the vc is not None."""
        result = vc_str
        if vc_str == VariantClassification.FIVE_PRIME_UTR and ref != alt:
            mutated_utr_region = self._mutate_exon(tx, ref, alt, variant_type, exon_start, buffer)
            atg_position = mutated_utr_region.find('ATG')
            if atg_position > -1:
                atg_exon_position = exon_start + atg_position - buffer
                cds_start_in_exon_space, cds_end_in_exon_space = TranscriptProviderUtils.determine_cds_in_exon_space(tx)
                if (cds_start_in_exon_space - atg_exon_position) % 3 == 0:
                    frameness = 'InFrame'
                else:
                    frameness = 'OutOfFrame'
                result = 'De_novo_Start_' + frameness

        return result

    def _determine_stranded_allele(self, unstranded_allele, strand):
        observed_allele_stranded = unstranded_allele
        if strand == '-':
            observed_allele_stranded = Bio.Seq.reverse_complement(unstranded_allele)
        return observed_allele_stranded

    def _determine_de_novo_old(self, vc, transcript_position_start, transcript_position_end, ref, alt, tx, variant_type):
        """Returns input vc if not de Novo.  Otherwise, returns updated variant classification.

        :param vc: Current variant classification.  Note that if this is not 5'UTR, this method will just return this input.
        :param transcript_position_start:
        :param transcript_position_end:
        :param ref: (str) Does not take into account strandedness (e.g. m.ref_allele)
        :param alt: (str) Does not take into account strandedness (e.g. m.alt_allele)
        :param tx: transcript
        :param variant_type:
         Will always return original vc if the vc is not None."""
        result = vc

        if vc == VariantClassification.FIVE_PRIME_UTR and ref != alt:
            observed_allele_stranded = self._determine_stranded_allele(alt, tx.get_strand())
            reference_allele_stranded = self._determine_stranded_allele(ref, tx.get_strand())
            tx_seq = tx.get_seq()

            if variant_type == VariantClassification.VT_INS:
                if tx.get_strand() == "-":
                    transcript_position_start = transcript_position_end
                else:
                    transcript_position_end = transcript_position_start
            utr_region_start, utr_region_end = transcript_position_start-2, transcript_position_end+2
            # TODO: This may not work for "+" strand.  Need unit test.
            utr_region_seq = tx_seq[utr_region_start:utr_region_end+1]

            mutated_utr_region_seq = TranscriptProviderUtils.mutate_reference_sequence(utr_region_seq, utr_region_start,
                transcript_position_start, transcript_position_end, observed_allele_stranded, variant_type)
            # Check for Denovo
            ATG_position = mutated_utr_region_seq.find('ATG')
            if ATG_position > -1:
                cds_start_in_exon_space, cds_end_in_exon_space = TranscriptProviderUtils.determine_cds_in_exon_space(tx)

                ATG_position = utr_region_start + ATG_position + 1
                if (cds_start_in_exon_space - ATG_position) % 3 == 0:
                    frameness = 'InFrame'
                else:
                    frameness = 'OutOfFrame'
                result = 'De_novo_Start_' + frameness
        return result

    def generate_protein_change_from_vc(self, vc):
        """

        :param vc: VariantClassification
        :return:
        """
        prot_position_start = vc.get_ref_protein_start()
        prot_position_end = vc.get_ref_protein_end()
        if prot_position_start == "" or prot_position_end == "":
            return ""
        ref_prot_allele = vc.get_ref_aa()
        alt_prot_allele = vc.get_alt_aa()
        result = TranscriptProviderUtils.render_protein_change(vc.get_vt(), vc.get_vc(), int(prot_position_start), int(prot_position_end), ref_prot_allele, alt_prot_allele, vc.get_secondary_vc())
        return result

    def generate_codon_change_from_vc(self, t, start, end, vc):
        """

        :param t: (Transcript)
        :param start: (int)
        :param end:  (int)
        :param vc:  (VariantClassification)

        :return:
        """
        dist_from_exon = self._get_splice_site_coordinates(t, start, end, vc.get_exon_i())
        exon_i = vc.get_exon_i()
        if vc.get_vc() == VariantClassification.SPLICE_SITE and vc.get_secondary_vc() == VariantClassification.INTRON:
            return TranscriptProviderUtils.render_intronic_splice_site_codon_change(dist_from_exon, exon_i)

        if vc.get_ref_codon_start_in_exon() == "" or vc.get_ref_codon_end_in_exon() == "":
            return ""

        codon_position_start_cds_space = int(vc.get_ref_codon_start_in_exon()) - int(vc.get_cds_start_in_exon_space())+1
        codon_position_end_cds_space = int(vc.get_ref_codon_end_in_exon()) - int(vc.get_cds_start_in_exon_space())+1

        ref_codon_seq = vc.get_ref_codon()
        alt_codon_seq = vc.get_alt_codon()

        result = TranscriptProviderUtils.render_codon_change(vc.get_vt(), vc.get_vc(), int(codon_position_start_cds_space), int(codon_position_end_cds_space), ref_codon_seq, alt_codon_seq, dist_from_exon, exon_i, vc.get_secondary_vc())
        return result

    def generate_transcript_change_from_tx(self, tx, variant_type, vc, start_genomic_space, end_genomic_space, ref_allele, alt_allele):
        """

        :param vc:
        :return:
        """

        if vc.get_vc() == VariantClassification.SPLICE_SITE and vc.get_secondary_vc() == VariantClassification.INTRON:
            return ""
            # dist_from_exon = self._get_splice_site_coordinates(tx, start_genomic_space, end_genomic_space, vc.get_exon_i())
            # exon_i = vc.get_exon_i()
            # return TranscriptProviderUtils.render_splice_site_transcript_change(tx, dist_from_exon, exon_i, vc.get_secondary_vc() == VariantClassification.INTRON)

        if vc.get_cds_start_in_exon_space() == "" or vc.get_cds_start_in_exon_space() < 0:
            return ""
        exon_position_start,exon_position_end = TranscriptProviderUtils.convert_genomic_space_to_exon_space(int(start_genomic_space), int(end_genomic_space), tx)

        if tx.get_strand() == "-":
            cds_position_start_cds_space = exon_position_start - int(vc.get_cds_start_in_exon_space())+1
            cds_position_end_cds_space = exon_position_end - int(vc.get_cds_start_in_exon_space())+1
        else:
            cds_position_start_cds_space = exon_position_start - int(vc.get_cds_start_in_exon_space())
            cds_position_end_cds_space = exon_position_end - int(vc.get_cds_start_in_exon_space())

        observed_allele_stranded, reference_allele_stranded = self._get_stranded_alleles(ref_allele, alt_allele, tx)
        result = TranscriptProviderUtils.render_transcript_change(variant_type, vc.get_vc(), cds_position_start_cds_space, cds_position_end_cds_space, reference_allele_stranded, observed_allele_stranded, vc.get_secondary_vc())
        return result


    def _get_splice_site_coordinates(self, t, start, end, exon_i):
        """Returns distance from exon."""

        left_diff, right_diff = TranscriptProviderUtils.determine_closest_distance_from_exon(start, end, exon_i,  t)

        if abs(left_diff) < abs(right_diff):
            dist_from_exon = left_diff * -1
            if dist_from_exon > -1: dist_from_exon = -1
        elif abs(right_diff) < abs(left_diff):
            dist_from_exon = right_diff * -1
            if dist_from_exon < 1: dist_from_exon = 1
        else:
            dist_from_exon = 0

        if t.get_strand() == "-":
            dist_from_exon *= -1
        return dist_from_exon

