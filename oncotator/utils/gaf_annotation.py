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


import itertools
import math
import Bio.Seq
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.utils.MutUtils import MutUtils


class GAFNonCodingTranscript(Exception):
    ####Custom Exception class for non-coding trancscipts
    def __init__(self):
        pass

#class MutReferenceTranscriptReferenceMismatch(Exception):
#    ####Custom Exception class for db lookup errors
#    def __init__(self, value):
#        self.value = value

class IntronUTRFlankMutation(Exception):
    ####Custom Exception class for db lookup errors
    def __init__(self, exon_affected, value):
        self.value = value
        self.exon_affected = exon_affected

class OutOfFrameIndel(Exception):
    def __init__(self, value):
        self.value = ''.join([value[0],value[1:].lower()])
 
class SpliceSiteMutation(Exception):
    ####Custom Exception class for non-coding trancscipts
    def __init__(self, exon_affected, cDNA_pos, dist_from_exon, prot_pos, prot_allele):
        self.exon_affected = exon_affected
        self.cDNA_pos = cDNA_pos
        self.dist_from_exon = dist_from_exon
        self.prot_pos = prot_pos
        self.prot_allele = prot_allele

class MicroRNA(Exception):
    ####Custom Exception class for non-coding trancscipts
    def __init__(self, transcript_pos_start, transcript_pos_end):
        self.transcript_pos_start = transcript_pos_start
        self.transcript_pos_end = transcript_pos_end
        
class DeNovoStart(Exception):
    ####Custom Exception class for non-coding trancscipts
    def __init__(self, denovo_frameness, transcript_position_start, transcript_position_end,
        utr_region_start, utr_region_end, utr_region_seq, mutated_utr_region_seq):
        self.denovo_frameness = denovo_frameness
        self.transcript_position_start = transcript_position_start
        self.transcript_position_end = transcript_position_end
        self.utr_region_start = utr_region_start
        self.utr_region_end = utr_region_end
        self.utr_region_seq = utr_region_seq
        self.mutated_utr_region_seq = mutated_utr_region_seq

def test_overlap(a_st, a_en, b_st, b_en, strand):
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

def find_overlapping_exon(start, end, t):
    ##Find exon that overlaps with alteration
    exon_affected = None
    overlap_type = None
    for i,g in enumerate(t['genomic_coords']):
        overlap_type = test_overlap(start, end, g[0], g[1], t['strand'])
        if overlap_type:
            exon_affected = i
            break
    return overlap_type, exon_affected
    
def find_overlapping_exon_with_splicesite_buffers(gaf, start, end, t):
    for i,g in enumerate(t['genomic_coords']):
        if i == 0:
            splice_dist_1, splice_dist_2 = 0, 2
        elif i + 1 == t['n_exons']:
            splice_dist_1, splice_dist_2 = 2, 0
        else:
            splice_dist_1, splice_dist_2 = 2, 2
        if t['strand'] == '+': h1, h2 = g[0]-splice_dist_1, g[1]+splice_dist_2
        elif t['strand'] == '-': h1, h2 = g[0]+splice_dist_2, g[1]-splice_dist_1
        overlap_type = test_overlap(start, end, h1, h2, t['strand'])
        if overlap_type:
            exon_affected = i
            cDNA_pos, dist_from_exon, prot_pos, prot_allele = get_splice_site_coordinates(gaf, t, start, end, exon_affected)
            exon_affected += 1 # +1 to account for 0-based python lists
            raise SpliceSiteMutation(exon_affected, cDNA_pos, dist_from_exon, prot_pos, prot_allele)


def determineClosestDistanceFromExon(end, exon_affected, start, t):
    left_start_diff = t['genomic_coords'][exon_affected][0] - start
    left_end_diff = t['genomic_coords'][exon_affected][0] - end
    right_start_diff = t['genomic_coords'][exon_affected][1] - start
    right_end_diff = t['genomic_coords'][exon_affected][1] - end
    left_diff = min(left_start_diff, left_end_diff)
    right_diff = max(right_start_diff, right_end_diff)
    return left_diff, right_diff


def get_splice_site_coordinates(gaf, t, start, end, exon_affected):
    """Returns tuple of cDNA_pos, distance from exon, protein position, and protein allele."""

    left_diff, right_diff = determineClosestDistanceFromExon(end, exon_affected, start, t)
    
    if abs(left_diff) < abs(right_diff):
        dist_from_exon = left_diff * -1
        if dist_from_exon > -1: dist_from_exon = -1
    elif abs(right_diff) < abs(left_diff):
        dist_from_exon = right_diff * -1
        if dist_from_exon < 1: dist_from_exon = 1
    else:
        dist_from_exon = 0
    
    if dist_from_exon < 1:
        cDNA_position = t['transcript_coords'][exon_affected][0]
    else:
        cDNA_position = t['transcript_coords'][exon_affected][1]
    
    if t['cds_start'] is None or cDNA_position < t['cds_start'] or cDNA_position > t['cds_stop']:
        prot_pos = 0
        prot_allele = ''
    else:
        prot_pos = get_protein_positions(cDNA_position, cDNA_position, t)[0]
        
        prot_seq = get_protein_sequence(t, gaf)
        prot_seq = ''.join([prot_seq, '*'])
        prot_allele = prot_seq[prot_pos-1]
    
    return cDNA_position, dist_from_exon, prot_pos, prot_allele


def infer_Intron_UTR_Flank_type(start, end, t):
    tx_overlap_type = test_overlap(start, end, t['tx_start'], t['tx_end'], t['strand'])
    if tx_overlap_type == 'a_within_b':
        vc = 'Intron'
    elif tx_overlap_type == 'a_overlaps_b_left_border':
        vc = "5'UTR"
    elif tx_overlap_type == 'a_overlaps_b_right_border':
        vc = "3'UTR"
    elif not tx_overlap_type:
        if start < t['tx_start'] and start < t['tx_end']: p = 5
        elif start > t['tx_start'] and start > t['tx_end']: p = 3
        if t['strand'] == '-': p = {5:3,3:5}[p]
        vc = "%d'Flank" % (p)
    exon_affected = ''
    raise IntronUTRFlankMutation(exon_affected, vc)

def convert_genomic_position_to_transcript_position(start, end, exon_affected, is_mirna, t):
    if t['strand'] == '+':
        idx_start = start - t['genomic_coords'][exon_affected][0]
        idx_end = end - t['genomic_coords'][exon_affected][0]
    else:
        idx_start = t['genomic_coords'][exon_affected][0] - end
        idx_end = t['genomic_coords'][exon_affected][0] - start
    transcript_position_start = t['transcript_coords'][exon_affected][0] + idx_start
    transcript_position_end = t['transcript_coords'][exon_affected][0] + idx_end
    if is_mirna:
        raise MicroRNA(transcript_position_start, transcript_position_end)
    else:
        return transcript_position_start, transcript_position_end

def get_protein_positions(transcript_position_start, transcript_position_end, t):
    protein_position_start = int(math.ceil((float(transcript_position_start) - \
        float(t['cds_start']) + 1) /3))
    if transcript_position_end == transcript_position_start:
        protein_position_end = protein_position_start
    else:
        protein_position_end = int(math.ceil((float(transcript_position_end) - \
            float(t['cds_start']) + 1) /3))
    return protein_position_start, protein_position_end

def get_cds_codon_positions(protein_position_start, protein_position_end, t):
    cds_codon_start = (protein_position_start * 3 - 2) + t['cds_start'] - 1
    cds_codon_end = (protein_position_end * 3) + t['cds_start'] - 1
    return cds_codon_start, cds_codon_end

def test_for_framshift_indel(m, start, end, reference_allele, observed_allele):
    if m['variant_type'] in ['DEL','INS']:
        if m['variant_type'] == 'DEL':
            dl = (end-start+1)%3
        elif  m['variant_type'] == 'INS':
            dl = len(observed_allele)%3
        if dl != 0:
            raise OutOfFrameIndel(m['variant_type'])
    elif len(reference_allele) != len(observed_allele):
        if abs(len(observed_allele) - len(reference_allele)) %3 != 0:
            if len(observed_allele) > len(reference_allele):
                type = 'INS'
            else:
                type = 'DEL'
            raise OutOfFrameIndel(type)

def mutate_reference_sequence(seq, seq_start_pos, mutation_start_pos, mutation_end_pos, observed_allele, m):
    mutated_seq = list(seq)
    if m['variant_type'] == 'DEL':
        del(mutated_seq[mutation_start_pos - seq_start_pos: \
            mutation_end_pos - seq_start_pos + 1])
    elif m['variant_type'] == 'INS':
        mutated_seq.insert(mutation_end_pos - seq_start_pos, observed_allele)
    else: #SNP or ONP
        mutated_seq[mutation_start_pos - seq_start_pos: \
            mutation_end_pos - seq_start_pos + 1] = str(observed_allele)
    mutated_seq = ''.join(mutated_seq)
    return mutated_seq


def determineIfSpliceSiteThrown(end, exon_affected, gaf, start, t, dist=2):

    """
     If overlap is detected, but the start or end is within dist bp, then this is a splice site.
    :param end:
    :param exon_affected:
    :param gaf:
    :param start:
    :param t:
    :raise: SpliceSiteMutation if the variant (start, end) should be a splice site.
    """
    ldist, rdist = determineClosestDistanceFromExon(end, exon_affected, start, t)
    if abs(ldist) < dist or abs(rdist) < dist:
        cDNA_pos, dist_from_exon, prot_pos, prot_allele = get_splice_site_coordinates(gaf, t, start, end, exon_affected)
        exon_affected += 1 # +1 to account for 0-based python lists
        raise SpliceSiteMutation(exon_affected, cDNA_pos, dist_from_exon, prot_pos, prot_allele)


def get_transcript_positions(gaf, start, end, is_mirna, t):
    if t['strand'] == '+': direction = 1
    elif t['strand'] == '-': direction = -1
    overlap_type, exon_affected = find_overlapping_exon(start, end, t)

    # This covers whether the variant overlaps the transcript, but should still be called a splice site
    if exon_affected is not None:
        determineIfSpliceSiteThrown(end, exon_affected, gaf, start, t)

    #!#If no exons found to overlap, search in splice_sites
    if exon_affected is None and t['n_exons'] > 1:
        find_overlapping_exon_with_splicesite_buffers(gaf, start, end, t)
    if exon_affected is None: #must be in intron or 5'/3' flank
        infer_Intron_UTR_Flank_type(start, end, t)
    if overlap_type == 'a_within_b':
        #transcript_position_start, transcript_position_end = \
        #    convert_genomic_position_to_transcript_position_deprecated(start, end, exon_affected, direction, is_mirna, t)
        transcript_position_start, transcript_position_end = \
            convert_genomic_position_to_transcript_position(start, end, exon_affected, is_mirna, t)
        # LTL April 9, 2013
        # exon_affected += 1
    else: # handle coding indels that overlap into adjacent introns
        if is_mirna: raise MicroRNA(None, None)
        exon_affected += 1
        if exon_affected == 1 and overlap_type == 'a_overlaps_b_left_border':
            raise IntronUTRFlankMutation(exon_affected, "5'UTR")
        elif exon_affected == t['n_exons'] and overlap_type == 'a_overlaps_b_right_border':
            raise IntronUTRFlankMutation(exon_affected, "3'UTR")
        else:
            cDNA_pos, dist_from_exon, prot_pos, prot_allele = get_splice_site_coordinates(gaf, t, start, end, exon_affected-1)
            raise SpliceSiteMutation(exon_affected, cDNA_pos, dist_from_exon, prot_pos, prot_allele)
    return transcript_position_start, transcript_position_end, exon_affected

def get_protein_sequence(t, gaf):
    protein_seq = gaf.get_protein_seq(t['transcript_id'])
    if protein_seq is None:
        raise GAFNonCodingTranscript
    else:
        return protein_seq

def is_non_coding_transcript(transcript, gaf):
    tx_seq = gaf.get_transcript_seq(transcript["transcript_id"])
    if not tx_seq:
        return True

    if 'cds_start' not in transcript or not transcript['cds_start']:
        return True

    return False
    
#    if t['transcript_id'] not in gaf.gaf_protein_sequences:
#        raise GAFNonCodingTranscript
#    else:
#        protein_seq = gaf.gaf_protein_sequences[t['transcript_id']]
#        return protein_seq

def infer_variant_classification(variant_type, reference_aa, observed_aa, reference_allele, observed_allele):
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

def annotate_mutations_not_fully_within_cds(transcript_seq, observed_allele, cds_overlap_type, transcript_position_start, transcript_position_end, m, t, exon_affected):
    if m['variant_type'] in ['DEL','INS']: vt = ''.join([m['variant_type'][0], m['variant_type'][1:].lower()])
    else: vt = m['variant_type']
    if cds_overlap_type == 'a_overlaps_b_left_border': raise IntronUTRFlankMutation(exon_affected, ''.join(['Start_Codon_', vt]))
    elif cds_overlap_type == 'a_overlaps_b_right_border': raise IntronUTRFlankMutation(exon_affected, ''.join(['Stop_Codon_', vt]))
    else:
        if transcript_position_start < t['cds_start'] and transcript_position_start < t['cds_stop']: p = 5
        elif transcript_position_start > t['cds_start'] and transcript_position_start > t['cds_stop']: p = 3
        vc = "%d'UTR" % (p)
        if vc == "5'UTR":
            utr_region_start, utr_region_end = transcript_position_start-2, transcript_position_end+2
            
            utr_region_seq = transcript_seq[utr_region_start-1:utr_region_end]
            mutated_utr_region_seq = mutate_reference_sequence(utr_region_seq, utr_region_start,
                transcript_position_start, transcript_position_end, observed_allele, m)
            
            ATG_position = mutated_utr_region_seq.find('ATG')
            if ATG_position > -1:
                ATG_position = utr_region_start + ATG_position
                if (t['cds_start'] - ATG_position) % 3 == 0:
                    frameness = 'InFrame'
                else:
                    frameness = 'OutOfFrame'
                raise DeNovoStart(frameness, transcript_position_start, transcript_position_end, utr_region_start, utr_region_end, utr_region_seq, mutated_utr_region_seq)
            
        raise IntronUTRFlankMutation(exon_affected, vc)

def chop(iterable, length=2):
    return itertools.izip(*(iter(iterable) , ) *length)

def adjust_protein_position_and_alleles(protein_seq, protein_position_start, protein_position_end, reference_aa, observed_aa):
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

def identify_best_effect_transcript(data, gaf):
    order_dict = TranscriptProviderUtils.retrieve_effect_dict()
    for m in data:
        if m['variant_type'] != 'ERR':
            
            best_transcript = None
            for i,t in enumerate(m['transcripts']):
                if t['variant_classification'] == 'ERR': continue
                if best_transcript is None:
                    best_transcript = i
                    best_rank = order_dict[t['variant_classification']]
                    try:
                        best_code_len = t['code_len']
                    except KeyError:
                        #for IGR mutations
                        best_code_len = 0
                else:
                    if order_dict[t['variant_classification']] < best_rank or (order_dict[t['variant_classification']] == best_rank and t['code_len'] > best_code_len):
                        best_transcript = i
                        best_rank = order_dict[t['variant_classification']]
                        best_code_len = t['code_len']
                        
            if best_transcript is None:
                m.createAnnotation('best_effect_transcript', 0, gaf.title)
            else:
                m.createAnnotation('best_effect_transcript', best_transcript, gaf.title)
            
        yield m


def identify_best_canonical_transcript(data, gaf):
    for m in data:
        if m['variant_type'] != 'ERR':
            genes = set()
            for t in m['transcripts']:
                try:
                    genes.add(t['gene'])
                except KeyError:
                    pass
            
            gene = ''
            if len(genes) > 1:
                # if mutation overlaps with more than one gene, gene chosen is that of the best effect transcript
                gene = m['transcripts'][m['best_effect_transcript']]['gene']
            elif len(genes) == 1:
                gene = genes.pop()
                
            if gene:
                try:
                    canonical_tx = gaf.get_gene(gene)['canonical_tx']
                except AttributeError:
                    m.createAnnotation('best_canonical_transcript', m['best_effect_transcript'], gaf.title)
                else:
                    for i,t in enumerate(m['transcripts']):
                        if t['transcript_id'] == canonical_tx:
                            m.createAnnotation('best_canonical_transcript', str(i), gaf.title)
                            
                    if 'best_canonical_transcript' not in m:
                        m.createAnnotation('best_canonical_transcript', m['best_effect_transcript'], gaf.title)
            else:
                m.createAnnotation('best_canonical_transcript', 0, gaf.title)
    
        yield m
    
def find_mut_in_gaf(data, gaf):
    for m in data:
        if m['variant_type'] != 'ERR':

            # Need to convert start and end positions into integers.
            start, end = int(m['start']), int(m['end'])
            transcripts = gaf.get_overlapping_transcripts(m['chr'],start,end)
            if not transcripts:
                transcripts.append(dict())
                transcripts[0]['variant_classification'] = 'IGR'
                nearest_genes = gaf.get_nearest_genes(m['chr'],start,end)
                if nearest_genes:
                    transcripts[0]['nearest_genes'] = '%s (%s upstream) : %s (%s downstream)' % (nearest_genes[0][0], nearest_genes[0][1], nearest_genes[1][0], nearest_genes[1][1])
                else:
                    transcripts[0]['nearest_genes'] = ''
                #transcripts[0]['code_len'] = 0
                m.createAnnotation('transcripts', transcripts, gaf.title)
            else:
                m.createAnnotation('transcripts', [dict() for i in xrange(len(transcripts))], gaf.title)
                for t_idx, t in enumerate(transcripts):
                    try:
                        reference_allele, observed_allele = str(m.ref_allele), str(m.alt_allele)
                        if t['class'] == 'miRNA': is_mirna = True
                        else: is_mirna = False
                        
                        transcript_position_start, transcript_position_end, exon_affected = get_transcript_positions(gaf, start, end, is_mirna, t)
                        
                        #transcript_seq = gaf.gaf_transcript_sequences[t['transcript_id']]
                        transcript_seq = gaf.get_transcript_seq(t['transcript_id'])
                        if t['strand'] == '-':
                            reference_allele, observed_allele = Bio.Seq.reverse_complement(reference_allele), Bio.Seq.reverse_complement(observed_allele)
                        
                        # Fix to correct reference transcript sequence when it differs from genomic reference
                        # -1 here to account for zero-base python lists    
                        if m['variant_type'] == 'SNP' and transcript_seq[transcript_position_start-1:transcript_position_end] != reference_allele:
                            new_transcript_seq = list(transcript_seq)
                            new_transcript_seq[transcript_position_start-1:transcript_position_end] = reference_allele
                            transcript_seq = ''.join(new_transcript_seq)
                            ref_tx_seq_has_been_changed = True
                        else:
                            ref_tx_seq_has_been_changed = False
                            
                        #if 'dbSNP_RS' not in m: #if dbsnp site(s), don't check for TX-REF-MISMATCH
                        #    if m['variant_type'] != 'INS' and transcript_seq[transcript_position_start-1:transcript_position_end] != reference_allele: #-1 here to account for zero-base python lists
                        #        #pdb.set_trace()
                        #        raise MutReferenceTranscriptReferenceMismatch(transcript_seq[transcript_position_start-1:transcript_position_end])
                        
                        protein_seq = get_protein_sequence(t, gaf)
                        protein_seq = ''.join([protein_seq, '*'])
        
                        cds_overlap_type = test_overlap(transcript_position_start, transcript_position_end,
                            t['cds_start'], t['cds_stop'], '+')  #always use '+' here because strand doesn't matter
                        
                        if cds_overlap_type == 'a_within_b':
                            protein_position_start, protein_position_end = get_protein_positions(transcript_position_start,
                                transcript_position_end, t)
                            
                            cds_codon_start, cds_codon_end = get_cds_codon_positions(protein_position_start,
                                protein_position_end, t)
                                
                            reference_codon_seq = transcript_seq[cds_codon_start-1:cds_codon_end]
                            reference_aa = protein_seq[protein_position_start-1:protein_position_end]
                            
                            test_for_framshift_indel(m, start, end, reference_allele, observed_allele)
                            
                            if m['variant_type'] == 'INS' and protein_position_start != protein_position_end:
                                #treat differently if insertion falls between codons
                                reference_codon_seq = ''
                                mutated_codon_seq = observed_allele
                                reference_aa = ''
                            else:
                                mutated_codon_seq = mutate_reference_sequence(reference_codon_seq,
                                    cds_codon_start, transcript_position_start, transcript_position_end, observed_allele, m)
                            
                            if ref_tx_seq_has_been_changed:
                                reference_aa = MutUtils.translate_sequence(reference_codon_seq)
                            
                            observed_aa = MutUtils.translate_sequence(mutated_codon_seq)
                            
                            variant_classification = infer_variant_classification(m['variant_type'],
                                reference_aa, observed_aa, reference_allele, observed_allele)

                            # If silent mutation w/in 2 bp of a splice junction, then change to splice site
                            if variant_classification.lower() == "silent":
                                determineIfSpliceSiteThrown(end, exon_affected,gaf,start,t, dist=2)

                            if m['variant_type'] != 'SNP':
                                reference_aa, observed_aa, protein_position_start, protein_position_end = \
                                    adjust_protein_position_and_alleles(protein_seq, protein_position_start,
                                        protein_position_end, reference_aa, observed_aa)
                        else:
                            annotate_mutations_not_fully_within_cds(transcript_seq, observed_allele, cds_overlap_type, transcript_position_start,
                                transcript_position_end, m, t, exon_affected)
                        
#                    except MutReferenceTranscriptReferenceMismatch as e:
#                        m['transcripts'][t_idx]['variant_classification'] = 'TX-REF-MISMATCH'
#                        m['transcripts'][t_idx]['error_message'] = 'Reference allele ("%s") of variant does not match transcript sequence ("%s")' % (reference_allele,e.value)
#                        m['transcripts'][t_idx]['reference_transcript_allele'] = reference_allele
#                        m['transcripts'][t_idx]['observed_transcript_allele'] = observed_allele
#                        m['transcripts'][t_idx]['transcript_position_start'] = transcript_position_start
#                        m['transcripts'][t_idx]['transcript_position_end'] = transcript_position_end
#                        m['transcripts'][t_idx]['exon_affected'] = exon_affected
                    except GAFNonCodingTranscript:
                        m['transcripts'][t_idx]['variant_classification'] = 'Non-coding_Transcript'
                        m['transcripts'][t_idx]['reference_transcript_allele'] = reference_allele
                        m['transcripts'][t_idx]['observed_transcript_allele'] = observed_allele
                        m['transcripts'][t_idx]['transcript_position_start'] = transcript_position_start
                        m['transcripts'][t_idx]['transcript_position_end'] = transcript_position_end
                        m['transcripts'][t_idx]['exon_affected'] = exon_affected
                    except IntronUTRFlankMutation as e:
                        m['transcripts'][t_idx]['variant_classification'] = e.value
                        m['transcripts'][t_idx]['exon_affected'] = e.exon_affected
                    except OutOfFrameIndel as e:
                        m['transcripts'][t_idx]['variant_classification'] = ''.join(['Frame_Shift_', e.value])
                        m['transcripts'][t_idx]['reference_transcript_allele'] = reference_allele
                        m['transcripts'][t_idx]['observed_transcript_allele'] = observed_allele
                        m['transcripts'][t_idx]['transcript_position_start'] = transcript_position_start
                        m['transcripts'][t_idx]['transcript_position_end'] = transcript_position_end
                        m['transcripts'][t_idx]['protein_position_start'] = protein_position_start
                        m['transcripts'][t_idx]['protein_position_end'] = protein_position_end
                        m['transcripts'][t_idx]['cds_codon_start'] = cds_codon_start
                        m['transcripts'][t_idx]['cds_codon_end'] = cds_codon_end
                        m['transcripts'][t_idx]['reference_codon_seq'] = reference_codon_seq
                        if len(reference_aa) == 0:
                            m['transcripts'][t_idx]['reference_protein_allele'] = ""
                        else:
                            m['transcripts'][t_idx]['reference_protein_allele'] = reference_aa[0]
                        m['transcripts'][t_idx]['exon_affected'] = exon_affected
                    except SpliceSiteMutation as e:
                        if m['variant_type'] in ['DEL','INS']: vt = m['variant_type'][0] + m['variant_type'][1:].lower()
                        else: vt = m['variant_type']
                        m['transcripts'][t_idx]['variant_classification'] = 'Splice_Site' #''.join(['Splice_Site_', vt])
                        m['transcripts'][t_idx]['exon_affected'] = e.exon_affected
                        m['transcripts'][t_idx]['transcript_position_start'] = e.cDNA_pos
                        m['transcripts'][t_idx]['transcript_position_end'] = e.cDNA_pos
                        m['transcripts'][t_idx]['dist_from_exon'] = e.dist_from_exon
                        m['transcripts'][t_idx]['prot_pos'] = e.prot_pos
                        m['transcripts'][t_idx]['prot_allele'] = e.prot_allele
                    except MicroRNA as e:
                        m['transcripts'][t_idx]['variant_classification'] = 'miRNA'
                        m['transcripts'][t_idx]['reference_transcript_allele'] = reference_allele
                        m['transcripts'][t_idx]['observed_transcript_allele'] = observed_allele
                        m['transcripts'][t_idx]['transcript_position_start'] = e.transcript_pos_start
                        m['transcripts'][t_idx]['transcript_position_end'] = e.transcript_pos_end
                    except DeNovoStart as e:
                        m['transcripts'][t_idx]['variant_classification'] = ''.join(['De_novo_Start_', e.denovo_frameness])
                        m['transcripts'][t_idx]['reference_transcript_allele'] = reference_allele
                        m['transcripts'][t_idx]['observed_transcript_allele'] = observed_allele
                        m['transcripts'][t_idx]['transcript_position_start'] = e.transcript_position_start
                        m['transcripts'][t_idx]['transcript_position_end'] = e.transcript_position_end
                        m['transcripts'][t_idx]['cds_codon_start'] = e.utr_region_start
                        m['transcripts'][t_idx]['cds_codon_end'] = e.utr_region_end
                        m['transcripts'][t_idx]['reference_codon_seq'] = e.utr_region_seq
                        m['transcripts'][t_idx]['mutated_codon_seq'] = e.mutated_utr_region_seq
                        m['transcripts'][t_idx]['exon_affected'] = exon_affected
                    else:
                        m['transcripts'][t_idx]['variant_classification'] = variant_classification
                        m['transcripts'][t_idx]['reference_transcript_allele'] = reference_allele
                        m['transcripts'][t_idx]['observed_transcript_allele'] = observed_allele
                        m['transcripts'][t_idx]['transcript_position_start'] = transcript_position_start
                        m['transcripts'][t_idx]['transcript_position_end'] = transcript_position_end
                        m['transcripts'][t_idx]['protein_position_start'] = protein_position_start
                        m['transcripts'][t_idx]['protein_position_end'] = protein_position_end
                        m['transcripts'][t_idx]['cds_codon_start'] = cds_codon_start
                        m['transcripts'][t_idx]['cds_codon_end'] = cds_codon_end
                        m['transcripts'][t_idx]['reference_codon_seq'] = reference_codon_seq
                        m['transcripts'][t_idx]['mutated_codon_seq'] = mutated_codon_seq
                        m['transcripts'][t_idx]['reference_protein_allele'] = reference_aa
                        m['transcripts'][t_idx]['observed_protein_allele'] = observed_aa
                        m['transcripts'][t_idx]['exon_affected'] = exon_affected
                    finally:
                        m['transcripts'][t_idx]['transcript_id'] = t['transcript_id']
                        m['transcripts'][t_idx]['gene'] = t['gene']
                        m['transcripts'][t_idx]['strand'] = t['strand']
                        m['transcripts'][t_idx]['code_len'] = t.get('code_len',0)
                        m['transcripts'][t_idx]['id'] = t.get('id',0)
                        if not is_mirna:
                            if 'refseq_mRNA_id' in t:
                                m['transcripts'][t_idx]['refseq_mRNA_id'] = t['refseq_mRNA_id']
                            if 'refseq_prot_id' in t:
                                m['transcripts'][t_idx]['refseq_prot_id'] = t['refseq_prot_id']
                            if 'description' in t:
                                m['transcripts'][t_idx]['description'] = t['description']
                            if 'tx_start' in t:
                                m['transcripts'][t_idx]['tx_start'] = t['tx_start']
                            if 'tx_end' in t:
                                m['transcripts'][t_idx]['tx_end'] = t['tx_end']
                            try:    
                                gene = gaf.get_gene(t['gene'])
                            except AttributeError:
                                pass
                            else:
                                if 'uniprot_accession' in gene:
                                    m['transcripts'][t_idx]['uniprot_accession'] = gene['uniprot_accession']
                                if 'uniprot_entry_name' in gene:
                                    m['transcripts'][t_idx]['uniprot_entry_name'] = gene['uniprot_entry_name']
                            
                m['transcripts'] = tuple(m['transcripts'])
        yield m
        
def correct_transcript_coordinates(data, gaf):
    for d in data:
        
        if 'transcripts' in d:
            for t in d['transcripts']:
                if 'transcript_id' in t:
                    tx_id = t['transcript_id']
                    transcript = gaf.get_transcript(tx_id)
                    
                    if not t['code_len']: continue
    
                    adjust_coord = lambda c: c - transcript['cds_start'] + 1
                    
                    for f in ['cds_codon_start', 'cds_codon_end', 'transcript_position_start', 'transcript_position_end']:
                        try:
                            if f.startswith('transcript_'):
                                t['absolute_' + f] = t[f]
                            t[f] = adjust_coord(t[f])
                        except KeyError:
                            continue
        
        yield d

def infer_output_fields(data, gaf):
    for i,m in enumerate(data):
        if m['variant_type'] != 'ERR':

            genome_change = TranscriptProviderUtils.determine_genome_change(m.chr, m.start, m.end, m.ref_allele, m.alt_allele, m['variant_type'])
            m.createAnnotation('genome_change', genome_change, gaf.title)
            
            for t in m['transcripts']:
                transcript_change, codon_change, protein_change = '','',''
                if t['variant_classification'].startswith('Splice_Site'):
                    transcript_change = 'c.%d_splice' % (t['transcript_position_start'])
                    
                    if t['dist_from_exon'] < 0:
                        dist_from_exon = str(t['dist_from_exon'])
                    else:
                        dist_from_exon = ''.join(['+' ,str(t['dist_from_exon'])])
                    codon_change = 'c.e%d%s' % (t['exon_affected'], dist_from_exon)
                    
                    if t['prot_pos'] > 0:
                        protein_change = 'p.%s%d_splice' % (t['prot_allele'], t['prot_pos'])
                        
                    
                    t['transcript_change'], t['codon_change'], t['protein_change'] = transcript_change, codon_change, protein_change
                    continue
                
                if 'transcript_position_start' in t:
                    try:
                        reference_allele, observed_allele = t['reference_transcript_allele'], t['observed_transcript_allele']
                        if m['variant_type'] == 'SNP':
                            transcript_change = 'c.%d%s>%s' % (t['transcript_position_start'], reference_allele, observed_allele)
                        elif m['variant_type'].endswith('NP'):
                            transcript_change = 'c.%d_%d%s>%s' % (t['transcript_position_start'],
                                t['transcript_position_end'], reference_allele, observed_allele)
                        elif m['variant_type'] == 'DEL':
                            if t['transcript_position_start'] == t['transcript_position_end']:
                                transcript_change = 'c.%ddel%s' % (t['transcript_position_start'], reference_allele)                            
                            else:
                                transcript_change = 'c.%d_%ddel%s' % (t['transcript_position_start'],
                                    t['transcript_position_end'], reference_allele)
                        elif m['variant_type'] == 'INS':
                            transcript_change = 'c.%d_%dins%s' % (t['transcript_position_start'],
                                t['transcript_position_end'], observed_allele)
                    except TypeError:
                        #exception handling for transcript positions with value 'None'
                        #e.g. insertion occurring between last position of transcript and 3' flank
                        pass
                
                if 'cds_codon_start' in t:
                    if t['variant_classification'].startswith('Frame_Shift'):
                        codon_change = 'c.(%d-%d)%sfs' % (t['cds_codon_start'], t['cds_codon_end'],
                            t['reference_codon_seq'])
                    elif m['variant_type'].endswith('NP'):
                        codon_change = 'c.(%d-%d)%s>%s' % (t['cds_codon_start'], t['cds_codon_end'],
                            t['reference_codon_seq'], t['mutated_codon_seq'])
                    elif m['variant_type'] == 'DEL':
                        if t['mutated_codon_seq'] == '': #full codon deleted
                            codon_change = 'c.(%d-%d)%sdel' % (t['cds_codon_start'], t['cds_codon_end'],
                                t['reference_codon_seq'])
                        else:
                            codon_change = 'c.(%d-%d)%s>%s' % (t['cds_codon_start'], t['cds_codon_end'],
                                t['reference_codon_seq'], t['mutated_codon_seq'])
                    elif m['variant_type'] == 'INS':
                        if t['reference_codon_seq'] == '': #insertion between codons
                            codon_change = 'c.(%d-%d)ins%s' % (t['cds_codon_start'], t['cds_codon_end'],
                                t['mutated_codon_seq'])
                        else:
                            codon_change = 'c.(%d-%d)%s>%s' % (t['cds_codon_start'], t['cds_codon_end'],
                                t['reference_codon_seq'], t['mutated_codon_seq'])

                if 'protein_position_start' in t:
                    if t['variant_classification'].startswith('Frame_Shift'):
                        protein_change = 'p.%s%dfs' % (t['reference_protein_allele'], t['protein_position_start'])
                    elif t['observed_protein_allele'] == '-':
                        protein_change = 'p.%s%ddel' % (t['reference_protein_allele'], t['protein_position_start'])
                    elif t['reference_protein_allele'] == '-':
                        protein_change = 'p.%d_%dins%s' % (t['protein_position_start'],
                            t['protein_position_end'], t['observed_protein_allele'])
                    elif len(t['reference_protein_allele']) == 1 and len(t['observed_protein_allele']) == 1:
                        protein_change = 'p.%s%d%s' % (t['reference_protein_allele'], t['protein_position_start'],
                            t['observed_protein_allele'])
                    else:
                        protein_change = 'p.%d_%d%s>%s' % (t['protein_position_start'], t['protein_position_end'],
                            t['reference_protein_allele'], t['observed_protein_allele'])
                
                t['transcript_change'], t['codon_change'], t['protein_change'] = transcript_change, codon_change, protein_change
                
                for f in ['cds_codon_start','cds_codon_end','mutated_codon_seq','reference_codon_seq']:
                    try:
                        del(t[f])
                    except KeyError:
                        continue
                
            if m['variant_type'] == 'ONP' and len(m.ref_allele) > len(m.alt_allele):
                m['variant_type'] = 'DEL'
                for t in m['transcripts']:
                    t['variant_classification'] = t['variant_classification'].replace('ONP','Del')
            elif m['variant_type'] == 'ONP' and len(m.ref_allele) < len(m.alt_allele):
                m['variant_type'] = 'INS'
                for t in m['transcripts']:
                    t['variant_classification'] = t['variant_classification'].replace('ONP','Ins')
                
        yield m
