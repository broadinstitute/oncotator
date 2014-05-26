import re
import logging

from Bio import Seq

from oncotator.datasources.ChangeTransformingDatasource import ChangeTransformingDatasource
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.utils.VariantClassifier import VariantClassifier
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils

AA_NAME_MAPPING = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val',
    '*': '*',
}

PROT_REGEXP = re.compile('p\.([A-Z*])(\d+)([A-Z*])')
PROT_FRAMESHIFT_REGEXP = re.compile('p\.([A-Z*-]+)(\d+)fs')
PROT_INFRAME_DEL_REGEXP = re.compile('p\.([A-Z*]+)(\d+)del')
PROT_INFRAME_DEL_REGEXP_2 = re.compile('p\.(\d+)_(\d+)([A-Z*]+)>([A-Z*]+)')
PROT_INFRAME_INS_REGEXP_1 = re.compile('p\.(\d+)_(\d+)ins([A-Z*]+)')
PROT_INFRAME_INS_REGEXP_2 = re.compile('p\.(\d+)_(\d+)([A-Z*]+)>([A-Z*]+)')
PROT_ONP_REGEXP = re.compile('p.(\d+)_(\d+)([A-Z*]+)>([A-Z*]+)')

class HgvsChangeTransformingDatasource(ChangeTransformingDatasource):

    def __init__(self, src_file, title='', version=None):
        super(ChangeTransformingDatasource, self).__init__(src_file, title=title, version=version)
        self.output_headers = ['HGVS_genomic_change', 'HGVS_coding_DNA_change', 'HGVS_protein_change']

        self.vcer = VariantClassifier()

        #### HACK #####
        self.gencode_ds = EnsemblTranscriptDatasource(src_file)

        #### HACK #####
        ensembl_id_mapping_fname = '/Users/aramos/Desktop/ot_hgvs/ensembl_id_mappingsGRCh37.p13.txt'
        import pandas as pd
        self.ensembl_id_mapping = dict()
        with open(ensembl_id_mapping_fname) as in_fh:
            for line in in_fh:
                tx_id, prot_id = line.split('\t')[1:]
                self.ensembl_id_mapping[tx_id] = prot_id.replace('\n', '')
        
    def annotate_mutation(self, mutation):
        # TODO
        # TODO: When annotating, you should probably give the annotation a name like:  self.title + "_" + output_annotation_name
        # You can assume that the annotations from the transript datasource have already been populated (such as protein_change and genome_change)

        mutation['start'], mutation['end'] = int(mutation['start']), int(mutation['end'])

        hgvs_genomic_change = self._adjust_genome_change(mutation)
        hgvs_coding_dna_change = self._adjust_coding_DNA_change(mutation)
        hgvs_protein_change = self._adjust_protein_change(mutation)

        mutation.createAnnotation('HGVS_genomic_change', hgvs_genomic_change, self.title)
        mutation.createAnnotation('HGVS_coding_DNA_change', hgvs_coding_dna_change, self.title)
        mutation.createAnnotation('HGVS_protein_change', hgvs_protein_change, self.title)

        return mutation

    def _adjust_genome_change(self, mutation):
        variant_type = mutation['variant_type']
        mut_start, mut_end = mutation['start'], mutation['end']
        chrn = mutation['chr']
        alt_allele = mutation['alt_allele']
        ref_allele = mutation['ref_allele']
        build = mutation['build']
        ref_context = mutation.get('ref_context')

        if variant_type == 'DNP':
            change = '%d_%ddelins%s' % (mut_start, mut_end, alt_allele)
        elif variant_type == 'INS':
            genome_pos_adjust = self._genome_pos_adjust_if_duplication(alt_allele, ref_context, variant_type)
            genome_pos_adjust = genome_pos_adjust - len(alt_allele)

            #BUG: duplications are only rendered if position needs to be adjusted! Implement duplication test like done for coding and protein position
            if self._determine_if_genomic_duplication(alt_allele, ref_context, variant_type, mutation):
                start_pos = mut_end + genome_pos_adjust
                end_pos = start_pos + len(alt_allele) - 1
                change = '%d_%ddup%s' % (start_pos, end_pos, alt_allele)
            else:
                change = '%d_%dins%s' % (mut_start, mut_end, alt_allele)
        elif variant_type == 'DEL':
            genome_pos_adjust = self._genome_pos_adjust_if_duplication(ref_allele, ref_context, variant_type)
            genome_pos_adjust = genome_pos_adjust - len(ref_allele) if genome_pos_adjust > 0 else genome_pos_adjust

            mut_start_changed = mut_start
            mut_start, mut_end = mut_start + genome_pos_adjust, mut_end + genome_pos_adjust
            if mut_start == mut_end:
                change = str(mut_start_changed) if mut_start == mut_end else '%d_%d' % (mut_start, mut_end)
                change = '%sdel%s' % (change, ref_allele)
            else:
                change = '%d_%ddel%s' % (mut_start, mut_end, ref_allele)
        else:
            change = '%d%s>%s' % (mut_start, ref_allele, alt_allele)
        
        if build == '':
            adjusted_genome_change = '%s:g.%s' % (chrn, change)
        else:
            if build == 'hg19' and not chrn.startswith('chr'):
                chrn = 'chr' + chrn
            adjusted_genome_change = '%s.%s:g.%s' % (chrn, build, change)

        return adjusted_genome_change

    def _adjust_coding_DNA_change(self, mutation):
        ### IF DELETION, THEN DUPLICATION CHECK AND POS ADJUSTMENT NOT BEING DONE!!
        if mutation['variant_type'] == 'DNP':
            return self._get_cdna_change_for_ONP(mutation)
        else:
            vc = mutation['variant_classification']
            if vc in ['Intron', 'Splice_Site']:
                return self._get_cdna_change_for_intron(mutation)
            elif vc in ["5'UTR", 'De_novo_Start_OutOfFrame']:
                return self._get_cdna_change_for_5_utr(mutation)
            elif vc == "3'UTR":
                return self._get_cdna_change_for_3_utr(mutation)
            elif vc == 'IGR':
                return ''
            elif vc.endswith('Ins'):
                return self._get_cdna_change_for_ins(mutation)
            elif vc in ['Stop_Codon_Del']:
                #need to make cdna str from scratch
                tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
                tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
                tx_stop_codon_genomic_coords = tx.get_stop_codon()
                tx_stop_start_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(tx_stop_codon_genomic_coords[0], tx_stop_codon_genomic_coords[1], tx)[0]

                tx_stop_start_pos += 1
                if mutation['transcript_strand'] == '+':
                    genome_stop_start_pos = tx_stop_codon_genomic_coords[0] + 1
                    tx_variant_pos_1_offset = mutation['start'] - genome_stop_start_pos
                    tx_variant_pos_2_offset = mutation['end'] - genome_stop_start_pos
                else:
                    genome_stop_start_pos = tx_stop_codon_genomic_coords[1]
                    tx_variant_pos_1_offset = genome_stop_start_pos - mutation['end']
                    tx_variant_pos_2_offset = genome_stop_start_pos - mutation['start']

                tx_variant_pos_1 = tx_stop_start_pos + tx_variant_pos_1_offset
                tx_variant_pos_2 = tx_stop_start_pos + tx_variant_pos_2_offset

                if tx_variant_pos_2 > tx_stop_start_pos + 2:
                    tx_variant_pos_2 = '*%d' % (tx_variant_pos_2 - (tx_stop_start_pos + 2))
                else:
                    tx_variant_pos_2 = str(tx_variant_pos_2)

                if mutation['start'] == mutation['end']:
                    change_str = 'c.%ddel%s' % (tx_variant_pos_1, tx_ref_allele)
                else:
                    change_str = 'c.%d_%sdel%s' % (tx_variant_pos_1, tx_variant_pos_2, tx_ref_allele)
                    
                return '%s:%s' % (mutation['annotation_transcript'], change_str)
            elif vc.endswith('Del'):
                return self._get_cdna_change_for_del(mutation)
            else:
                # just use the transcript change from the TranscriptProvider
                return '%s:%s' % (mutation['annotation_transcript'], mutation['transcript_change'])

    def _adjust_protein_change(self, mutation):
        vc = mutation['variant_classification']
        if vc in ['Intron', 'IGR', "3'UTR", "5'UTR", 'RNA',
            'lincRNA', 'Silent', 'Splice_Site', 'De_novo_Start_OutOfFrame', 'De_novo_Start_InFrame',
            'Start_Codon_Ins', 'Start_Codon_Del']:
            return ''
        elif vc in ['Frame_Shift_Ins', 'Frame_Shift_Del']:
            adjusted_prot_change = self._get_prot_change_for_frame_shift(mutation)
        elif vc == 'In_Frame_Del':
            adjusted_prot_change = self._get_prot_change_for_in_frame_del(mutation)
        elif vc == 'In_Frame_Ins':
            adjusted_prot_change = self._get_prot_change_for_in_frame_ins(mutation)
        elif vc == 'Nonstop_Mutation' or vc.startswith('Stop_Codon'):
            adjusted_prot_change = self._get_prot_change_for_stop_codon_variant(mutation)
        elif mutation['variant_type'] in ['DNP', 'TNP', 'ONP'] and '_' in mutation['protein_change']:
            #e.g. p.3589_3590QI>HV NEED TO MAKE TEST FOR THIS VARIANT!!!
            regx_res = PROT_ONP_REGEXP.match(mutation['protein_change'])
            aa_pos_1, aa_pos_2, ref_aa, alt_aa = [regx_res.group(i) for i in range(1, 5)]
            aa_pos_1, aa_pos_2 = int(aa_pos_1), int(aa_pos_2)
            adjusted_prot_change = 'p.%s%d_%s%ddelins%s' % (ref_aa[0], aa_pos_1, ref_aa[-1], aa_pos_2, alt_aa)
        else:
            regx_res = PROT_REGEXP.match(mutation['protein_change'])
            #try:
            ref_aa, aa_pos, alt_aa = [regx_res.group(i) for i in range(1, 4)]
            #except AttributeError:
            #    from IPython import embed; embed()
            adjusted_prot_change = 'p.%s%s%s' % (AA_NAME_MAPPING[ref_aa], aa_pos, AA_NAME_MAPPING[alt_aa])

        prot_id = self._get_ensembl_prot_id_from_tx_id(mutation['annotation_transcript'])
        adjusted_prot_change = '%s:%s' % (prot_id, adjusted_prot_change) if adjusted_prot_change else ''
        return adjusted_prot_change

    def _determine_overlap_type(self, mutation, tx_stop_codon_genomic_coords):
        """For determining the overlap betwen a deletion and the stop codon."""

        if mutation['transcript_strand'] == '+':
            if mutation['start'] >= tx_stop_codon_genomic_coords[0] + 1 and mutation['end'] <= tx_stop_codon_genomic_coords[1]:
                return 'del_within_stop_codon'
            elif mutation['start'] < tx_stop_codon_genomic_coords[0] + 1:
                return 'del_extends_into_cds'
            elif mutation['end'] > tx_stop_codon_genomic_coords[1]: #del extends into utr
                return 'del_extends_into_utr'
            else:
                raise Exception('Unexpected!!')
        else:
            if mutation['start'] >= tx_stop_codon_genomic_coords[0] + 1 and mutation['end'] <= tx_stop_codon_genomic_coords[1]:
                return 'del_within_stop_codon'
            elif mutation['end'] > tx_stop_codon_genomic_coords[1]:
                return 'del_extends_into_cds'
            elif mutation['start'] < tx_stop_codon_genomic_coords[0] + 1:
                return 'del_extends_into_utr'
            else:
                raise Exception('Unexpected!!')
 
    def _get_extension_amt(self, new_prot_seq):
        for i, extension_residue in enumerate(new_prot_seq):
            if extension_residue == '*':
                break
        return i

    def _get_prot_change_for_frame_shift(self, mutation):
        regx_res = PROT_FRAMESHIFT_REGEXP.match(mutation['protein_change'])
        ref_aa, aa_pos = [regx_res.group(i) for i in range(1, 3)]
        aa_pos = int(aa_pos)
        if ref_aa == '-':
            #will need to retrieve actual ref_aa if current ref_aa is "-"
            tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
            prot_seq = tx.get_protein_seq()
            ref_aa = prot_seq[aa_pos-1]
        ref_aa = ref_aa[0] #only need first residue
        return 'p.%s%dfs' % (AA_NAME_MAPPING[ref_aa], aa_pos)

    def _get_ensembl_prot_id_from_tx_id(self, tx_id):
        if '.' in tx_id:
            tx_id = tx_id[:tx_id.index('.')]
        return self.ensembl_id_mapping.get(tx_id)

    def _get_cdna_change_for_intron(self, mutation):
        #### HACK #####
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        nearest_exon = self.vcer._determine_closest_exon(tx, int(mutation['start']),
            int(mutation['end']))
        dist_to_exon = self.vcer._get_splice_site_coordinates(tx, int(mutation['start']),
            int(mutation['end']), nearest_exon)
        tx_exons = tx.get_exons()
    
        if dist_to_exon < 0 and mutation['transcript_strand'] == '+':
            cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(
                tx_exons[nearest_exon][0], tx_exons[nearest_exon][0], tx)
            if mutation['variant_classification'] == 'Intron':
                #do not do for splice sites
                dist_to_exon = dist_to_exon - 1 #why substract 1? why only for positive strand transcript
        elif dist_to_exon < 0 and mutation['transcript_strand'] == '-':
            cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(
                tx_exons[nearest_exon][1], tx_exons[nearest_exon][1], tx)
        elif dist_to_exon > 0 and mutation['transcript_strand'] == '+':
            cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(
                tx_exons[nearest_exon][1], tx_exons[nearest_exon][1], tx)
        elif dist_to_exon > 0 and mutation['transcript_strand'] == '-':
            cds_position_of_nearest_exon = TranscriptProviderUtils.convert_genomic_space_to_cds_space(
                tx_exons[nearest_exon][0], tx_exons[nearest_exon][0], tx)

        if cds_position_of_nearest_exon[0] == 0:
            #this means intron occurs before start codon and thus cds position needs to be adjusted to a negative number
            exon_coords = TranscriptProviderUtils.convert_genomic_space_to_exon_space(mutation['start'],
                tx.determine_cds_start() + 1, tx)
            cds_position_of_nearest_exon = exon_coords[0] - exon_coords[1]
        else:
            cds_position_of_nearest_exon = cds_position_of_nearest_exon[0]

        if dist_to_exon < 0:   
            cds_position_of_nearest_exon += 1 #why add 1?
        
        sign = '-' if dist_to_exon < 0 else '+'
        dist_to_exon = abs(dist_to_exon)

        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)

        adjusted_tx_change = 'c.%d%s%d%s>%s' % (cds_position_of_nearest_exon, sign, dist_to_exon,
            tx_ref_allele, tx_alt_allele)

        return '%s:%s' % (mutation['annotation_transcript'], adjusted_tx_change)

    def _get_cdna_change_for_5_utr(self, mutation):
        #### HACK #####
        ## Will not work correctly if variant is in exon that does not contain CDS start
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        dist_from_cds_start = abs(mutation['start'] - tx.determine_cds_start())
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        adjusted_tx_change = 'c.-%d%s>%s' % (dist_from_cds_start, tx_ref_allele, tx_alt_allele)
        return '%s:%s' % (mutation['annotation_transcript'], adjusted_tx_change)

    def _get_cdna_change_for_3_utr(self, mutation):
        #### HACK #####
        ## Will not work correctly if variant is in exon that does not contain CDS stop
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        dist_from_cds_stop = abs(mutation['start'] - tx.determine_cds_stop() + 2)
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        adjusted_tx_change = 'c.*%d%s>%s' % (dist_from_cds_stop, tx_ref_allele, tx_alt_allele)
        return '%s:%s' % (mutation['annotation_transcript'], adjusted_tx_change)

    def _get_tx_alleles(self, mutation):
        tx_ref_allele, tx_alt_allele = mutation['ref_allele'], mutation['alt_allele']
        if mutation['transcript_strand'] == '-':
            tx_ref_allele, tx_alt_allele = Seq.reverse_complement(tx_ref_allele), Seq.reverse_complement(tx_alt_allele)
        return tx_ref_allele, tx_alt_allele

    def _get_cdna_change_for_ONP(self, mutation):
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        tx_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(mutation['start'], mutation['end'], tx)
        tx_pos = tuple(t + 1 for t in tx_pos)

        tx_stop_codon_genomic_coords = tx.get_stop_codon()
        if (mutation['transcript_strand'] == '+' and mutation['end'] >= tx_stop_codon_genomic_coords[0] + 1) or \
            (mutation['transcript_strand'] == '-' and mutation['start'] <= tx_stop_codon_genomic_coords[1]):
        # bumping up against stop codon, need to determine tx coords this way
            tx_stop_start_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(tx_stop_codon_genomic_coords[0], tx_stop_codon_genomic_coords[1], tx)[0]
            tx_stop_start_pos += 1
            if mutation['transcript_strand'] == '+':
                genome_stop_start_pos = tx_stop_codon_genomic_coords[0] + 1
                tx_variant_pos_1_offset = mutation['start'] - genome_stop_start_pos
                tx_variant_pos_2_offset = mutation['end'] - genome_stop_start_pos
            else:
                genome_stop_start_pos = tx_stop_codon_genomic_coords[1]
                tx_variant_pos_1_offset = genome_stop_start_pos - mutation['end']
                tx_variant_pos_2_offset = genome_stop_start_pos - mutation['start']
    
            tx_pos_1 = tx_stop_start_pos + tx_variant_pos_1_offset
            tx_pos_2 = tx_stop_start_pos + tx_variant_pos_2_offset
            tx_pos = (tx_pos_1, tx_pos_2)

        adjusted_tx_change = 'c.%d_%ddelins%s' % (tx_pos[0], tx_pos[1], mutation['alt_allele'])
        return '%s:%s' % (mutation['annotation_transcript'], adjusted_tx_change)

    def _determine_if_genomic_duplication(self, alt_allele, ref_context, variant_type, mutation):
        allele_len = len(alt_allele)
        if variant_type == 'INS':
            half_len = len(ref_context) / 2
        elif variant_type == 'DEL':
            half_len = (len(ref_context) - len(alt_allele)) / 2
        
        downstream_seq = ref_context[half_len:half_len + allele_len].upper()
        upstream_seq = ref_context[half_len - allele_len:half_len].upper()

        return alt_allele == downstream_seq or alt_allele == upstream_seq

    def _genome_pos_adjust_if_duplication(self, alt_allele, ref_context, variant_type):
        if variant_type == 'INS':
            half_len = len(ref_context) / 2
        elif variant_type == 'DEL':
            half_len = (len(ref_context) - len(alt_allele)) / 2
        
        downstream_seq = ref_context[half_len:].upper()
        return self._adjust_pos_if_repeated_in_seq(alt_allele, downstream_seq)

    def _determine_if_prot_duplication(self, aa_pos_1, aa_pos_2, alt_allele, prot_seq, variant_type):
        allele_len = len(alt_allele)
        aa_pos_1, aa_pos_2 = int(aa_pos_1) - 1, int(aa_pos_2) - 1 #convert to 0-based indexing

        if variant_type == 'INS':
            downstream_seq = prot_seq[aa_pos_2:aa_pos_2 + allele_len]
            upstream_seq = prot_seq[aa_pos_1 - allele_len + 1:aa_pos_1 + 1]
        elif variant_type == 'DEL':
            downstream_seq = prot_seq[aa_pos_2 + 1:aa_pos_2 + 1 + allele_len]
            upstream_seq = prot_seq[aa_pos_1 - allele_len:aa_pos_1]

        return alt_allele == downstream_seq or alt_allele == upstream_seq

    def _determine_if_cdna_duplication(self, tx_pos_1, tx_pos_2, alt_allele, tx_seq, tx_strand, variant_type):
        allele_len = len(alt_allele)
        tx_pos_1, tx_pos_2 = int(tx_pos_1) - 1, int(tx_pos_2) - 1 #convert to 0-based indexing

        if variant_type == 'INS':
            if tx_strand == '+':
                wut = tx_seq[tx_pos_2:tx_pos_2 + allele_len] == alt_allele or \
                    tx_seq[tx_pos_1 - allele_len + 1:tx_pos_1 + 1] == alt_allele

                downstream_seq = tx_seq[tx_pos_2:tx_pos_2 + allele_len]
                upstream_seq = tx_seq[tx_pos_1 - allele_len + 1:tx_pos_1 + 1]

            elif tx_strand == '-':
                wut = tx_seq[tx_pos_2 + 1:tx_pos_2 + allele_len + 1] == alt_allele or \
                    tx_seq[tx_pos_2 - allele_len + 1:tx_pos_2 + 1] == alt_allele

                downstream_seq = tx_seq[tx_pos_2 + 1:tx_pos_2 + allele_len + 1]
                upstream_seq = tx_seq[tx_pos_2 - allele_len + 1:tx_pos_2 + 1]

        elif variant_type == 'DEL':
            if tx_strand == '+':
                downstream_seq = tx_seq[tx_pos_2 + 1:tx_pos_2 + allele_len + 1]
                upstream_seq = tx_seq[tx_pos_1 - allele_len:tx_pos_1]
            elif tx_strand == 'INS':
                downstream_seq = tx_seq[tx_pos_2 + 2:tx_pos_2 + allele_len + 2]
                upstream_seq = tx_seq[tx_pos_2 - allele_len:tx_pos_2]

        return alt_allele == downstream_seq or alt_allele == upstream_seq

    def _prot_pos_adjust_if_duplication(self, aa_pos_1, aa_pos_2, alt_allele, prot_seq):
        aa_pos_1, aa_pos_2 = int(aa_pos_1) - 1, int(aa_pos_2) - 1 #convert to 0-based indexing
        downstream_seq = prot_seq[aa_pos_2 + 1:]
        return self._adjust_pos_if_repeated_in_seq(alt_allele, downstream_seq)

    def _get_cdna_change_for_del(self, mutation):
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        len_alt_allele = len(tx_alt_allele)
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        tx_pos = TranscriptProviderUtils.convert_genomic_space_to_exon_space(mutation['start'], mutation['end'], tx)
        downstream_seq = tx.get_seq()[tx_pos[0]:]
        coding_pos_adjust = self._adjust_pos_if_repeated_in_seq(tx_alt_allele, downstream_seq)
        coding_pos_adjust = coding_pos_adjust - len_alt_allele

        coding_tx_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(mutation['start'], mutation['end'], tx)
        is_duplication = self._determine_if_cdna_duplication(tx_pos[0], tx_pos[1], tx_alt_allele, tx.get_seq(), mutation['transcript_strand'], 'INS')

        if is_duplication:
            coding_tx_pos = tuple(t + 1 for t in coding_tx_pos)
            start_pos = coding_tx_pos[0] + coding_pos_adjust
            start_pos = start_pos + 1 if mutation['transcript_strand'] == '-' else start_pos
            if len_alt_allele == 1:
                change_str = 'c.%ddel%s' % (start_pos, tx_ref_allele)
            else:
                change_str = 'c.%d_%ddel%s' % (start_pos, end_pos, tx_ref_allele)
        else:
            change_str = mutation['transcript_change']

        return '%s:%s' % (mutation['annotation_transcript'], change_str)

    def _get_cdna_change_for_ins(self, mutation):
        tx_ref_allele, tx_alt_allele = self._get_tx_alleles(mutation)
        len_alt_allele = len(tx_alt_allele)
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        tx_pos = TranscriptProviderUtils.convert_genomic_space_to_exon_space(mutation['start'], mutation['end'], tx)
        downstream_seq = tx.get_seq()[tx_pos[0]:]
        coding_pos_adjust = self._adjust_pos_if_repeated_in_seq(tx_alt_allele, downstream_seq)
        coding_pos_adjust = coding_pos_adjust - len_alt_allele

        coding_tx_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(mutation['start'], mutation['end'], tx)
        is_duplication = self._determine_if_cdna_duplication(tx_pos[0], tx_pos[1], tx_alt_allele, tx.get_seq(), mutation['transcript_strand'], 'INS')

        if is_duplication:
            coding_tx_pos = tuple(t + 1 for t in coding_tx_pos)
            start_pos = coding_tx_pos[0] + coding_pos_adjust
            start_pos = start_pos + 1 if mutation['transcript_strand'] == '-' else start_pos
            if len_alt_allele > 1:
                end_pos = start_pos + len_alt_allele - 1
                change = 'c.%d_%ddup%s' % (start_pos, end_pos, tx_alt_allele)
            else:
                change = 'c.%ddup%s' % (start_pos, tx_alt_allele)
        else:
            if mutation['transcript_strand'] == '-':
                coding_tx_pos = tuple(t + 1 for t in coding_tx_pos)

            tx_stop_codon_genomic_coords = tx.get_stop_codon()
            if (mutation['transcript_strand'] == '+' and mutation['end'] >= tx_stop_codon_genomic_coords[0] + 1) or \
                (mutation['transcript_strand'] == '-' and mutation['start'] <= tx_stop_codon_genomic_coords[1]):
            # bumping up against stop codon, need to determine tx coords this way
                tx_stop_start_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(tx_stop_codon_genomic_coords[0], tx_stop_codon_genomic_coords[1], tx)[0]
                tx_stop_start_pos += 1
                if mutation['transcript_strand'] == '+':
                    genome_stop_start_pos = tx_stop_codon_genomic_coords[0] + 1
                    tx_variant_pos_1_offset = mutation['start'] - genome_stop_start_pos
                    tx_variant_pos_2_offset = mutation['end'] - genome_stop_start_pos
                else:
                    genome_stop_start_pos = tx_stop_codon_genomic_coords[1]
                    tx_variant_pos_1_offset = genome_stop_start_pos - mutation['end']
                    tx_variant_pos_2_offset = genome_stop_start_pos - mutation['start']
    
                tx_variant_pos_1 = tx_stop_start_pos + tx_variant_pos_1_offset
                tx_variant_pos_2 = tx_stop_start_pos + tx_variant_pos_2_offset
                coding_tx_pos = (tx_variant_pos_1, tx_variant_pos_2)

            change = 'c.%d_%dins%s' % (coding_tx_pos[0], coding_tx_pos[1], tx_alt_allele)            

        return '%s:%s' % (mutation['annotation_transcript'], change)

    def _adjust_pos_if_repeated_in_seq(self, allele, downstream_seq):
        allele_len = len(allele)
        adjust_amt = 0
        for i in range(0, len(downstream_seq), allele_len):
            if allele == downstream_seq[i:i + allele_len]:
                adjust_amt += allele_len
            else:
                break

        return adjust_amt

    def _get_prot_change_for_in_frame_ins(self, mutation):
        is_insdel = False
        prot_change = mutation['protein_change']
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        prot_seq = tx.get_protein_seq()

        if 'ins' in prot_change:
            regx_res = PROT_INFRAME_INS_REGEXP_1.match(mutation['protein_change'])
            aa_pos_1, aa_pos_2, alt_allele = [regx_res.group(i) for i in range(1, 4)]
            aa_pos_1, aa_pos_2 = int(aa_pos_1), int(aa_pos_2)
            ref_allele_1, ref_allele_2 = prot_seq[aa_pos_1 - 1], prot_seq[aa_pos_2 - 1]
        elif '>' in prot_change:
            regx_res = PROT_INFRAME_INS_REGEXP_2.match(mutation['protein_change'])
            aa_pos_1, aa_pos_2, ref_allele, alt_allele = [regx_res.group(i) for i in range(1, 5)]
            aa_pos_1, aa_pos_2 = int(aa_pos_1), int(aa_pos_2)
            if aa_pos_1 == aa_pos_2 and alt_allele.startswith(ref_allele): #e.g. p.54_54A>AA
                alt_allele = alt_allele[1:]
            elif aa_pos_2 == aa_pos_1 + 1 and ref_allele[0] == alt_allele[0] and ref_allele[1] == alt_allele[-1]:
                #e.g. p.1326_1327KT>KET
                ref_allele_1 = ref_allele[0]
                ref_allele_2 = ref_allele[1]
                alt_allele = alt_allele[1:-1]
            else:
                #e.g. p.1327_1327T>RA
                is_insdel = True

        else:
            raise Exception('Unexpected protein_change str for in-frame insertion: %s' % prot_change)

        len_alt_allele = len(alt_allele)
        is_duplication = self._determine_if_prot_duplication(aa_pos_1, aa_pos_2, alt_allele, prot_seq, mutation['variant_type'])
        prot_pos_adjust_amnt =  self._prot_pos_adjust_if_duplication(aa_pos_1, aa_pos_2, alt_allele, prot_seq)

        alt_allele = ''.join(AA_NAME_MAPPING[aa] for aa in alt_allele)

        if is_duplication:
            #aa_pos_1, aa_pos_2 = int(aa_pos_1) + prot_pos_adjust_amnt, int(aa_pos_2) + prot_pos_adjust_amnt
            if alt_allele == '*':
                return '' # duplication of stop does not change protein
            elif aa_pos_1 == aa_pos_2:
                aa_pos_1 = aa_pos_1 + prot_pos_adjust_amnt
                adjusted_prot_change = 'p.%s%ddup' % (alt_allele, aa_pos_1)
            else: #means that prot positions describe residues that alt allele gets inserted between
                if prot_pos_adjust_amnt == 0: #duplication is already at 3'-most position
                    aa_pos_2 = aa_pos_1
                    aa_pos_1 = aa_pos_1 - len_alt_allele + 1
                else:
                    raise Exception('Need to implement this!!')

                adjusted_prot_change = 'p.%s%d_%s%ddup' % (alt_allele[:3], aa_pos_1, alt_allele[-3:], aa_pos_2)
        elif is_insdel:
            adjusted_prot_change = 'p.%s%ddelins%s' % (AA_NAME_MAPPING[ref_allele], aa_pos_1, alt_allele)
        else:
            ref_allele_1, ref_allele_2 = AA_NAME_MAPPING[ref_allele_1], AA_NAME_MAPPING[ref_allele_2]
            if alt_allele == '*':
                adjusted_prot_change = 'p.%s%d*' % (ref_allele_2, aa_pos_2)
            else:
                adjusted_prot_change = 'p.%s%d_%s%dins%s' % (ref_allele_1, aa_pos_1, ref_allele_2, aa_pos_2, alt_allele)

        return adjusted_prot_change

    def _get_prot_change_for_in_frame_del(self, mutation):
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        prot_seq = tx.get_protein_seq()

        if 'del' in mutation['protein_change']:
            regx_res = PROT_INFRAME_DEL_REGEXP.match(mutation['protein_change'])
            ref_aa, aa_pos = [regx_res.group(i) for i in range(1, 3)]
            aa_pos_1 = int(aa_pos)

            if len(ref_aa) > 1:
                ref_aa_1, ref_aa_2 = ref_aa[0], ref_aa[-1]
                
                aa_pos_2 = aa_pos_1 + len(ref_aa) - 1
                
            else:
                aa_pos_2 = aa_pos_1
            
            prot_pos_adjust_amnt = self._prot_pos_adjust_if_duplication(aa_pos_1, aa_pos_2, ref_aa, prot_seq)
            aa_pos_1, aa_pos_2 = aa_pos_1 + prot_pos_adjust_amnt, aa_pos_2 + prot_pos_adjust_amnt
            if len(ref_aa) > 1:
                return 'p.%s%d_%s%ddel' % (AA_NAME_MAPPING[ref_aa_1], aa_pos_1, AA_NAME_MAPPING[ref_aa_2], aa_pos_2)
            else:
                return 'p.%s%ddel' % (AA_NAME_MAPPING[ref_aa], aa_pos_1)

        elif '>' in mutation['protein_change']:
            regx_res = PROT_INFRAME_DEL_REGEXP_2.match(mutation['protein_change'])
            aa_pos_1, aa_pos_2, ref_aa, alt_aa = [regx_res.group(i) for i in range(1, 5)]
            aa_pos_1, aa_pos_2 = int(aa_pos_1), int(aa_pos_2)
    
            if ref_aa.endswith(alt_aa):
                #e.g. p.100_101EG>G
                new_alt = ref_aa[:-len(alt_aa)]
                aa_pos_2 = aa_pos_1 + len(new_alt) - 1
            elif ref_aa.startswith(alt_aa):
                #e.g. p.100_101EG>E
                new_alt = ref_aa[len(alt_aa):]
                aa_pos_1 = aa_pos_2 - len(new_alt) + 1
            else:
                #e.g. p.100_101EG>K
                ref_aa_1, ref_aa_2 = ref_aa[0], ref_aa[-1]
                return 'p.%s%d_%s%ddelins%s' % (AA_NAME_MAPPING[ref_aa_1], aa_pos_1,
                    AA_NAME_MAPPING[ref_aa_2], aa_pos_2, AA_NAME_MAPPING[alt_aa])
        
            prot_pos_adjust_amnt = self._prot_pos_adjust_if_duplication(aa_pos_1, aa_pos_2, new_alt, prot_seq)
            aa_pos_1, aa_pos_2 = aa_pos_1 + prot_pos_adjust_amnt, aa_pos_2 + prot_pos_adjust_amnt

            if len(new_alt) == 1:
                return 'p.%s%sdel' % (AA_NAME_MAPPING[new_alt], aa_pos_1)
            else:
                new_alt_1, new_alt_2 = new_alt[0], new_alt[-1]
                return 'p.%s%d_%s%ddel' % (AA_NAME_MAPPING[new_alt_1], aa_pos_1, AA_NAME_MAPPING[new_alt_2], aa_pos_2)

    def _get_prot_change_for_stop_codon_variant(self, mutation):
        tx = self.gencode_ds.transcript_db[mutation['annotation_transcript']]
        tx_seq = tx.get_seq()
        tx_stop_codon_genomic_coords = tx.get_stop_codon()
        tx_stop_codon_tx_coords = TranscriptProviderUtils.convert_genomic_space_to_exon_space(tx_stop_codon_genomic_coords[0],
            tx_stop_codon_genomic_coords[1], tx)
        prot_seq = tx.get_protein_seq()

        stop_codon_seq = tx_seq[tx_stop_codon_tx_coords[0]:tx_stop_codon_tx_coords[1]]
        utr_seq = tx_seq[tx_stop_codon_tx_coords[1]:]

        variant_tx_coords = TranscriptProviderUtils.convert_genomic_space_to_exon_space(mutation['start'], mutation['end'], tx)
        variant_tx_pos_1 = variant_tx_coords[0]
        variant_tx_pos_1 = variant_tx_pos_1 - 1 if mutation['transcript_strand'] == '+' else variant_tx_pos_1
        variant_codon_pos = variant_tx_pos_1 - tx_stop_codon_tx_coords[0]

        prot_stop_pos = TranscriptProviderUtils.convert_genomic_space_to_cds_space(tx_stop_codon_genomic_coords[0], tx_stop_codon_genomic_coords[0], tx)[0]
        prot_stop_pos = (prot_stop_pos / 3) + 1

        if mutation['variant_type'] == 'SNP':
            new_stop_codon_seq = list(stop_codon_seq)
            new_stop_codon_seq[variant_codon_pos] = mutation['alt_allele']
            new_stop_codon_seq = ''.join(new_stop_codon_seq)
            new_tx_seq = new_stop_codon_seq + utr_seq
            new_prot_seq = Seq.translate(new_tx_seq)
        elif mutation['variant_type'] == 'DEL':
            overlap_type = self._determine_overlap_type(mutation, tx_stop_codon_genomic_coords)

            if overlap_type == 'del_within_stop_codon':
                #del_in_just_stop_codon
                #if mutation['start'] == 29416090: from IPython import embed; embed()

                new_stop_codon_seq = list(stop_codon_seq)
                n_bases_to_del = len(mutation['ref_allele'])
                codon_positions_to_del = list()
                for i in range(n_bases_to_del):
                    codon_positions_to_del.append(variant_codon_pos + i)
                
                for i in codon_positions_to_del[::-1]:
                    del new_stop_codon_seq[i]

                new_stop_codon_seq = ''.join(new_stop_codon_seq)
                new_tx_seq = new_stop_codon_seq + utr_seq
                new_stop_codon_seq = new_tx_seq[:3]
                new_prot_seq = Seq.translate(new_tx_seq)

            elif overlap_type == 'del_extends_into_cds':
                if mutation['transcript_strand'] == '+':
                    n_bases_into_cds = (tx_stop_codon_genomic_coords[0] + 1) - mutation['start']
                else:
                    n_bases_into_cds = mutation['end'] - tx_stop_codon_genomic_coords[1]

                n_codons_into_cds = 0
                for i in range(1, n_bases_into_cds + 1):
                    if i % 3 == 0:
                        n_codons_into_cds += 1

                n_codons_into_cds = 1 if n_codons_into_cds == 0 else n_codons_into_cds
                if i > 3 and i % 3 != 0:
                    n_codons_into_cds += 1 #to get num of codons correct
                
                bases_into_cds_deleted = n_bases_into_cds
                bases_into_stop_deleted = len(mutation['ref_allele']) - n_bases_into_cds

                new_ref_codon_seq = tx_seq[tx_stop_codon_tx_coords[0] - (3 * n_codons_into_cds): tx_stop_codon_tx_coords[0]]
                new_alt_codon_seq = new_ref_codon_seq[:-bases_into_cds_deleted]

                if bases_into_stop_deleted > 3:
                    #del extends into UTR
                    bases_into_utr = bases_into_stop_deleted - 3
                    bases_into_stop_deleted = 3
                else:
                    bases_into_utr = 0

                new_cds_seq = new_alt_codon_seq
                new_alt_codon_seq = new_cds_seq + stop_codon_seq[bases_into_stop_deleted:]
                new_ref_codon_seq = new_ref_codon_seq + stop_codon_seq

                ref_aa = Seq.translate(new_ref_codon_seq)
                alt_aa = Seq.translate(new_alt_codon_seq)

                aa_pos_2 = prot_stop_pos
                aa_pos_1 = aa_pos_2 - n_codons_into_cds

                if '*' in alt_aa:
                    #no extension
                    alt_aa = ''.join(AA_NAME_MAPPING[aa] for aa in alt_aa)
                    return 'p.%s%d_%s%ddelins%s' % (aa_pos_1, AA_NAME_MAPPING[ref_aa[1]], aa_pos_2,
                        AA_NAME_MAPPING[ref_aa[-1]], alt_aa)
                else:
                    new_tx_seq = new_cds_seq + stop_codon_seq[bases_into_stop_deleted:] + utr_seq[bases_into_utr:]
                    new_prot_seq = Seq.translate(new_tx_seq)
                    #new_tx_seq = stop_codon_seq[bases_into_stop_deleted:] + utr_seq[bases_into_utr:]
                    #new_prot_seq = Seq.translate(new_tx_seq[n_bases_into_cds:])

                    if new_prot_seq[0] == '*': #no change to stop codon seq
                        if ref_aa == alt_aa + '*':
                            return '' #no change on protein sequence
                        else:
                            alt_aa = alt_aa + '*'
                            alt_aa = ''.join(AA_NAME_MAPPING[aa] for aa in alt_aa)
                            'p.%s%d_%s%ddelins%s' % (AA_NAME_MAPPING[ref_aa[1]], aa_pos_1,
                                AA_NAME_MAPPING[ref_aa[-1]], aa_pos_2, alt_aa)
                    
                    elif ref_aa == alt_aa + '*': #means that protein seq upstream of stop is not being affected
                        if '*' not in new_prot_seq: # is extension with no new stop codon:
                            return 'p.*%d%sext*?' % (prot_stop_pos, AA_NAME_MAPPING[new_prot_seq[len(alt_aa)]])
                        else:
                            extension_amt = self._get_extension_amt(new_prot_seq)
                            return 'p.*%d%sext*%d' % (prot_stop_pos, AA_NAME_MAPPING[new_prot_seq[len(alt_aa)]], extension_amt)
                    else:
                        alt_aa = ''.join(AA_NAME_MAPPING[aa] for aa in alt_aa)
                        if alt_aa == '':
                            alt_aa = AA_NAME_MAPPING[new_prot_seq[0]]
                        if '*' not in new_prot_seq: # is extension with no new stop codon:
                            return 'p.%s%d_%s%ddelins%sext*?' % (AA_NAME_MAPPING[ref_aa[0]], aa_pos_1,
                                AA_NAME_MAPPING[ref_aa[-1]], aa_pos_2, alt_aa)
                        else:
                            extension_amt = self._get_extension_amt(new_prot_seq)
                            return 'p.%s%d_%s%ddelins%sext*%d' % (AA_NAME_MAPPING[ref_aa[0]], aa_pos_1,
                                AA_NAME_MAPPING[ref_aa[-1]], aa_pos_2, alt_aa, extension_amt)


            elif overlap_type == 'del_extends_into_utr':
                del_ref_allele = self._get_tx_alleles(mutation)[0]

                for i in range(3):
                    if del_ref_allele.startswith(stop_codon_seq[i:]):
                        break
                n_bases_into_stop_deleted = 3 - i
                n_bases_into_utr_deleted = len(del_ref_allele) - n_bases_into_stop_deleted
                new_tx_seq = stop_codon_seq[:-n_bases_into_stop_deleted] + utr_seq[n_bases_into_utr_deleted:]
                new_prot_seq = Seq.translate(new_tx_seq)
                if new_prot_seq[0] == '*':
                    return '' #no change on protein sequence
                elif '*' not in new_prot_seq: # is extension with no new stop codon
                    return 'p.*%d%sext*?' % (prot_stop_pos, AA_NAME_MAPPING[new_prot_seq[0]])
                else:
                    extension_amt = self._get_extension_amt(new_prot_seq)
                    return 'p.*%d%sext*%d' % (prot_stop_pos, AA_NAME_MAPPING[new_prot_seq[0]], extension_amt)


        elif mutation['variant_type'] == 'INS':
            alt_allele = self._get_tx_alleles(mutation)[1]
            if variant_codon_pos == 2:
                return '' # we don't care about insertions after stop codon, no change in protein
            else:
                new_tx_seq = stop_codon_seq[:variant_codon_pos-2] + alt_allele + stop_codon_seq[variant_codon_pos+1:] + utr_seq
                new_prot_seq = Seq.translate(new_tx_seq)
                alt_aa = new_prot_seq[0]

                if alt_aa == '*':
                    return '' #no change in protein
                elif '*' not in new_prot_seq: # is extension with no new stop codon
                    return 'p.*%d%sext*?' % (prot_stop_pos, AA_NAME_MAPPING[alt_aa])
                else:
                    extension_amt = self._get_extension_amt(new_prot_seq)
                    return 'p.*%d%sext*%d' % (prot_stop_pos, AA_NAME_MAPPING[alt_aa], extension_amt)

        else: # ONPs
            regx_res = PROT_ONP_REGEXP.match(mutation['protein_change'])
            aa_pos_1, aa_pos_2, ref_aa, alt_aa = [regx_res.group(i) for i in range(1, 5)]
            aa_pos_1, aa_pos_2 = int(aa_pos_1), int(aa_pos_2)

            new_prot_seq = Seq.translate(utr_seq)
            alt_aa = [AA_NAME_MAPPING[aa] for aa in alt_aa]
            alt_aa = ''.join(alt_aa)
            if alt_aa[0] == '*':
                return 'p.%s%d_%s%ddelins%s' % (AA_NAME_MAPPING[ref_aa[0]], aa_pos_1,
                    AA_NAME_MAPPING[ref_aa[-1]], aa_pos_2, '*')
            elif '*' in alt_aa:
                return 'p.%s%d_%s%ddelins%s' % (AA_NAME_MAPPING[ref_aa[0]], aa_pos_1,
                    AA_NAME_MAPPING[ref_aa[-1]], aa_pos_2, alt_aa)
            elif '*' not in new_prot_seq:
                return 'p.%s%d_%s%ddelins%sext*?' % (AA_NAME_MAPPING[ref_aa[0]], aa_pos_1,
                    AA_NAME_MAPPING[ref_aa[-1]], aa_pos_2, alt_aa)
            else:
                extension_amt = self._get_extension_amt(new_prot_seq)
                extension_amt += 1 #to account for aa in stop codon
                return 'p.%s%d_%s%ddelins%sext*%d' % (AA_NAME_MAPPING[ref_aa[0]], aa_pos_1,
                    AA_NAME_MAPPING[ref_aa[-1]], aa_pos_2, alt_aa, extension_amt)


        alt_aa = Seq.translate(new_stop_codon_seq)
        if new_prot_seq[0] == '*': #no change to stop codon seq
            return ''
        elif '*' not in new_prot_seq: # is extension with no new stop codon
            return 'p.*%d%sext*?' % (prot_stop_pos, AA_NAME_MAPPING[alt_aa])
        else: # is extension with downstream stop codon
            extension_amt = self._get_extension_amt(new_prot_seq)
            return 'p.*%d%sext*%d' % (prot_stop_pos, AA_NAME_MAPPING[alt_aa], extension_amt)
