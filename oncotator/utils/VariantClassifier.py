import Bio
from Bio import Seq
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils


class VariantClassifier(object):

    def __init__(self):
        pass

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
            # TODO: Handle a frameshift here.

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

            variant_classification = infer_variant_classification(m['variant_type'],
                reference_aa, observed_aa, reference_allele, observed_allele)

            # If silent mutation w/in 2 bp of a splice junction, then change to splice site
            if variant_classification.lower() == "silent":
                determineIfSpliceSiteThrown(end, exon_affected,gaf,start,t, dist=2)

            if variant_type != 'SNP':
                reference_aa, observed_aa, protein_position_start, protein_position_end = \
                    adjust_protein_position_and_alleles(protein_seq, protein_position_start,
                        protein_position_end, reference_aa, observed_aa)
        else:
            annotate_mutations_not_fully_within_cds(transcript_seq, observed_allele, cds_overlap_type, transcript_position_start,
                transcript_position_end, m, t, exon_affected)

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