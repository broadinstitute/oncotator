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
import shove
from oncotator.Annotation import Annotation
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources.Datasource import Datasource
from oncotator.datasources.SegmentDatasource import SegmentDatasource
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.index.gaf import region2bins
from oncotator.utils.Hasher import Hasher
from oncotator.utils.HgvsChangeTransformer import HgvsChangeTransformer
from oncotator.utils.VariantClassification import VariantClassification
from oncotator.utils.VariantClassifier import VariantClassifier
from oncotator.utils.txfilter.TranscriptFilterFactory import TranscriptFilterFactory


class EnsemblTranscriptDatasource(TranscriptProvider, Datasource, SegmentDatasource):
    """ Similar to a GAF datasource, but uses ensembl transcripts.
        Also, supports gencode

        Though all transcripts for GENCODE can be loaded, it is currently set to ignore any transcripts that are not
            "basic"


        Any transcripts in the custom canonical transcript list will be selected before the tx-mode parameter.
    """
    """This is the list of annotations that get populated by this datasource"""
    POPULATED_ANNOTATION_NAMES = {'transcript_exon', 'variant_type', 'variant_classification', 'other_transcripts',
                                  'gene', 'gene_id', 'annotation_transcript', 'genome_change', 'strand',
                                  'transcript_id', 'secondary_variant_classification', 'protein_change', 'codon_change',
                                  'transcript_change', 'transcript_strand', 'gene', 'gene_type',
                                  'gencode_transcript_tags', 'gencode_transcript_status', 'havana_transcript',
                                  'ccds_id', 'gencode_transcript_type', 'transcript_position',
                                  'gencode_transcript_name'}

    def __init__(self,  src_file, title='ENSEMBL', version='', tx_mode=TranscriptProvider.TX_MODE_CANONICAL, protocol="file", is_thread_safe=False, tx_filter="dummy", custom_canonical_txs=None):
        super(EnsemblTranscriptDatasource, self).__init__(src_file=src_file, title=title, version=version)

        ensembl_index_fname = src_file + ".transcript.idx"
        ensembl_gene_to_transcript_index_fname = src_file + ".transcript_by_gene.idx"
        ensembl_genomic_position_bins_to_transcript_index_fname = src_file + ".transcript_by_gp_bin.idx"

        # Seconds before a cache entry should be cleared out
        timeout = 1000
        max_entries = 25000
        cache_protocol = "memory"
        if not is_thread_safe:
            logging.getLogger(__name__).info("%s %s is being set up in faster, NOT thread-safe mode (for annotation).  " % (title, version))
            cache_protocol = "simple"

        # Contains a key of transcript id and value of a Transcript class, with sequence data where possible.
        # By specifying "memory" for the cache, this is thread safe.  Otherwise, use "simple"
        self.transcript_db = shove.Shove(protocol + '://%s' % ensembl_index_fname, cache_protocol + "://", timeout=timeout, max_entries=max_entries)
        self.gene_db = shove.Shove(protocol + '://%s' % ensembl_gene_to_transcript_index_fname, cache_protocol + "://", timeout=timeout, max_entries=max_entries)
        self.gp_bin_db = shove.Shove(protocol + '://%s' % ensembl_genomic_position_bins_to_transcript_index_fname, cache_protocol + "://", timeout=timeout, max_entries=max_entries)

        tmp = self.gp_bin_db.keys()

        logging.getLogger(__name__).info("%s %s is being set up with default tx-mode: %s.  " % (title, version, tx_mode))
        self.set_tx_mode(tx_mode)

        logging.getLogger(__name__).info("%s %s is being set up with %s filtering.  " % (title, version, tx_filter))
        self._tx_filter = TranscriptFilterFactory.create_instance(tx_filter)

        self._hgvs_xformer = HgvsChangeTransformer()

        # Store a list of the custom canonical transcripts
        self._custom_canonical_txs = custom_canonical_txs or []

        # IMPORTANT: Any new attributes that can change the results of annotations and, therefore, should invalidate the
        #  cache, should be added to the list in get_hashcode.  There should be a way to do this dynamically, but that
        # has not been implemented yet.

    def set_tx_mode(self, tx_mode):
        if tx_mode == TranscriptProvider.TX_MODE_CANONICAL:
            logging.getLogger(__name__).warn("Attempting to set transcript mode of CANONICAL for ensembl.  This operation is only supported for GENCODE.  Otherwise, will be the same as EFFECT.")
        self.tx_mode = tx_mode

    def _create_basic_annotation(self, value):
        return Annotation(value=value, datasourceName=self.title)

    def _create_blank_set_of_annotations(self):
        final_annotation_dict = dict()
        for k in EnsemblTranscriptDatasource.POPULATED_ANNOTATION_NAMES:
            final_annotation_dict[k] = self._create_basic_annotation('')

        return final_annotation_dict

    def _retrieve_gencode_tag_value(self, tx, attribute_name):
        """
        If transcript is not gencode, no error is thrown.  Just a blank value.
        Note that gencode have other attributes, but not plain ol' ENSEMBL
        :param tx:
        :param attribute_name:
        :return: "" if other attributes are not present.  "" if specified tag is not present.  Otherwise, tag value.
        """
        attribute_dict = tx.get_other_attributes()
        if attribute_dict is None:
            return ""
        return str(attribute_dict.get(attribute_name, ""))

    def get_transcripts_by_pos(self, chr, start, end):
        """
        Returns filtered list of transcripts that overlap the given genomic position.

        :rtype : list
        :param str chr:
        :param str|int start:
        :param str|int end:
        """
        txs_unfiltered = self.get_overlapping_transcripts(chr, start, end)
        txs = self._filter_transcripts(txs_unfiltered)
        return txs

    def _create_hgvs_dict(self, chosen_tx, mutation):
        hgvs_dict = dict()
        if self._hgvs_xformer is not None:
            hgvs_dict = self._hgvs_xformer.hgvs_annotate_mutation_given_tx(mutation, chosen_tx)
        return hgvs_dict

    def _create_hgvs_annotation_dict(self, mutation, chosen_tx):
        hgvs_dict_keys = HgvsChangeTransformer.HEADERS
        hgvs_dict = self._create_hgvs_dict(chosen_tx, mutation)
        hgvs_dict_annotations = dict()
        for k in hgvs_dict_keys:
            hgvs_dict_annotations[k] = self._create_basic_annotation(hgvs_dict.get(k, ""))
        return hgvs_dict_annotations

    def annotate_mutation(self, mutation):
        chr = mutation.chr
        start = int(mutation.start)
        end = int(mutation.end)
        txs = self.get_transcripts_by_pos(chr, start, end)
        final_annotation_dict = self._create_blank_set_of_annotations()
        final_annotation_dict['variant_type'] = Annotation(value=TranscriptProviderUtils.infer_variant_type(mutation.ref_allele, mutation.alt_allele), datasourceName=self.title)
        chosen_tx = None

        # We have hit IGR if no transcripts come back.  Most annotations can just use the blank set.
        if len(txs) == 0:
            final_annotation_dict['variant_classification'] = self._create_basic_annotation(VariantClassification.IGR)
            nearest_genes = self._get_nearest_genes(chr, int(start), int(end))
            final_annotation_dict['other_transcripts'] = self._create_basic_annotation(value='%s (%s upstream) : %s (%s downstream)' % (nearest_genes[0][0], nearest_genes[0][1], nearest_genes[1][0], nearest_genes[1][1]))
            final_annotation_dict['gene'] = self._create_basic_annotation('Unknown')
            final_annotation_dict['gene_id'] = self._create_basic_annotation('0')
            final_annotation_dict['genome_change'] = self._create_basic_annotation(TranscriptProviderUtils.determine_genome_change(mutation.chr, mutation.start, mutation.end, mutation.ref_allele, mutation.alt_allele, final_annotation_dict['variant_type'].value))
        else:
            # Choose the best effect transcript
            chosen_tx = self._choose_transcript(txs, self.get_tx_mode(), final_annotation_dict['variant_type'].value, mutation.ref_allele, mutation.alt_allele, start, end)
            vcer = VariantClassifier()

            final_annotation_dict['annotation_transcript'] = self._create_basic_annotation(chosen_tx.get_transcript_id())
            final_annotation_dict['genome_change'] = self._create_basic_annotation(TranscriptProviderUtils.determine_genome_change(mutation.chr, mutation.start, mutation.end, mutation.ref_allele, mutation.alt_allele, final_annotation_dict['variant_type'].value))
            final_annotation_dict['strand'] = self._create_basic_annotation(chosen_tx.get_strand())

            final_annotation_dict['transcript_position'] = self._create_basic_annotation(TranscriptProviderUtils.render_transcript_position(int(start), int(end), chosen_tx))

            final_annotation_dict['transcript_id'] = self._create_basic_annotation(chosen_tx.get_transcript_id())

            variant_classfication = vcer.variant_classify(tx=chosen_tx, variant_type=final_annotation_dict['variant_type'].value,
                                             ref_allele=mutation.ref_allele, alt_allele=mutation.alt_allele, start=mutation.start, end=mutation.end)
            final_annotation_dict['transcript_exon'] = self._create_basic_annotation(str(variant_classfication.get_exon_i()+1))
            final_annotation_dict['variant_classification'] = self._create_basic_annotation(variant_classfication.get_vc())
            final_annotation_dict['secondary_variant_classification'] = self._create_basic_annotation(variant_classfication.get_secondary_vc())
            final_annotation_dict['protein_change'] = self._create_basic_annotation(vcer.generate_protein_change_from_vc(variant_classfication))
            final_annotation_dict['codon_change'] = self._create_basic_annotation(vcer.generate_codon_change_from_vc(chosen_tx, start, end, variant_classfication))
            final_annotation_dict['transcript_change'] = self._create_basic_annotation(vcer.generate_transcript_change_from_tx(chosen_tx, final_annotation_dict['variant_type'].value, variant_classfication, start, end, mutation.ref_allele, mutation.alt_allele))

            final_annotation_dict['transcript_strand'] = self._create_basic_annotation(chosen_tx.get_strand())
            final_annotation_dict['gene'] = self._create_basic_annotation(chosen_tx.get_gene())
            final_annotation_dict['gene_type'] = self._create_basic_annotation(chosen_tx.get_gene_type())
            final_annotation_dict['gencode_transcript_tags'] = self._create_basic_annotation(self._retrieve_gencode_tag_value(chosen_tx, 'tag'))
            final_annotation_dict['gencode_transcript_status'] = self._create_basic_annotation(self._retrieve_gencode_tag_value(chosen_tx, 'transcript_status'))
            final_annotation_dict['havana_transcript'] = self._create_basic_annotation(self._retrieve_gencode_tag_value(chosen_tx, 'havana_transcript'))
            final_annotation_dict['ccds_id'] = self._create_basic_annotation(self._retrieve_gencode_tag_value(chosen_tx, 'ccdsid'))
            final_annotation_dict['gencode_transcript_type'] = self._create_basic_annotation(self._retrieve_gencode_tag_value(chosen_tx, 'transcript_type'))
            final_annotation_dict['gencode_transcript_name'] = self._create_basic_annotation(self._retrieve_gencode_tag_value(chosen_tx, 'transcript_name'))

            other_transcript_value = self._render_other_transcripts(txs, [txs.index(chosen_tx)], final_annotation_dict['variant_type'].value, mutation.ref_allele, mutation.alt_allele, mutation.start, mutation.end)
            final_annotation_dict['other_transcripts'] = self._create_basic_annotation(other_transcript_value)
            # final_annotation_dict['gene_id'].value

        mutation.addAnnotations(final_annotation_dict)

        # Add the HGVS annotations ... setting to "" if not available.
        hgvs_dict_annotations = self._create_hgvs_annotation_dict(mutation, chosen_tx)
        mutation.addAnnotations(hgvs_dict_annotations)

        return mutation

    def _filter_transcripts(self, txs):
        return self._tx_filter.filter(txs)

    def _choose_transcript(self, txs, tx_mode, variant_type, ref_allele, alt_allele, start, end):
        """Given a list of transcripts and a transcript mode (e.g. CANONICAL), choose the transcript to use.
        :param list txs: a list of transcripts that presumably overlap the variant
        :param tx_mode:
        :param str variant_type:
        :param str ref_allele:
        :param str alt_allele:
        :param start:
        :param end:
        :return Transcript : chosen transcript given tx-mode
        """
        if len(txs) == 1:
            return txs[0]
        if tx_mode == TranscriptProvider.TX_MODE_CANONICAL:
            return self._choose_canonical_transcript(txs, variant_type, ref_allele, alt_allele, start, end)
        return self._choose_best_effect_transcript(txs, variant_type, ref_allele, alt_allele, start, end)

    @staticmethod
    def _get_best_scores(txs, scoring_function, comparator):
        scores = {tx: scoring_function(tx) for tx in txs}
        best = comparator(scores.itervalues())
        return [k for (k, v) in scores.iteritems() if v == best]

    @staticmethod
    def _select_best_with_multiple_criteria(txs, scoring_functions):
        """Sort using multiple scoring functions
        :param txs: transcripts to sort
        :param scoring_functions: a tuple of the form ( tx -> B, [B] -> B)

        an example is (lambda x: x.seq_len(), min)
        which will apply seq_len to each transcript, and then use the minimum resulting length
        """
        best_txs = txs
        for (f, cmp) in scoring_functions:
            if len(best_txs) == 1:
                return best_txs
            else:
                best_txs = EnsemblTranscriptDatasource._get_best_scores(best_txs, f, cmp)
        return best_txs

    @staticmethod
    def _calculate_canonical_score(tx):
        """

        For score, higher is better.

        "Protein coding" is worth four points.
        4 - (gencode transcript level) are added.

        For example, a protein coding level 2 (automatic curation) transcript will score 6.
        For example, a protein coding level 1 (automatic and manual curation) transcript will score 7.
        Level 1 is validated
        Level 2 is manual annotation
        Level 3 is automated annotation.
        :param tx: Transcript
        :return int: score as described above
        """
        # higher ranks are more important.
        lvl_rank = 0
        lvl = tx.get_other_attributes().get('level', [None])[0]
        if lvl is None:
            lvl_score = 0
        else:
            lvl_score = 4 - int(lvl)

        type_rank = 2
        type_score = 0
        if tx.get_gene_type() == "protein_coding":
            type_score = 1

        return (lvl_score << lvl_rank) + (type_score << type_rank)

    @staticmethod
    def _get_appris_rank(tx):
        """Get the appris ranking from the transcript's tag field"""
        appris_ranks = TranscriptProviderUtils.APPRIS_RANKING_DICT
        tags =  set(tx.get_other_attributes().get('tag', "").split("|"))
        for tag in TranscriptProviderUtils.APPRIS_TAGS:
            if tag in tags:
                if "CCDS" in tags:
                    if tag == "appris_candidate":
                        return appris_ranks["appris_candidate_ccds"]
                    elif tag in {'appris_candidate_longest_seq', 'appris_candidate_longest'}:
                        return appris_ranks["appris_candidate_longest_ccds"]
                return appris_ranks[tag]
        else:
            return TranscriptProviderUtils.NO_APPRIS_VALUE

    @staticmethod
    def _calculate_effect_score(tx, start, end, alt_allele, ref_allele, variant_type):
        """Compute the effect score"""
        effect_dict = TranscriptProviderUtils.retrieve_effect_dict()
        vcer = VariantClassifier()
        if (ref_allele == "" or ref_allele == "-") and (alt_allele == "" or alt_allele == "-"):
            vc = VariantClassification.SILENT
        else:
            vc = vcer.variant_classify(tx, ref_allele, alt_allele, start, end, variant_type).get_vc()
        effect_score = effect_dict.get(vc, 25)
        return effect_score

    @staticmethod
    def _is_tx_in_tx_list(tx, tx_id_list):
        """
        Check if a transcript ID (with or without version number) is in the given list of transcripts.

        NOTE: The tx_id_list is assumed to not include version numbers.
        :param tx:
        :param tx_id_list:
        :return bool:
        """
        tx_id = tx.get_transcript_id()
        if tx_id.rsplit('.', 1)[0] in tx_id_list:
            return True
        return False

    def _choose_best_effect_transcript(self, txs, variant_type, ref_allele, alt_allele, start, end):
        """Choose the transcript with the most detrimental effect.
         The rankings are in TranscriptProviderUtils.
         Ties are broken by which transcript has the longer coding length.
         Additional ties are broken with appris rank

        0.  membership in custom canonical transcript list
        1.  most detrimental effect
        2.  curation level
        3.  appris rank
        4.  longest protein change
        5.  lexicographical sort on transcript ID

        :param list txs: list of Transcript
        :param str variant_type:
        :param str ref_allele:
        :param str alt_allele:
        :param str start:
        :param str end:
        :return Transcript:
         """
        if len(txs) == 0:
            return None
        best_effect_txs = EnsemblTranscriptDatasource._select_best_with_multiple_criteria(txs,
                        [(lambda x: self._is_tx_in_tx_list(x, self._custom_canonical_txs), max),
                        (lambda x: self._calculate_effect_score(x, start, end, alt_allele, ref_allele, variant_type), min),
                        (self._calculate_canonical_score, max),
                        (self._get_appris_rank, min),
                        (lambda x: len(x.get_seq()), max),
                        (lambda x: x.get_transcript_id(), min)]
                    )
        return best_effect_txs[0]

    def _choose_canonical_transcript(self, txs, variant_type, ref_allele, alt_allele, start, end):
        """Use the level tag to choose canonical transcript.
        Choose highest canonical score.
        The following order of preference is used:
        0.  membership in custom canonical transcript list
        1.  curation level
        2.  appris rank
        3.  most detrimental effect
        4.  longest protein change
        5.  lexicographical on transcript ID

        :param list txs: list of Transcript
        :param str variant_type:
        :param str ref_allele:
        :param str alt_allele:
        :param str start:
        :param str end:
        :return Transcript:
        """
        if len(txs) == 0:
            return None
        highest_scoring_tx = EnsemblTranscriptDatasource._select_best_with_multiple_criteria(txs, [
            (lambda x: self._is_tx_in_tx_list(x, self._custom_canonical_txs), max),
            (self._calculate_canonical_score, max),
            (self._get_appris_rank, min),
            (lambda x: self._calculate_effect_score(x, start, end, alt_allele, ref_allele, variant_type), min),
            (lambda x: len(x.get_seq()), max),
            (lambda x: x.get_transcript_id(), min)])
        return highest_scoring_tx[0]

    def get_overlapping_transcripts(self, chr, start, end, padding=0):
        new_start = str(int(start) - padding)
        new_end = str(int(end) + padding)
        records = self._get_binned_transcripts(chr, new_start, new_end)
        return self._get_overlapping_transcript_records(records, new_start, new_end)

    def get_overlapping_genes(self, chr, start, end):
        txs = self.get_overlapping_transcripts(chr, start, end)
        txs = self._filter_transcripts(txs)
        return set([tx.get_gene() for tx in txs])

    def _get_binned_transcripts_given_index(self, chr, start, end, index_dict):
        bins = region2bins(int(start), int(end))
        records = list()

        for b in bins:
            key = chr + "_" + str(b)
            try:
                txs = index_dict[key]
                records.extend(txs)
            except KeyError:
                pass
        return set(records)

    def _get_binned_genes(self, chr, start, end):
        return self._get_binned_transcripts_given_index(chr, start, end, self.gene_db)

    def _get_binned_transcripts(self, chr, start, end):
        return self._get_binned_transcripts_given_index(chr, start, end, self.gp_bin_db)

    def _get_overlapping_transcript_records(self, records, start, end):
        return [r for r in records if TranscriptProviderUtils.test_overlap(int(start), int(end), r.get_start(), r.get_end())]

    def _get_nearest_genes(self, chr, start, end):
        size_extensions = [1000, 10000, 100000, 1000000]

        left_gene, left_dist = None, None
        for s in size_extensions:
            new_start = start - s
            if new_start < 0: new_start = 1
            txs = self.get_transcripts_by_pos(chr, new_start, end)
            nearest_gene_border = 0
            for tx in txs:
                if tx.get_strand() == "-":
                    highest_genome_position = tx.determine_transcript_start()
                else:
                    highest_genome_position = tx.determine_transcript_stop()
                if highest_genome_position > nearest_gene_border:
                    nearest_gene_border = highest_genome_position
                    nearest_gene = tx.get_gene()
            if nearest_gene_border:
                left_dist = start - nearest_gene_border
                left_gene = nearest_gene
                break

        right_gene, right_dist = None, None
        for s in size_extensions:
            new_end = end + s
            txs = self.get_transcripts_by_pos(chr, start, new_end)
            nearest_gene_border = int(1e9)
            for tx in txs:
                if tx.get_strand() == "-":
                    lowest_genome_position = tx.determine_transcript_stop()
                else:
                    lowest_genome_position = tx.determine_transcript_start()
                if lowest_genome_position < nearest_gene_border:
                    nearest_gene_border = lowest_genome_position
                    nearest_gene = tx.get_gene()
            if nearest_gene_border < int(1e9):
                right_dist = nearest_gene_border - end
                right_gene = nearest_gene
                break

        return ((str(left_gene), str(left_dist)), (str(right_gene), str(right_dist)))

    def _render_other_transcripts(self, txs, transcriptIndicesToSkip, variant_type, ref_allele, alt_allele, start, end):
        """
        Create a list of transcripts that are not being chosen.

        Other transcripts are formatted <gene>_<transcript_id>_<variant_classification>_<protein_change>
            Note:  There are other areas of Oncotator (e.g. Generic_GeneProteinPositionDatasource) that depend
                on this format.  Changing it here may introduce bugs in other pieces of code.

                Also, do not include any transcript that would render as IGR.

        txs -- a list of transcripts to render.
        transcriptIndicesToSkip -- a list of transcripts that are being used (i.e. not an "other transcript").  This will usually be the canonical or any transcript chosen by tx_mode.
        """
        vcer = VariantClassifier()
        other_transcripts = list()
        for i, ot in enumerate(txs):
            if i not in transcriptIndicesToSkip:
                vc = vcer.variant_classify(tx=ot, variant_type=variant_type, ref_allele=ref_allele, alt_allele=alt_allele, start=start, end=end)
                if vc.get_vc() == VariantClassification.IGR:
                    continue
                o = '_'.join([ot.get_gene(), ot.get_transcript_id(),
                              vc.get_vc(), vcer.generate_protein_change_from_vc(vc)])
                o = o.strip('_')
                other_transcripts.append(o)

        return '|'.join(other_transcripts)

    def retrieve_transcripts_by_gene(self, gene):
        """
        Given a gene, return all of the transcripts (filtered) tagged with the gene.

        :param str gene:
        :rtype : list
        """
        txs_unfiltered = self.gene_db.get(gene, None)
        txs = self._filter_transcripts(txs_unfiltered)
        if txs is None:
            return []
        return txs

    def retrieveExons(self, gene, padding=10, isCodingOnly=False):
        """Return a list of (chr, start, end) tuples for each exon"""
        result = set()
        txs = self.gene_db.get(gene, None)
        if txs is None:
            return result
        txs = self._filter_transcripts(txs)

        for tx in txs:
            # If tx is coding
            if isCodingOnly and tx.get_gene_type() != "protein_coding":
                continue
            if isCodingOnly:
                exons = tx.get_cds()
            else:
                exons = tx.get_exons()

            for exon in exons:
                start = min(exon[0], exon[1])
                end = max(exon[0], exon[1])
                result.add((gene, tx.get_contig(), str(start - padding), str(end + padding)))

        return result

    def getTranscriptDict(self):
        return self.transcript_db

    def get_transcript(self, tx_id):
        if tx_id is None:
            return None
        return self.transcript_db.get(tx_id, None)

    def get_tx_mode(self):
        return self.tx_mode

    def _extract_segment_start_overlap(self, seg):
        """
        Given a segment, return the gene and exons (e.g. 6+) overlapping the start of the segment.

        :param MutationData seg:
         :return tuple:
            [0]: start_exon start exon index (0-based) and whether it is all previous exons
            ("-") or downstream exons ("+") in the coding direction.
             [1]: start_gene -- gene symbol on the canonical transcript
        """
        start_txs = self.get_transcripts_by_pos(chr=seg.chr, start=str(seg.start), end=str(seg.start))
        if start_txs is None or len(start_txs) == 0:
            start_gene = ""
            start_exon = ""
        else:
            start_chosen_tx = self._choose_transcript(start_txs, self.get_tx_mode(), VariantClassification.VT_SNP, "",
                                                      "", str(seg.start), str(seg.start))
            result_tuple = self._determine_exons_affected_by_start(int(seg.start), start_chosen_tx)
            start_gene = start_chosen_tx.get_gene()
            start_exon = str(result_tuple[0]) + result_tuple[1]
        return start_exon, start_gene

    def _extract_segment_end_overlap(self, seg):
        pos = seg.end
        end_txs = self.get_transcripts_by_pos(chr=seg.chr, start=str(pos), end=str(pos))
        if end_txs is None or len(end_txs) == 0:
            end_gene = ""
            end_exon = ""
        else:
            end_chosen_tx = self._choose_transcript(end_txs, self.get_tx_mode(), VariantClassification.VT_SNP, "",
                                                      "", str(pos), str(pos))
            result_tuple = self._determine_exons_affected_by_end(pos, end_chosen_tx)
            end_gene = end_chosen_tx.get_gene()
            end_exon = str(result_tuple[0]) + result_tuple[1]
        return end_exon, end_gene

    def annotate_segment(self, seg):
        """
        Akin to annotate_mutation, but for segments.

        Generates the following annotations:
        genes -- a comma-separated list of the genes found in a given region.
        start_gene -- gene symbol overlapped by the segment start position
        end_gene -- gene symbol overlapped by the segment end position
        start_exon -- exon overlap for the start gene.  Includes start exon index (0-based) and whether it is all previous exons
            ("-") or downstream exons ("+") in the coding direction.
            For example:
                6+ sixth exon and on
                6- sixth exon and previous

                Reminder that the exons are 0-based

        :returns MutationData seg: Annotated segment/region
        """
        txs = self.get_transcripts_by_pos(seg.chr, seg.start, seg.end)
        genes = set(([tx.get_gene() for tx in txs]))
        genes_annotation_value = ",".join(sorted(list(genes)))
        seg.createAnnotation("genes", genes_annotation_value, annotationSource=self.title, annotationDataType="String", annotationDescription="List of genes in the region.")

        # See if we can determine which gene and exon that is overlapped by start
        start_exon, start_gene = self._extract_segment_start_overlap(seg)

        seg.createAnnotation("start_gene", start_gene, annotationSource=self.title, annotationDataType="String", annotationDescription="Gene overlapping start of the region.")
        seg.createAnnotation("start_exon", start_exon, annotationSource=self.title, annotationDataType="String", annotationDescription="Exon index (0-based) that overlaps start with '+' or '-'.  '+' indicates all further exons of the gene are in the region.  '-' indicates all previous exons of the gene are in the region")

        end_exon, end_gene = self._extract_segment_end_overlap(seg)

        seg.createAnnotation("end_gene", end_gene, annotationSource=self.title, annotationDataType="String", annotationDescription="Gene overlapping end of the region.")
        seg.createAnnotation("end_exon", end_exon, annotationSource=self.title, annotationDataType="String", annotationDescription="Exon index (0-based) that overlaps end with '+' or '-'.  '+' indicates all further exons of the gene are in the region.  '-' indicates all previous exons of the gene are in the region")

        return seg

    def _extract_exon_info(self, position, tx):
        """
        Create basic information about the given position relative to the transcript.

        :param int position: in genomic space
        :param Transcript tx:
         :return tuple:
            [0]: closest exon index of the position (0-based),
             [1]: whether the distance was left in genomic space (false for overlap)
             [2]: whether the position overlaps an exon

        """
        exon_index = TranscriptProviderUtils.determine_closest_exon(tx, position, position)
        if exon_index is None:
            return exon_index, None, None, None
        left_distance, right_distance = TranscriptProviderUtils.determine_closest_distance_from_exon(position, position,
                                                                                                     exon_index, tx)
        is_in_exon = (left_distance <= 0) and (right_distance >= 0)
        is_diff_is_positive = (left_distance > 0) and (right_distance > 0)
        is_negative_strand = (tx.get_strand() == "-")
        return exon_index, is_diff_is_positive, is_in_exon, is_negative_strand

    def _determine_exons_affected_by_start(self, start, tx):
        """
        Return the exons affected by start position for the given transcript.
        The exon returned is always affected.

        Description (position is ___ of the nearest exon is )	l_diff vs r_diff	both are	seg	tx strand	result
        "to the right"	ldiff < rdiff	negative	start	"-"	< exon_i
        "to the left"	ldiff < rdiff	positive	start	"-"	<= exon_i
        "to the right"	ldiff < rdiff	negative	end	"-"	>= exon_i
        "to the left"	ldiff < rdiff	positive	end	"-"	> exon_i

        "to the right"	ldiff < rdiff	negative	start	"+"	> exon_i
        "to the left"	ldiff < rdiff	positive	start	"+"	>= exon_i
        "to the right"	ldiff < rdiff	negative	end	"+"	<= exon_i
        "to the left"	ldiff < rdiff	positive	end	"+"	< exon_i

        # When in an exon, include that exon as being affected
        "in an exon"    ldiff is negative, rdiff positive    both    start  "-" <= exon_i
        "in an exon"    ldiff is negative, rdiff positive    both    end "-" >= exon_i
        "in an exon"    ldiff is negative, rdiff positive    both    start  "+" >= exon_i
        "in an exon"    ldiff is negative, rdiff positive    both    end "+" <= exon_i

        :param int start: start position in genomic space
        :param Transcript tx: transcript affected by start and end
        :return tuple: Tuple of (exon_id, {"+","-"}) where the second is whether all higher number exons ("+") or lower number exons ("-")
        Please note that it is possible to get a result that should be interpreted as "no exon affected", such as (-1, "-")
        All "no exon affected" will have a [0] of -1.
        """
        exon_index, is_diff_is_positive, is_in_exon, is_negative_strand = self._extract_exon_info(int(start), tx)
        if exon_index is None:
            return None
        result_list = [-2, "X"]
        if is_in_exon and is_negative_strand:
            result_list = [exon_index, "-"]
        if is_in_exon and not is_negative_strand:
            result_list = [exon_index, "+"]

        # The rest are not in an exon
        if not is_in_exon:
            if is_negative_strand and is_diff_is_positive:
                result_list = [exon_index, "-"]
            if is_negative_strand and not is_diff_is_positive:
                result_list = [exon_index-1, "-"]
            if not is_negative_strand and is_diff_is_positive:
                result_list = [exon_index, "+"]
            if not is_negative_strand and not is_diff_is_positive:
                result_list = [exon_index + 1, "+"]

        if (result_list[0] < 0 and result_list[1] == "-") or (result_list[0] >= len(tx.get_exons()) and result_list[1] == "+"):
            # This gene is unaffected by the start position
            result_list[0] = -1

        return tuple(result_list)

    def _determine_exons_affected_by_end(self, end, tx):
        """
        Table of calculation is in the docs for _determine_exons_affected_by_start

        :param int end:
        :param Transcript tx:
        :return tuple: Tuple of (exon_id, {"+","-"}) where the second is whether all higher number exons ("+") or lower number exons ("-")
        Please note that it is possible to get a result that should be interpreted as "no exon affected", such as (-1, "-")
        All "no exon affected" will have a [0] of -1.
        """
        exon_index, is_diff_is_positive, is_in_exon, is_negative_strand = self._extract_exon_info(int(end), tx)
        if exon_index is None:
            return None

        result_list = [-2, "X"]
        if is_in_exon and is_negative_strand:
            result_list = [exon_index, "+"]
        if is_in_exon and not is_negative_strand:
            result_list = [exon_index, "-"]

        # The rest are not in an exon
        if not is_in_exon:
            if is_negative_strand and is_diff_is_positive:
                result_list = [exon_index + 1, "+"]
            if is_negative_strand and not is_diff_is_positive:
                result_list = [exon_index, "+"]
            if not is_negative_strand and is_diff_is_positive:
                result_list = [exon_index - 1, "-"]
            if not is_negative_strand and not is_diff_is_positive:
                result_list = [exon_index, "-"]

        if (result_list[0] < 0 and result_list[1] == "-") or (result_list[0] >= len(tx.get_exons()) and (result_list[1] == "+")):
            # This gene is unaffected by the end position
            result_list[0] = -1

        return tuple(result_list)

    def get_gene_symbols(self):
        """Return all of the gene symbols recognized by this datasource.
                """
        return self.gene_db.keys()

    def set_custom_canonical_txs(self, tx_list):
        """

        :param tx_list: list of transcript IDs, with or without version numbers.
        :return:
        """
        self._custom_canonical_txs = tx_list

    def get_custom_canonical_txs(self):
        return self._custom_canonical_txs

    def get_hashcode(self):
        """ Since this class can change annotation values depending on certain state attributes (e.g. tx-mode), we need the
        hashcode to change.  The super class hashcode attribute is treated like an initial hashcode here.

        In other words, hashcode is not a simple attribute for this datasource class.


        :return: hashcode including state information
        """
        hasher = Hasher()
        attrs_relevant_for_caching = [self.hashcode, self.get_tx_mode(), str(self._custom_canonical_txs)]
        for attr in attrs_relevant_for_caching:
            hasher.update(attr)
        return Hasher.md5_hash(hasher.hexdigest())
