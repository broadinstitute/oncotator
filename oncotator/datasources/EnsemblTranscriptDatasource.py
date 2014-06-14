# LICENSE_GOES_HERE
import logging
import shove
from oncotator.Annotation import Annotation
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources.Datasource import Datasource
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.index.gaf import region2bins
from oncotator.utils.HgvsChangeTransformer import HgvsChangeTransformer
from oncotator.utils.VariantClassification import VariantClassification
from oncotator.utils.VariantClassifier import VariantClassifier
from oncotator.utils.txfilter.TranscriptFilterFactory import TranscriptFilterFactory



class EnsemblTranscriptDatasource(TranscriptProvider, Datasource):
    """ Similar to a GAF datasource, but uses ensembl transcripts.
        Also, supports gencode

        Though all transcripts for GENCODE can be loaded, it is currently set to ignore any transcripts that are not
            "basic"
    """
    """This is the list of annotations that get populated by this datasource"""
    POPULATED_ANNOTATION_NAMES = set(['transcript_exon', 'variant_type', 'variant_classification', 'other_transcripts', 'gene', 'gene_id', 'annotation_transcript', 'genome_change', 'strand', 'transcript_id', 'secondary_variant_classification', 'protein_change', 'codon_change', 'transcript_change', 'transcript_strand', 'gene', 'gene_type', 'gencode_transcript_tags', 'gencode_transcript_status', 'havana_transcript', 'ccds_id', 'gencode_transcript_type', 'transcript_position', 'gencode_transcript_name'])
    def __init__(self,  src_file, title='ENSEMBL', version='', tx_mode=TranscriptProvider.TX_MODE_CANONICAL, protocol="file", is_thread_safe=False, tx_filter="dummy"):
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

        logging.getLogger(__name__).info("%s %s is being set up with default tx-mode: %s.  " % (title, version, tx_mode))
        self.set_tx_mode(tx_mode)

        logging.getLogger(__name__).info("%s %s is being set up with %s filtering.  " % (title, version, tx_filter))
        self._tx_filter = TranscriptFilterFactory.create_instance(tx_filter)

        self._hgvs_xformer = HgvsChangeTransformer()

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
        """Given a list of transcripts and a transcript mode (e.g. CANONICAL), choose the transcript to use. """
        if len(txs) == 1:
            return txs[0]
        if tx_mode == TranscriptProvider.TX_MODE_CANONICAL:
            return self._choose_canonical_transcript(txs, variant_type, ref_allele, alt_allele, start, end)
        return self._choose_best_effect_transcript(txs, variant_type, ref_allele, alt_allele, start, end)

    def _choose_best_effect_transcript(self, txs, variant_type, ref_allele, alt_allele, start, end):
        """Choose the transcript with the most detrimental effect.
         The rankings are in TranscriptProviderUtils.
         Ties are broken by which transcript has the longer coding length.
         """
        vcer = VariantClassifier()
        effect_dict = TranscriptProviderUtils.retrieve_effect_dict()
        best_effect_score = 100000000 # lower score is more likely to get picked
        best_effect_tx = None
        for tx in txs:
            vc = vcer.variant_classify(tx, ref_allele, alt_allele, start, end, variant_type).get_vc()
            effect_score = effect_dict.get(vc, 25)
            if effect_score < best_effect_score:
                best_effect_score = effect_score
                best_effect_tx = tx
            elif (effect_score == best_effect_score) and (len(best_effect_tx.get_seq()) < len(tx.get_seq())):
                best_effect_score = effect_score
                best_effect_tx = tx

        return best_effect_tx

    def _choose_canonical_transcript(self, txs, variant_type, ref_allele, alt_allele, start, end):
        """Use the level tag to choose canonical transcript.

        Choose highest canonical score.

        Ties are broken by whichever transcript has the most detrimental effect.
        """
        if len(txs) == 0:
            return None
        scores = dict()
        for tx in txs:
            score = self._calculate_canonical_score(tx)
            if score not in scores.keys():
                scores[score] = set()
            scores[score].add(tx)

        highest_score = max(scores.keys())
        highest_scoring_txs = scores[highest_score]
        if len(highest_scoring_txs) == 1:
            return list(highest_scoring_txs)[0]
        else:
            return self._choose_best_effect_transcript(highest_scoring_txs, variant_type, ref_allele, alt_allele, start, end)

    def _calculate_canonical_score(self, tx):
        """
        Level 1 is validated
        Level 2 is manual annotation
        Level 3 is automated annotation.
        :param tx: Transcript
        :return:
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
