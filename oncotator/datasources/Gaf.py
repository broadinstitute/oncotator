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

import cPickle
import logging
import os
import shove
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.datasources.Datasource import Datasource
from oncotator.datasources.GafInvalidChromosomeValue import GafInvalidChromosomeValue
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.index.gaf import region2bins
from oncotator.utils import gaf_annotation
from oncotator.utils.Hasher import Hasher
from oncotator.utils.MutUtils import MutUtils
from oncotator.utils.TagConstants import TagConstants
from oncotator.datasources.GafDatasourceException import GafDatasourceException


class Gaf(Datasource, TranscriptProvider):
    """
    Annotate from the Gaf 3.0.

    Annotations added to the mutations:
        transcripts -- a list of transcript dictionaries.  This annotation can be useful for this datasource, but generally should not be rendered.

        gene -- Hugo Symbol
        genome_change --
        annotation_transcript --
        transcript_strand --
        transcript_exon --
        transcript_position --
        cDNA_change --
        codon_change --
        protein_change --
        other_transcripts --
        variant_classification --


        Please see http://docs.sqlalchemy.org/en/rel_0_8/dialects/sqlite.html regarding multiprocess support for SQLite
    """

    def __init__(self, gaf_fname, gaf_transcript_sequences_fname, title='Gaf', version='3.0', tx_mode="CANONICAL", protocol="sqlite"):
        super(Gaf, self).__init__(src_file=gaf_fname, title=title, version=version)
        self.logger = logging.getLogger(__name__)
        if os.path.exists(gaf_fname):
            if not gaf_fname.endswith('.idx') and os.path.exists(gaf_fname + '.idx'):
                gaf_fname = gaf_fname + '.idx'
            else:
                raise Exception('Missing index for gaf file.  Index file with oncotator-index.py first.')
        else:
            raise Exception('Gaf file does not exist! -- %s' % gaf_fname)

        if os.path.exists(gaf_transcript_sequences_fname):
            if not gaf_transcript_sequences_fname.endswith('.idx') and os.path.exists(
                            gaf_transcript_sequences_fname + '.idx'):
                gaf_transcript_sequences_fname = gaf_transcript_sequences_fname + '.idx'
            else:
                raise Exception('Missing index for gaf file.  Index file with oncotator-index.py first.')
        else:
            raise Exception('Gaf transcript sequences file does not exist! -- %s' % gaf_transcript_sequences_fname)

        # 'Loading GAF...'
        self.logger.info("Loading GAF...")
        self.Transcripts, self.Genes = cPickle.load(open(gaf_fname, 'rb'))

        # 'Loading transcript sequences...'
        self.logger.info("Loading transcript sequences (" + protocol + ")...")
        self.gaf_transcript_sequences = shove.Shove(protocol + ':///%s' % gaf_transcript_sequences_fname, "memory://")

        # "Indexing Transcript IDs..."
        self.logger.info("Indexing transcript IDs...")
        self.transcript_id_idx = dict()
        for k in self.Transcripts:
            for b in self.Transcripts[k]:
                for i, t in enumerate(self.Transcripts[k][b]):
                    self.transcript_id_idx[t['transcript_id']] = (k, b, i)
        self.add_padding_to_GAF_transcripts(fiveprime_padding=3000, threeprime_padding=0)

        # "Indexing Gene IDs..."
        self.logger.info("Indexing gene IDs...")
        self.gene_id_idx = dict()
        for k in self.Genes:
            for b in self.Genes[k]:
                for i, t in enumerate(self.Genes[k][b]):
                    self.gene_id_idx[t['gene']] = (k, b, i)

        self.logger.info("Datasource " + self.title + " " + self.version + " finished initialization")

        # TODO: Check for valid values.
        self.tx_mode = tx_mode

    def get_transcript(self, tx_id):
        """Throws NotImplementedError """
        raise NotImplementedError("Cannot get Transcript objects from a GAF 3.0 datasource.")

    def get_transcripts_by_pos(self, chr, start, end):
        """ Return Transcript instances for the given chromosome start and end.
        Not implemented for GAF datasource.
        """
        raise NotImplementedError("Cannot get Transcript objects from a GAF 3.0 datasource.")

    def get_tx_mode(self):
        return self.tx_mode

    def set_tx_mode(self, value):
        self.tx_mode = value

    def get_hashcode(self):
        """The GAF datasource has to adjust  its key based on the internal tx mode.  set_hashcode sends in
         an initial hashcode, which is then adjusted by tx-mode

         """
        hasher = Hasher()
        hasher.update(self.hashcode)
        hasher.update(self.get_tx_mode())
        return hasher.hexdigest()

    def retrieveExons(self, gene, padding=10, isCodingOnly=False):
        """Return a list of (chr, start, end) tuples for each exon"""
        result = set()
        geneTuple = self.gene_id_idx.get(gene, None)
        if geneTuple is None:
            return result
        ctr = 0
        contig = MutUtils.convertChromosomeStringToMutationDataFormat(geneTuple[0])
        for b in self.Transcripts.get(contig, []):
            for i in self.Transcripts[contig][b]:
                if i['gene'] == gene:
                    if isCodingOnly and gaf_annotation.is_non_coding_transcript(i, self):
                        ctr += 1
                        continue

                    if isCodingOnly:
                        genomic_coords = self.getCodingTranscriptCoords(i)
                    else:
                        genomic_coords = i['genomic_coords']

                    for coord in genomic_coords:
                        start = min(coord[0], coord[1])
                        end = max(coord[0], coord[1])
                        result.add((gene, i['chr'], str(start - padding), str(end + padding)))
        return result

    def getCodingTranscriptCoords(self, t):
        """Gets the GAF exons, but adjusts to only include exons that code.  And truncates the UTR portions from GAF
         exons that do both."""
        exons = t['genomic_coords']
        t_coords = t['transcript_coords']
        cds_start = t['cds_start']
        cds_stop = t['cds_stop']
        strand = t['strand']
        is_neg_strand = (strand == "-")
        result = []
        t_start_coord_of_interest_idx = -1
        t_end_coord_of_interest_idx = -1

        # Get the start transcript_coord tuple
        for i in xrange(0,len(t_coords)):
            c = t_coords[i]

            # Strand is always positive, since the transcript_coords and cds_start and stop are in coding direction.
            overlap = gaf_annotation.test_overlap(c[0], c[1], cds_start, cds_stop, "+")
            if overlap == 'a_overlaps_b_left_border' or overlap == 'a_within_b' or overlap == 'a_encompasses_b':
                t_start_coord_of_interest_idx = i
                break

        # Get the end transcript_coord tuple
        for i in xrange(0,len(t_coords)):
            c = t_coords[i]

            # Strand is always positive, since the transcript_coords and cds_start and stop are in coding direction.
            overlap = gaf_annotation.test_overlap(c[0], c[1], cds_start, cds_stop, "+")
            if overlap == 'a_overlaps_b_right_border' or (overlap is None and c[0] > cds_stop) or overlap == 'a_encompasses_b':
                t_end_coord_of_interest_idx = i
                break
        start_t_coords = t_coords[t_start_coord_of_interest_idx]
        end_t_coords = t_coords[t_end_coord_of_interest_idx]
        codingTuples = exons[t_start_coord_of_interest_idx:t_end_coord_of_interest_idx+1]
        # if is_neg_strand:
        #     start_index = len(exons) - t_end_coord_of_interest_idx - 1
        #     end_index = len(exons) - t_start_coord_of_interest_idx
        #     codingTuples = exons[start_index:end_index]

        for tup in codingTuples:
            result.append([tup[0], tup[1]])


        if len(result) == 0:
            return result
        if not is_neg_strand:
            #Adjust start and end for overlap
            result[0][0] += (cds_start - start_t_coords[0])
            result[-1][1] -= (end_t_coords[1] - cds_stop)
        else:
            #Adjust start and end for overlap
            result[0][0] -= (cds_start - start_t_coords[0])
            result[-1][1] += (end_t_coords[1] - cds_stop)
        return result

    def getTranscriptDict(self):
        result = dict()
        for k in self.transcript_id_idx.keys():
            tid = self.transcript_id_idx[k]
            result[k] = self.Transcripts[tid[0]][tid[1]][tid[2]]
        return result

    def annotate_mutation(self, mutation, upstream_padding=3000, downstream_padding=0):
        mutation.createAnnotation('variant_type', TranscriptProviderUtils.infer_variant_type(mutation.ref_allele, mutation.alt_allele), self.title)
        data = [mutation]
        data = gaf_annotation.find_mut_in_gaf(data, self)
        data = gaf_annotation.identify_best_effect_transcript(data, self)
        data = gaf_annotation.identify_best_canonical_transcript(data, self)
        data = gaf_annotation.correct_transcript_coordinates(data, self)
        data = gaf_annotation.infer_output_fields(data, self)

        data = self._annotateMutationFromTranscripts(data)

        annotated_mutation = data.next()
        return annotated_mutation

    def _renderOtherTranscripts(self, m, transcriptIndicesToSkip):
        """
        Create a list of transcripts that are not being chosen.

        Other transcripts are formatted <gene>_<transcript_id>_<variant_classification>_<protein_change>
            Note:  There are other areas of Oncotator (e.g. Generic_GeneProteinPositionDatasource) that depend
                on this format.  Changing it here may introduce bugs in other pieces of code.

        m -- a mutation data object
        transcriptIndicesToSkip -- a list of transcripts that are being used (i.e. not an "other transcript").  This will usually be the canonical or any transcript chosen by tx_mode.
        """
        other_transcripts = list()
        for i, ot in enumerate(m['transcripts']):
            if i not in transcriptIndicesToSkip:
                o = '_'.join([ot['gene'], ot['transcript_id'],
                              ot['variant_classification'], ot.get('protein_change', '')])
                o = o.strip('_')
                other_transcripts.append(o)

        return '|'.join(other_transcripts)

    def _determineTranscriptIndices(self, m, tx_mode='CANONICAL'):
        if tx_mode == 'ALL':
            transcript_idxs = range(len(m['transcripts']))
            self.logger.warn('transcript mode of ALL is not supported yet.')
            raise NotImplementedError('transcript mode of ALL is not supported yet.')

        elif tx_mode == 'CANONICAL':
            transcript_idxs = [m['best_canonical_transcript']]
        elif tx_mode == 'EFFECT':
            transcript_idxs = [m['best_effect_transcript']]
        else:
            raise NotImplementedError("Cannot get here.")

        # Convert to int
        result = []
        for idx in transcript_idxs:
            result.append(int(idx))
        return result

    def _determineTranscriptPosition(self, transcript):
        absolute_transcript_pos = [transcript.get('absolute_transcript_position_start', '')]
        if transcript.get('absolute_transcript_position_start') != transcript.get('absolute_transcript_position_end'):
            absolute_transcript_pos.append(transcript.get('absolute_transcript_position_end'))
        absolute_transcript_pos = '_'.join(map(str, absolute_transcript_pos))
        return absolute_transcript_pos

    def _annotateMutationFromTranscripts(self, muts):
        """ Assumes that the transcripts annotation has already been populated. """
        #TODO: We should remove the transcript annotation.  This should not be necessary anymore, since the datasources are all encapsulated from each other.  What is worse, the transcript annotation may be blocking multiprocess support
        for m in muts:
            if (m['transcripts'] is None) or (m['transcripts'] == ''):
                errStr = "Transcript annotation was blank for mutation: " + str(
                    m) + ".  Was it populated before calling _annotateMutationFromTranscripts?"
                self.logger.error(errStr)
                raise GafDatasourceException(errStr)
            if len(m['transcripts']) == 0:
                self.logger.warn('No transcripts found for mutation: ' + m.chr + " " + m.start + " " + m.end)
            transcriptIndices = self._determineTranscriptIndices(m, self.tx_mode)

            if len(transcriptIndices) > 1:
                self.logger.warn('More than one transcript found for annotation.  This is not supported yet.')
            if len(transcriptIndices) == 0:
                self.logger.warn(
                    'No transcript found for mutation.  This is not supported yet.  A failure may be imminent.')

            # TODO: Support more than one transcript at a time.  This would support tx-mode=ALL
            transcriptIndex = transcriptIndices[0]
            transcript = m['transcripts'][transcriptIndex]

            if transcript['variant_classification'] == 'IGR':
                m.createAnnotation('gene', 'Unknown', self.title, tags=[TagConstants.NOT_SPLIT])
                m.createAnnotation('gene_id', "0", self.title, tags=[TagConstants.NOT_SPLIT])
                m.createAnnotation('other_transcripts', transcript['nearest_genes'], self.title, tags=[TagConstants.SPLIT])

                #TODO: Refactor to reduce duplicate code
                m.createAnnotation('annotation_transcript', '', self.title, tags=[TagConstants.NOT_SPLIT])
                m.createAnnotation('codon_change', '', self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('strand', '', self.title, tags=[TagConstants.NOT_SPLIT])
                m.createAnnotation('protein_change', '', self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('transcript_exon', '', self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('transcript_position', '', self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('transcript_change', '', self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('transcript_id', '', self.title, tags=[TagConstants.NOT_SPLIT]),
                m.createAnnotation('transcript_strand', '', self.title, tags=[TagConstants.NOT_SPLIT])

            else:
                m.createAnnotation('gene', transcript['gene'], self.title, tags=[TagConstants.NOT_SPLIT])
                m.createAnnotation('gene_id', transcript['id'], self.title, tags=[TagConstants.NOT_SPLIT])
                m.createAnnotation('annotation_transcript', transcript['transcript_id'], self.title, tags=[TagConstants.NOT_SPLIT])
                m.createAnnotation('transcript_exon', transcript['exon_affected'], self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('transcript_position', self._determineTranscriptPosition(transcript), self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('transcript_change', transcript.get('transcript_change', ''), self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('transcript_id', transcript.get('transcript_id', ''), self.title, tags=[TagConstants.NOT_SPLIT]),
                m.createAnnotation('transcript_strand', transcript.get('strand', ''), self.title, tags=[TagConstants.NOT_SPLIT]),
                m.createAnnotation('codon_change', transcript['codon_change'], self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('protein_change', transcript['protein_change'], self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('other_transcripts', self._renderOtherTranscripts(m, transcriptIndices), self.title, tags=[TagConstants.SPLIT])
                m.createAnnotation('strand', transcript['strand'], self.title, tags=[TagConstants.NOT_SPLIT])
            m.createAnnotation('tx_start', transcript.get('tx_start', ''), self.title, tags=[TagConstants.SPLIT])
            m.createAnnotation('tx_end', transcript.get('tx_end', ''), self.title, tags=[TagConstants.SPLIT])
            m.createAnnotation('variant_classification', transcript['variant_classification'], self.title, tags=[TagConstants.SPLIT])
            m.createAnnotation('transcript_description', transcript.get('description', ''), self.title, tags=[TagConstants.SPLIT])
            m.createAnnotation('transcript_protein_position_start', str(transcript.get('protein_position_start', '')), self.title, tags=[TagConstants.SPLIT])
            m.createAnnotation('transcript_protein_position_end', str(transcript.get('protein_position_end', '')), self.title, tags=[TagConstants.SPLIT])

            # TODO: Low Priority -- rename transcripts to _transcripts instead of deleting
            del m['transcripts']
            yield m

    def determineGenomeChange(self, variant_type, chrom, start, end, ref, alt):
        """ Returns genome_change as a string.

        Inputs:
        variant_type -- string.  Usually 'SNP', 'INS', or 'DEL'.  Though xNP is also supported.
        chrom -- chromosome as string
        start -- start position as integer
        end -- end position as integer
        ref -- reference allele as string
        alt -- alternate/observed allele as string.
    """
        genome_change = ''

        if variant_type == 'SNP':
            genome_change = 'g.chr%s:%d%s>%s' % (chrom, start, ref, alt)
        elif variant_type.endswith('NP'):
            genome_change = 'g.chr%s:%d_%d%s>%s' % (chrom, start, end, ref, alt)
        elif variant_type == 'DEL':
            if start == end:
                genome_change = 'g.chr%s:%ddel%s' % (chrom, start, ref)
            else:
                genome_change = 'g.chr%s:%d_%ddel%s' % (chrom, start, end, ref)
        elif variant_type == 'INS':
            genome_change = 'g.chr%s:%d_%dins%s' % (chrom, start, end, alt)
        return genome_change

    def close(self):
        self.gaf_transcript_sequences.close()

    def get_transcript_seq(self, transcript_id):
        return self.gaf_transcript_sequences.get(transcript_id, None)

    def get_protein_seq(self, transcript_id):
        gaf_record = self.get_transcript(transcript_id)
        tx_seq = self.get_transcript_seq(transcript_id)
        if not gaf_record or not tx_seq:
            return None

        if 'cds_start' not in gaf_record or not gaf_record['cds_start']:
            return None

        prot_seq = MutUtils.translate_sequence(tx_seq[gaf_record['cds_start']-1:gaf_record['cds_stop']])
        if prot_seq[-1] == '*':
            prot_seq = prot_seq[:-1]

        return prot_seq

    def get_overlapping_transcripts(self, chr, start, end):
        records = self.__get_binned_data(chr,start, end, 'transcript')
        return self.__get_overlapping_records(records, start, end, 'transcript')

    def get_overlapping_genes(self, chr, start, end):
        records = self.__get_binned_data(chr,start, end, 'gene')
        return self.__get_overlapping_records(records, start, end, 'gene')

    def get_gene(self, symbol):
        try:
            keys = self.gene_id_idx[symbol]
            return self.Genes[keys[0]][keys[1]][keys[2]]
        except KeyError:
            return {}

    def get_transcript(self, id):
        try:
            keys = self.transcript_id_idx[id]
            return self.Transcripts[keys[0]][keys[1]][keys[2]]
        except KeyError:
            return {}

    def get_nearest_genes(self, chr, start, end):
        if not self.Genes:
            return ((str(None), str(None)), (str(None), str(None)))

        size_extensions = [1000, 10000, 100000, 1000000]

        left_gene, left_dist = None, None
        for s in size_extensions:
            new_start = start - s
            if new_start < 0: new_start = 1
            records = self.get_overlapping_genes(chr, new_start, end)
            neareast_gene_border = 0
            for r in records:
                if r['end'] > neareast_gene_border:
                    neareast_gene_border = r['end']
                    nearest_gene = r['gene']
            if neareast_gene_border:
                left_dist = start - neareast_gene_border
                left_gene = nearest_gene
                break

        right_gene, right_dist = None, None
        for s in size_extensions:
            new_end = end + s
            records = self.get_overlapping_genes(chr, start, new_end)
            neareast_gene_border = int(1e9)
            for r in records:
                if r['start'] < neareast_gene_border:
                    neareast_gene_border = r['start']
                    nearest_gene = r['gene']
            if neareast_gene_border < int(1e9):
                right_dist = neareast_gene_border - end
                right_gene = nearest_gene
                break

        return ((str(left_gene), str(left_dist)), (str(right_gene), str(right_dist)))

    def add_padding_to_GAF_transcripts(self, fiveprime_padding, threeprime_padding):
        for chr in self.Transcripts:
            for bin in self.Transcripts[chr]:
                for t in self.Transcripts[chr][bin]:
                    if t['strand'] == '+':
                        t['footprint_start'] = t['footprint_start'] - fiveprime_padding
                        t['footprint_end'] = t['footprint_end'] + threeprime_padding
                    elif t['strand'] == '-':
                        t['footprint_start'] = t['footprint_start'] - threeprime_padding
                        t['footprint_end'] = t['footprint_end'] + fiveprime_padding

    def __get_overlapping_records(self, records, start, end, type):
        if type == 'gene':
            st_key, en_key = 'start', 'end'
        elif type == 'transcript':
            st_key, en_key = 'footprint_start', 'footprint_end'

        out_records = list()
        for r in records:
            if TranscriptProviderUtils.test_overlap(start, end, r[st_key], r[en_key]):
                out_records.append(r)

        return out_records

    def __get_binned_data(self, chr, start, end, type):
        if type == 'gene':
            data = self.Genes
        elif type == 'transcript':
            data = self.Transcripts

        # GAF uses M_rCRS whereas mutations are often using M for the chromosome field.
        if chr == 'M' and ('M' not in data) and ('M_rCRS' not in data):
            raise GafInvalidChromosomeValue("Unable to process mitochondria mutation with chr: %s   .... data keys are %s" % (str(chr)), str(data.keys()))
        if chr == 'M' and ('M' not in data):

            # TODO: Verify that it is okay to do this.
            chr = 'M_rCRS'

        if chr not in data:
            self.logger.warn("Invalid chromosome value for Gaf search: %s" % (str(chr)))
            return list()
            #raise GafInvalidChromosomeValue("Invalid chromosome value: %s" % (str(chr)))

        bins = region2bins(start, end)
        records = list()
        for b in bins:
            records.extend(data[chr].get(b, []))

        return records