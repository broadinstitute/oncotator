"""
# By downloading the PROGRAM you agree to the following terms of use:
#
# BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
# FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
#
# This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
# WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
# WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
# NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
#
# 1. DEFINITIONS
# 1.1	"PROGRAM" shall mean copyright in the object code and source code known as Oncotator and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.
#
# 2. LICENSE
# 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
# LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
# The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
# 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
# 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
#
# 3. OWNERSHIP OF INTELLECTUAL PROPERTY
# LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
#
# Copyright 2012 Broad Institute, Inc.
# Notice of attribution:  The Oncotator program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc.
#
# LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
#
# 4. INDEMNIFICATION
# LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
#
# 5. NO REPRESENTATIONS OR WARRANTIES
# THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
# IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
#
# 6. ASSIGNMENT
# This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
#
# 7. MISCELLANEOUS
# 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
# 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
# 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
# 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
# 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
# 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
# 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
#"""


import abc
import cPickle
from operator import mod
import os
import shove
from Bio import Seq
from shove.core import Shove
from oncotator.Annotation import Annotation
from oncotator.MissingAnnotationException import MissingAnnotationException
from oncotator.TranscriptProviderUtils import TranscriptProviderUtils
from oncotator.utils import gaf_annotation
from oncotator.index.gaf import region2bins
import vcf
import logging
import re
import copy
from oncotator.utils.Hasher import Hasher
from oncotator.utils.MutUtils import MutUtils
from oncotator.utils.VariantClassifier import VariantClassifier
from oncotator.utils.gaf_annotation import GAFNonCodingTranscript
from oncotator.utils.TagConstants import TagConstants
from oncotator.utils.dbNSFP import dbnsfp_fieldnames

try:
    import pysam
except ImportError:
    warningString = (
        'ERROR: Could not load pysam.  Some features will be disabled (e.g. COSMIC annotations) and may cause Oncotator to fail.')
    print(warningString)
    logging.getLogger('').error(warningString)

import collections
from GafDatasourceException import GafDatasourceException
from oncotator.utils.db import get_db_data, get_binned_data, get_overlapping_records, get_summary_output_string


class GafInvalidChromosomeValue(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class TranscriptProvider(object):

    TX_MODE_CANONICAL = "CANONICAL"
    TX_MODE_BEST_EFFECT = "EFFECT"
    TX_MODE_CHOICES = [TX_MODE_CANONICAL, TX_MODE_BEST_EFFECT] # Simply list the ones above here.

    @abc.abstractmethod
    def getTranscriptDict(self):
        """ Return a dict containing all transcripts where key is the transcript ID.
        """
        return

    @abc.abstractmethod
    def retrieveExons(self, gene, padding=10, isCodingOnly=False):
        """Return a list of (chr, start, end) tuples for each exon in each transcript"""
        return

    @abc.abstractmethod
    def set_tx_mode(self, tx_mode):
        # TODO: Throw exception if not in TX_MODE_CHOICES
        return

    @abc.abstractmethod
    def get_tx_mode(self):
        return self.tx_mode


class Datasource(object):
    """
    An individual datasource used for annotation.  This serves as a base class for attributes
    and methods that are common to all types of datasources.  Subclasses of Datasource will
    define behavior for more specific data types.
    
    src_file
        Absolute path to the annotation datasource file.
        
    title
        Title string to be appended to each datasource header.
        
    version
        Version of datasource. Will be displayed in comments of output file.


    Note for developers:  When creating new types of datasources, you may need to be aware of the sorting position.
    For example, datasources that rely on the gene annotation must be placed after the GAF (or other transcript
        datasource).  See DatasourceCreator::sortDatasources()
    """

    def __init__(self, src_file, title='', version=None, default_missing_value=''):
        """
        This method is used to load/parse the annotation datasource.  Should be overridden
        for custom parsing/loading
        
        The title option should be used to label output headers. Please see the Generic_Gene_DataSource
        and Generic_GenomicPosition_DataSource Datasource subclasses as examples.
        
        If added annotation fields are to be outputted by AnnotationSet.write_data method,
        then a output_headers list must be instantiated.  Please see the Generic_Gene_DataSource
        and Generic_GenomicPosition_DataSource Datasource subclasses as examples.
        
        TODO: Update this documentation
        
        """
        self.src_file = src_file
        self.title = title
        self.version = str(version)
        self.output_headers = []
        self.hashcode = ""

    def attach_to_class(self, cls, name, options):
        """
        Do not touch.  This method is necessary for the AnnotationSetMeta class to work.
        
        This coded lets AnnotationSet become aware of their constituent Datasource classes.
        """
        self.cls = cls
        self.name = name
        self.options = options
        if self.title is None:
            self.title = name
        options.add_datasource(self)

    @abc.abstractmethod
    def annotate_mutation(self, mutation):
        # The default implementation raise a NotImplementedError
        # to ensure that any subclasses must override this method.
        raise NotImplementedError

    def set_hashcode(self, value):
        self.hashcode = value

    def get_hashcode(self):
        return self.hashcode

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
        for b in self.Transcripts[geneTuple[0]]:
            for i in self.Transcripts[geneTuple[0]][b]:
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

        prot_seq = Seq.translate(tx_seq[gaf_record['cds_start']-1:gaf_record['cds_stop']])
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

class GenericGeneDataSourceException(Exception):
    def __init__(self, str):
        """
        
        """
        self.value = str

    def __str__(self):
        return repr(self.value)

class Generic_Gene_DataSource(Datasource):
    """
    A datasource derived from a generic TSV file in which the first column is a HUGO gene 
    symbol.  First header value must be specified in the constructor (default: 'gene').  All other columns will be used for
    annotation.
    TODO: This is no longer true.  Actually, you can specify any column as the gene column in the config file.  Update documentation.
    use_binary
        if True, existing indexed binary will be used or created for future use.
        
    """
    def __init__(self, src_file, title='', version=None, use_binary=True, geneColumnName='gene'):
        super(Generic_Gene_DataSource, self).__init__(src_file, title=title, version=version)
        
        index_mode = 'gene'
        
        self.db_obj, self.output_headers = get_db_data(src_file, title, use_binary, index_mode,indexColumnNames=geneColumnName)
            
    def annotate_mutation(self, mutation, index_field='gene'):
        
        if index_field not in mutation:
            raise GenericGeneDataSourceException("Index field (" + index_field + ") not found.  Remember that datasources must be ordered.  Please put a datasource that provides the 'gene' annotation in front of this datasource.")
        
        #if any([c in mutation for c in self.output_headers]):
        for c in self.output_headers:
            if c in mutation:
                raise Exception('Error: Non-unique header value in annotation table (%s)' % (c))
        
        gene = mutation[index_field]
        if gene in self.db_obj:
            annotations = self.db_obj[gene]
            for k in annotations.keys():
                mutation.createAnnotation(k, self.db_obj[gene][k], annotationSource=self.title)
        else:
            for h in self.output_headers:
                mutation.createAnnotation(h, '', annotationSource=self.title)

        return mutation

class Generic_VariantClassification_Datasource(Generic_Gene_DataSource):
    """ Used for generic TSV that is indexed by variant classification. """
    def __init__(self, src_file, title='', version=None, use_binary=True, geneColumnName='variant_classification'):
        super(Generic_VariantClassification_Datasource,self).__init__(src_file, title, version, use_binary, geneColumnName)

    def annotate_mutation(self, mutation):
        return super(Generic_VariantClassification_Datasource,self).annotate_mutation(mutation,'variant_classification')

class Generic_Transcript_Datasource(Generic_Gene_DataSource):
    """ Used for generic TSV that is indexed by transcript ID. """
    def __init__(self, src_file, title='', version=None, use_binary=True, geneColumnName='transcript_id'):
        super(Generic_Transcript_Datasource,self).__init__(src_file, title, version, use_binary, geneColumnName)

    def annotate_mutation(self, mutation):
        return super(Generic_Transcript_Datasource,self).annotate_mutation(mutation,'transcript_id')
   
        
class Generic_GenomicPosition_DataSource(Datasource):
    """
    A datasource derived from a generic TSV file in which the first three columns are 'chr',
    'start', and 'end'.  All other columns will be used for annotation.
    
    use_binary
        if True, existing indexed binary will be used or created for future use.
        
    """
    def __init__(self, src_file, title='', version=None, use_binary=True, gpColumnNames="chr,start,end"):
        super(Generic_GenomicPosition_DataSource, self).__init__(src_file, title=title, version=version)
        
        index_mode = 'genomic_pos'
        self.db_obj, self.output_headers = get_db_data(src_file, title, use_binary, index_mode, gpColumnNames)
    
    def annotate_mutation(self, mutation):
        #if any([c in mutation for c in self.output_headers]):
        for c in self.output_headers:
            if c in mutation:
                raise Exception('Error: Non-unique header value in annotation table (%s)' % (c))

        if all(field in mutation for field in ['chr','start','end']):
            chr, start, end = mutation.chr, mutation.start, mutation.end
                
            records = get_binned_data(self.db_obj, chr,int(start), int(end))
            records = get_overlapping_records(records, int(start), int(end))
            if records:
#                for r in records:
                for c in self.output_headers:
                    summarized_results = get_summary_output_string([r[c].strip() for r in records])
                    mutation.createAnnotation(c, summarized_results, annotationSource=self.title)
        
        for header in self.output_headers:
            if header not in mutation:
                mutation.createAnnotation(header, '', annotationSource=self.title)        
        
        return mutation

class Generic_GenomicMutation_Datasource(Generic_GenomicPosition_DataSource):
    """
    This datasource extends Generic_GenomicPosition_DataSource to also match on ref_allele
    and alt_allele columns.  All other columns will be used for annotation.
    """
    def __init__(self, src_file, title='', version=None, **kwargs):
        super(Generic_GenomicMutation_Datasource, self).__init__(src_file, title=title, version=version, **kwargs)

        self.ref_allele_fieldname, self.alt_allele_fieldname = None, None
        for header in self.output_headers:
            if header.endswith('_ref_allele'):
                self.ref_allele_fieldname = header
            if header.endswith('_alt_alele'):            
                self.alt_allele_fieldname = header

        if self.ref_allele_fieldname is None or self.alt_allele_fieldname is None:
            raise Exception('Unable to determine ref_allele or alt_allele column name.')

        self.output_headers = [h for h in self.output_headers if not h.endswith('_ref_allele') and not h.endswith('_alt_alele')]

    def annotate_mutation(self, mutation):
        #if any([c in mutation for c in self.output_headers]):
        for c in self.output_headers:
            if c in mutation:
                raise Exception('Error: Non-unique header value in annotation table (%s)' % (c))

        if all(field in mutation for field in ['chr','start','end', 'ref_allele', 'alt_allele']):
            chr, start, end = mutation.chr, mutation.start, mutation.end
            ref_allele, alt_allele = mutation.ref_allele, mutation.alt_allele
                
            records = get_binned_data(self.db_obj, chr,int(start), int(end))
            records = get_overlapping_records(records, int(start), int(end))
            records = [rec for rec in records if rec[self.ref_allele_fieldname] == ref_allele and rec[self.alt_allele_fieldname] == alt_allele]
            if records:
                for c in self.output_headers:
                    summarized_results = get_summary_output_string([r[c].strip() for r in records])
                    mutation.createAnnotation(c, summarized_results, annotationSource=self.title)
        
        for header in self.output_headers:
            if header not in mutation:
                mutation.createAnnotation(header, '', annotationSource=self.title)        
        
        return mutation


class IndexedVCF_DataSource(Datasource):
    """
    A datasource derived from a VCF file.  Expects a bgzipped vcf using Tabix.
    
    Instructions on how to index file using Tabix prior Oncotator:
    
    bgzip foo.vcf
    tabix -p vcf foo.vcf.gz

    Please see http://samtools.sourceforge.net/tabix.shtml for more info.
    
    
    The following VCF columns will be added to the output.
        CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO
    
    
    Multiple annotation data for the same mutation will be delimited by "|".
    
    TODO: Support for FORMAT and sample columns in output
    
    """
    def __init__(self, src_file, title='', version=None):
        super(IndexedVCF_DataSource, self).__init__(src_file, title=title, version=version)  # all of the info is coming from config file

        self.vcf_reader = vcf.Reader(filename=src_file, strict_whitespace=True)
        self.vcf_headers = self.vcf_reader.infos.keys()
        self.output_vcf_headers = dict([(header, '_'.join([self.title, header])) for header in self.vcf_headers])
        self.logger = logging.getLogger(__name__)

    def _determine_tags(self):
        """
        Determine tag constants associated with annotation IDs.

        :return: map of tag IDs and respective tag constants
        """
        tagsMap = {}
        for ID in self.vcf_headers:
            num = self.vcf_reader.infos[ID].num
            if num == -1:
                tagsMap[ID] = [TagConstants.INFO, TagConstants.SPLIT]
            else:
                tagsMap[ID] = [TagConstants.INFO, TagConstants.NOT_SPLIT]
        return tagsMap

    def _determine_annotationValue(self, vcf_record, ID):
        """
        Determine the appropriate value that corresponds to a given ID.
        This method checks whether ID is present in the input Vcf record or not. If it is not found, then an appropriate
        missing value is computed.

        :param vcf_record: input Vcf record
        :param ID: ID
        :return: value that corresponds to a given ID
        """
        if ID in vcf_record.INFO:
            vals = vcf_record.INFO[ID]
            if isinstance(vals, list):
                val = ",".join(map(str, ["" if val is None else val for val in vals]))
            else:
                val = str(vals)
        else:
            val = self._determine_missing_annotationValue(vcf_record, ID)
        return val

    def _determine_missing_annotationValue(self, vcf_record, ID):
        """
        Determine the missing value that corresponds to a given ID.
        This method takes into account the number field and computes an appropriate missing value for a given ID.

        :param vcf_record: input Vcf record used to determine the number
        :param ID: ID
        :return: value that corresponds to a given ID
        """
        nsamples = len(self.vcf_reader.samples)
        num = self.vcf_reader.infos[ID].num
        if num == -2:
            val = ",".join([""]*nsamples)
        elif num == -1:
            val = ""
            if vcf_record is not None:
                val = ",".join([""]*len(vcf_record.ALT))
        elif num == 0:
            val = "False"
        elif num is None:
            val = ""
        else:
            val = ",".join([""]*num)
        return val

    def annotate_mutation(self, mutation):
        """
        Annotate mutation with appropriate annotation value pairs from Tabix indexed Vcf file.

        :param mutation: mutation to annotate
        :return: annotated mutation
        """
        vcf_records = None
        try:
            vcf_records = self.vcf_reader.fetch(mutation.chr, int(mutation.start), int(mutation.end))
        except ValueError as ve:
            self.logger.warn("Exception when looking for vcf records. Empty set of records being returned: " +
                             repr(ve))

        tagsMap = self._determine_tags()
        valsMap = dict()

        if vcf_records is not None:
            for vcf_record in vcf_records:
                for ID in self.vcf_headers:
                    val = self._determine_annotationValue(vcf_record, ID)
                    if ID not in valsMap:
                        valsMap[ID] = [val]
                    else:
                        valsMap[ID] += [val]

        if len(valsMap) == 0:
            for ID in self.vcf_headers:
                vcf_record = None
                val = self._determine_missing_annotationValue(vcf_record, ID)
                valsMap[ID] = [val]

        for ID in self.vcf_headers:
            # multiple values are delimited by "|"
            val = "|".join(valsMap[ID])
            tags = copy.copy(tagsMap[ID])
            mutation.createAnnotation(self.output_vcf_headers[ID], val, self.title, self.vcf_reader.infos[ID].type,
                                      self.vcf_reader.infos[ID].desc, tags=tags, number=self.vcf_reader.infos[ID].num)
        return mutation
        


class IndexedTSV_Datasource(Datasource):
    """
    A datasource derived from a TSV file.  Expects a bgzipped vcf using Tabix.

    Instructions on how to index file using Tabix prior Oncotator:

    bgzip foo.vcf
    tabix -p vcf foo.vcf.gz

    Please see http://samtools.sourceforge.net/tabix.shtml for more info.

    Multiple annotation data for the same mutation will be delimited by "|".
    """
    def __init__(self, src_file, title='', version=None, colnames=[], indexColnames=[]):
        super(IndexedTSV_Datasource, self).__init__(src_file, title=title, version=version)

        self.tsv_reader = pysam.Tabixfile(filename=src_file)
        self.tsv_headers = dict([(x, i) for (i, x) in enumerate(colnames)])  # index of the tsv header
        self.output_tsv_headers = dict([(colname, '_'.join([self.title, colname])) for colname in colnames])
        for colname in indexColnames:
            if colname in self.output_tsv_headers:
                self.output_tsv_headers.pop(colname, None)
        self.logger = logging.getLogger(__name__)

    def annotate_mutation(self, mutation):
        """
        Annotate mutation with appropriate annotation value pairs from a Tabix indexed TSV file.

        :param mutation: mutation to annotate
        :return: annotated mutation
        """
        valsMap = dict()
        try:
            tsv_records = self.tsv_reader.fetch(mutation.chr, int(mutation.start)-1, int(mutation.end),
                                                parser=pysam.asTuple())
            for tsv_record in tsv_records:
                for colname in self.output_tsv_headers:
                    val = ""
                    if tsv_record is not None:
                        val = tsv_record[self.tsv_headers[colname]]

                    # multiple records are delimited by "|"
                    if colname not in valsMap:
                        valsMap[colname] = val
                    else:
                        valsMap[colname] += "|" + val
        except ValueError as ve:
            msg = "Exception when looking for tsv records. Empty set of records being returned: " + repr(ve)
            self.logger.warn(msg)

        for colname in self.output_tsv_headers:
            val = ""
            if len(valsMap) != 0:
                val = valsMap[colname]
            else:
                msg = "Exception when looking for tsv records. Empty set of records being returned."
                self.logger.warn(msg)

            mutation.createAnnotation(self.output_tsv_headers[colname], val, self.title,
                                      tags=[TagConstants.INFO, TagConstants.NOT_SPLIT])

        return mutation


class dbSNP(Datasource):
    """
    dbSNP Datasource expecting a Tabix compressed and indexed dbsnp vcf file.
    
    Output columns are:

        dbSNP_RS
            RS id of overlapping dbSNPs.  Multiple records are delimited by '|'
        dbSNP_Val_Status
            'byFrequency' if dbSNP is validated by frequency.
            'b10000genomes' if dbSNP is found in 1000genomes dataset.
            
    """
    def __init__(self, src_file, title='dbSNP', version=None):
        self.title = title
    
        super(dbSNP, self).__init__(src_file, title=self.title, version=version)
        self.vcf_reader = vcf.Reader(filename=src_file, strict_whitespace=True)
        self.output_headers = ['dbSNP_RS', 'dbSNP_Val_Status']
        self.logger = logging.getLogger(__name__)
    
    def annotate_mutation(self, mutation):
        chr, start, end = mutation.chr, int(mutation.start), int(mutation.end)
        
        #TODO: Do not annotate with a set.  Create set variables and then convert to string and annotate.  This should give a speedup.
        #TODO: This block of code is causing issue 53
        for header in self.output_headers:
            mutation.createAnnotation(header, set(), self.title)
        
        overlapping_vcf_records = []
        try:
            overlapping_vcf_records = self.vcf_reader.fetch(chr, start, end)
        except ValueError as ve:
            self.logger.warn("Exception when looking for vcf records.  Empty set of records being returned: " + repr(ve))
        
        for vcf_rec in overlapping_vcf_records:
            for header in self.output_headers:
                mutation['dbSNP_RS'].add(vcf_rec.ID)
                if 'VLD' in vcf_rec.INFO:
                    mutation['dbSNP_Val_Status'].add('byFrequency')
                if 'KGVAL' in vcf_rec.INFO or 'KGPROD' in vcf_rec.INFO or 'KGPilot1' in vcf_rec.INFO or 'KGPilot123' in vcf_rec.INFO:
                    mutation['dbSNP_Val_Status'].add('by1000genomes')    

        for header in self.output_headers:
            mutation[header] = '|'.join(mutation[header])
        
        return mutation

class dbNSFP(Datasource):
    """
    dbNSFP datasource expecting a directory of Tabix compressed and indexed tsv files.

    Example Tabix cmds:
        bgzip dbNSFP2.0_variant.chr1
        tabix -s 1 -b 2 -e 2 dbNSFP2.0_variant.chr1.gz

    This datasource requires the following annotations to be populated:
        variant_classification
        protein_change

    See the 'dbNSFP2.0.readme.txt' readme file in the datasource directory for field descriptions.
    """
    def __init__(self, src_dir, title='dbNSFP', version=None):
        self.title = title
        self.prot_regexp = re.compile('p\.([A-Z])\d+([A-Z])')
        
        super(dbNSFP, self).__init__(src_dir, title=self.title, version=version)
        self.db_obj = self._load_tabix_file_objs(src_dir)
        self.output_headers = [dbnsfp_fieldnames[14], dbnsfp_fieldnames[17]] + dbnsfp_fieldnames[21:] #only fields of interest
        self.logger = logging.getLogger(__name__)

    def annotate_mutation(self, mutation):
        if mutation['variant_classification'] == 'Missense_Mutation':
            chrn = mutation.chr
            start = mutation.start
            ref_allele, alt_allele = mutation.ref_allele, mutation.alt_allele
            regex_res = self.prot_regexp.match(mutation['protein_change'])

            if regex_res is None:
                logging.getLogger(__name__).error("Could not find a protein change for this missense mutation.")
                dbnsfp_results =[]
            else:
                ref_prot_allele, alt_prot_allele = regex_res.group(1), regex_res.group(2)
                dbnsfp_results = self._query_dbnsfp(chrn, start, start, self.db_obj)
    
            dbnsfp_data = None
            for res in dbnsfp_results:
                if all((res['ref'] == ref_allele, res['alt'] == alt_allele, res['aaref'] == ref_prot_allele, res['aaalt'] == alt_prot_allele)):
                    dbnsfp_data = res
    
            if dbnsfp_data is not None:
                for output_field in self.output_headers:
                    mutation.createAnnotation(output_field, dbnsfp_data[output_field], self.title)
            else:
                for output_field in self.output_headers:
                    mutation.createAnnotation(output_field, '', self.title)
        else:
            for output_field in self.output_headers:
                mutation.createAnnotation(output_field, '', self.title)

        return mutation

    def _load_tabix_file_objs(self, dbnsfp_dir):
        """ Returns dictionary of Tabixfile objects keyed by chromosome. """
        tabix_fnames = [f for f in os.listdir(dbnsfp_dir) if f.endswith('.gz')]
        tabix_objs = dict()
        for tabix_fname in tabix_fnames:
            contig = tabix_fname.replace('.gz', '').replace('dbNSFP2.0_variant.chr', '')
            tabix_objs[contig] = pysam.Tabixfile(os.path.join(dbnsfp_dir, tabix_fname))
    
        return tabix_objs

    def _get_dbnsfp_data_from_rows(self, dbnsfp_rows):
        return [dict(zip(dbnsfp_fieldnames, r.strip().split('\t'))) for r in dbnsfp_rows]
    
    def _query_dbnsfp(self, chromosome, start, end, tabix_objs):
        tabix_file = tabix_objs[chromosome]
        query_str = '{}:{}-{}'.format(chromosome, start, end)
        res = tabix_file.fetch(region=query_str)
        return self._get_dbnsfp_data_from_rows(res)

class Cosmic(Datasource):
    """
    Cosmic Datasource expecting a tabix compressed tsv file from Cosmic.

    Retrieves annotations based on specific genomic position and entire transcript protein sequence.

    Preferred input annotations:
    transcript_protein_position_start -- used to get COSMIC annotations by transcript protein seq
    transcript_protein_position_end -- used to get COSMIC annotations by transcript protein seq

    Output columns are:
        COSMIC_n_overlapping_mutations
            Number of COSMIC entries that overlap by genomic position.

        COSMIC_overlapping_mutation_AAs
            Protein change of overlapping COSMIC entries.  Number of entries with protein change is in parentheses.  Multiple protein
            changes are delimited with "|" character.

        COSMIC_overlapping_mutation_descriptions
            COSMIC mutation description of overlapping COSMIC entries.  Number of entries with mutation description is in parentheses.  Multiple mutation
            descriptions are delimited with "|" character.

        COSMIC_overlapping_primary_sites
            Primary site of overlapping COSMIC entries.  Number of entries with primary site is in parentheses.  Multiple primary
            sites are delimited with "|" character.
            
        
        NOTE: fusion genes are handled in a different datasource.

    """
    def __init__(self, src_file, title='COSMIC', version=None, gpp_tabix_file=None):
        self.title = title
        
        super(Cosmic, self).__init__(src_file, title=self.title, version=version)

        if gpp_tabix_file is None:
            raise ValueError("A second index by gene protein position must be specified.")

        self.db_genomePos = pysam.Tabixfile(src_file)
        header = self.db_genomePos.header.next()
        self.src_headers = header.lstrip('#').strip().split('\t')

        self.db_geneProteinPos = pysam.Tabixfile(gpp_tabix_file)
        gppHeader = self.db_geneProteinPos.header.next()
        self.gpp_headers = gppHeader.lstrip('#').strip().split('\t')

        self.output_headers = ['COSMIC_n_overlapping_mutations', 'COSMIC_overlapping_mutation_AAs',
            'COSMIC_overlapping_mutation_descriptions', 'COSMIC_overlapping_primary_sites']
        self.logger = logging.getLogger(__name__)
        
    def annotate_mutation(self, mutation):
        # transcriptPos = [mutation['tx_start'],mutation['tx_end']]
        overlapping_cosmic_entries = self.fetch_overlapping_cosmic_records(mutation.chr, mutation.start, mutation.end)

        overlappingGeneProteinPositionEntries = []
        if "transcript_protein_position_start" in mutation.keys() and "transcript_protein_position_end" in mutation.keys():
            overlappingGeneProteinPositionEntries = self.fetch_overlapping_gene_proteinPos_records(mutation['gene'],mutation["transcript_protein_position_start"], mutation["transcript_protein_position_end"])

        overlapping_cosmic_entries.extend(overlappingGeneProteinPositionEntries)

        n_overlapping_mutations = len(overlapping_cosmic_entries)
        
        mutation_AAs = collections.Counter([entry['Mutation AA'] for entry in overlapping_cosmic_entries])
        mutation_descriptions = collections.Counter([entry['Mutation Description'] for entry in overlapping_cosmic_entries])
        primary_sites = collections.Counter([entry['Primary site'] for entry in overlapping_cosmic_entries])
        
        get_output_str = lambda counter_obj: '|'.join(['%s(%d)' % s for s in sorted(counter_obj.items(), key=lambda x: x[1], reverse=True)])
        
        overlapping_mutations_AA_str = get_output_str(mutation_AAs)
        overlapping_mutation_descriptions_str = get_output_str(mutation_descriptions)
        overlapping_primary_sites_str = get_output_str(primary_sites)

        mutation.createAnnotation('COSMIC_n_overlapping_mutations', str(n_overlapping_mutations), self.title)
        mutation.createAnnotation('COSMIC_overlapping_mutation_AAs', overlapping_mutations_AA_str, self.title)
        mutation.createAnnotation('COSMIC_overlapping_mutation_descriptions', overlapping_mutation_descriptions_str, self.title)
        mutation.createAnnotation('COSMIC_overlapping_primary_sites', overlapping_primary_sites_str, self.title)
        
        return mutation

    def _fetchRecords(self, db, headers, chromosome, start, end):
        finalResult = []
        results = None
        try:
            results = self._fetch_from_tabix_file(db, chromosome, start, end)
        except ValueError as ve:
            self.logger.debug(
                "Exception when looking for COSMIC records.  Empty set of records being returned: " + repr(ve))
        if results is not None:
            for res in results:
                finalResult.append(dict(zip(headers, res.strip().split('\t'))))
        return finalResult

    def fetch_overlapping_cosmic_records(self, chromosome, start, end):
        headers = self.src_headers
        db = self.db_genomePos
        return self._fetchRecords(db, headers, chromosome, start, end)

    def fetch_overlapping_gene_proteinPos_records(self, gene, startAA, endAA):
        headers = self.gpp_headers
        db = self.db_geneProteinPos


        if (startAA.strip() == "") or (endAA.strip() == ""):
            return []
        startAA = str(int(startAA) - 1)
        return self._fetchRecords(db, headers, gene, startAA, endAA)

    def _fetch_from_tabix_file(self, tabix_file, chromosome, start, end):
        #TODO: Low priority: This only supports human, if that matters.
        if chromosome == 'X':
            chromosome = "23"
        if chromosome == 'Y':
            chromosome = "24"
        
        return tabix_file.fetch(region='%s:%s-%s' % (chromosome, str(start), str(end)))


class ReferenceDatasource(Datasource):
    """ Reference annotations.  A custom datasource initialized by genome flat files.

    Genome flat files are a simple format.  TODO: Finish this blurb.

    All annotations are with respect to the reference.

    Provides the following annotations:
        ref_context -- Small window into the reference at the variant.  The center position should be the same as Reference_Allele.  The total string should be of odd length and have a minimum length of 3.
            For example (SNV): Reference Allele is G, Chromosome is 1, Start_position and End_position are 120906037:  ref_context is CTTTTTTCGCGCAAAAATGCC  (string size is 21, in this case)


        gc_content --

        TODO: Finish the documentation.
    """

    def __init__(self, src_dir, title="Flat File Reference", version='hg19', windowSizeRef=10, windowSizeGCContent=100):
        """ Constructor
        src_dir is the parent directory of the chrXXX.txt files.
        """

        self.directoryName = src_dir
        self.windowSizeRef = windowSizeRef
        self.windowSizeGC = windowSizeGCContent
        self.logger = logging.getLogger(__name__)
        super(ReferenceDatasource, self).__init__(src_dir, title=title, version=version)
        self._filePointers = dict()
    
    def annotate_mutation(self, m):
        
        # TODO: Error checking
        iStart = int(m['start'])
        iEnd = int(m['end'])
        m.createAnnotation('ref_context', self.getRange(m['chr'], iStart - (self.windowSizeRef+1), iEnd + (self.windowSizeRef-1)), annotationSource=self.title)
        
        # Populate gc content
        gcWindow = self.getRange(m['chr'], iStart - (self.windowSizeGC+1), iEnd + (self.windowSizeGC-1))
        gcCountInWindow = gcWindow.count("C") + gcWindow.count("G") + gcWindow.count("c") + gcWindow.count("g") 
        m.createAnnotation('gc_content', "0", self.title)
        if len(gcWindow) <> 0:
            gc_content = "%0.3f" % (float(gcCountInWindow)/float(len(gcWindow)))
            m['gc_content'] = gc_content
        return m
    
    def _getFilePointer(self, fullChrFilename):
        """ Implements lazy loading of file pointers. """
        if fullChrFilename not in self._filePointers.keys():
            self._filePointers[fullChrFilename] = file(fullChrFilename, 'r')
        return self._filePointers[fullChrFilename]

    def getRange(self, chr, start, end):
        """ Returns string of the reference genome in the range requested. """
        
        chrFilename = self.convertMutationChrToFilename(chr)
        
        # TODO: We need a master lookup table.  Chip Stewart may have one for all reference builds.

        iEnd = int(end)
        iStart = max(int(start),0)
        if iEnd < iStart:
            return ""
        
        fullChrFilename = self.directoryName + "/" + chrFilename
        
        if not (os.path.exists(fullChrFilename)):
            self.logger.warn(fullChrFilename + " not found.  Please add it.")
            return ""
        else:
            fp = self._getFilePointer(fullChrFilename)
        
        fp.seek(iStart, 0) # Second parameter indicates "from the beginning of the file"
        result = fp.read(iEnd-iStart+1)
        return result

    # TODO: Need cleanup of file pointers

    def convertMutationChrToFilename(self, chr):
        """ Convert the standard mutation chromosome convention to the convention used for the filenames. 
        
        Examples (for hg19):
        GL000209.1 --> chr19_gl000209_random.txt
        X --> chrX.txt
        20 --> chr20.txt
        """
        result = "chr" + str(chr) + ".txt"
        
        # If we have a chr that starts with GL:
        #    1) Make all lowercase
        #    2) Strip away anything after '.'
        #    3) Find the file that contains the result of step 1 and 2. 
        if chr.startswith("GL"):
            tmp = chr.lower()
            if tmp.find(".") != -1:
                tmp = tmp[0:tmp.find(".")]
            
            allFiles = os.listdir(self.directoryName)
            for f in allFiles:
                if f.find(tmp)<> -1:
                    result = f 
                    break
            
        return result


class Generic_GeneProteinPositionDatasource(Generic_GenomicPosition_DataSource):
    """ For annotating protein positions and changes within a gene.
        In order for this datasource to function properly, input mutations must already
          be annotated with gene and protein_change.
          (For example, GafDatasource does this).

        Additionally, this datasource assumes that a range of positions are going to be in one annotation
            separated by '_'

    TODO: allow required annotations from the config file.
    TODO: Range is still unsupported.  53_54
    """
    def __init__(self, src_file, title='', version=None, use_binary=True, proteinPositionAnnotation="protein_change", gpColumnNames="gene,start_AA,end_AA"):
        # In this case, we want to initialize with the Datasource class
        super(Generic_GenomicPosition_DataSource, self).__init__(src_file, title=title, version=version)
        self.proteinPositionAnnotation = proteinPositionAnnotation
        self.proteinRegexp = re.compile("[A-Z\*a-z]*([0-9]+)")
        index_mode = 'gene_protein_pos'
        self.db_obj, self.output_headers = get_db_data(src_file, title, use_binary, index_mode, gpColumnNames)

    def _extractProteinPositions(self, otherTranscriptList):
        transcriptProteinPos = []
        # Parse the other_transcript annotation
        for ot in otherTranscriptList:
            # There was no protein change here (Intron, etc)
            otList = ot.split("p.")
            if len(otList) < 2:
                continue
            if len(otList) > 2:
                logging.getLogger(__name__).warn("More than one protein change detected in an other_transcript annotation.")
            proteinChange = otList[1]
            proteinPosition = self.proteinRegexp.match(proteinChange)
            if proteinPosition is not None:
                transcriptProteinPos.append(proteinPosition.group(1))

        return transcriptProteinPos

    def annotate_mutation(self, mutation):
        requiredAnnotations = ['gene', self.proteinPositionAnnotation]
        if all(field in mutation for field in requiredAnnotations):
            gene = mutation['gene']
            p = None
            transcriptProteinPos = []
            proteinChange = mutation[self.proteinPositionAnnotation].replace("p.", "")
            proteinPosition = self.proteinRegexp.match(proteinChange)
            if proteinPosition is not None:
                p = proteinPosition.group(1)

            if p is not None:
                records = get_binned_data(self.db_obj, gene, int(p), int(p))
                records = get_overlapping_records(records, int(p), int(p))

                for c in self.output_headers:
                    summarized_results = get_summary_output_string([r[c].strip() for r in records])
                    mutation.createAnnotation(c, summarized_results, annotationSource=self.title)
        else:
            missingList = []
            for r in requiredAnnotations:
                if r not in mutation:
                    missingList.append(r)
            raise MissingAnnotationException("Missing required annotation: " + str(missingList))

        for header in self.output_headers:
            if header not in mutation:
                mutation.createAnnotation(header, '', annotationSource=self.title)

        return mutation

class PositionTransformingDatasource(Datasource):
    """ Given a coordinate from system A, translates it to system B.

    Useful for wrapping mappings, such as GAF transcript protein position to uniprot protein position.

    """

    def __init__(self, src_file, title='', version=None, use_binary=True):
        raise NotImplementedError("Not implemented yet")

    def annotate_mutation(self, mutation):
        raise NotImplementedError("Not implemented yet")

class TranscriptToUniProtProteinPositionTransformingDatasource(PositionTransformingDatasource):
    """ Given a transcript protein sequence, map it to the proper place in the uniprot gene.

    This datasource requires the following annotations to be populated:
        transcript_id
        protein_change

    transcript_id is required as the key into the Shove (sqlite) database.
    protein_change is configurable (i.e. a different annotation can be used), but this is the default.

    src_file is a URL.  Such as sqlite:///db.sqlite

    Annotates gene with [self.title]_aapos

    self.title + "_" will be prepended to the given outputPositionAnnotationName

    Requires a shove database url (as src_file).  The database itself can be created using
        the scripts/uniprot_utils/createUniprotProteinSeqsAlignments.py

    If a transcript is not found in the backing database, the newposition will be ""
    """
    def __init__(self, src_file, title='', version=None, use_binary=True, inputPositionAnnotationName="protein_change", outputPositionAnnotationName="aapos"):
        super(PositionTransformingDatasource, self).__init__(src_file, title=title, version=version)
        self.db = Shove(src_file, "memory://")
        self.inputAnnotationName = inputPositionAnnotationName
        self.outputAnnotationName = self.title + "_" + outputPositionAnnotationName
        self.proteinRegexp = re.compile("[A-Z\*a-z]*([0-9]+)")

        logging.getLogger(__name__).info("Loading keys for aa xform uniprot...")
        # Since this is a readonly datasource, cache all of the available keys ahead of time.  This saves a lot of time.
        self.dbKeys = self.db.keys()
        logging.getLogger(__name__).info("Keys loaded aa xform uniprot...")

    def _parsePosition(self, position):
        """ Utility class to strip decoration from the position itself.
        """
        tmp = position
        if position.find("p.") != -1:
            tmp = position.replace("p.", "")
        proteinPosition = self.proteinRegexp.match(tmp)
        if proteinPosition is not None:
            return proteinPosition.group(1)
        else:
            return ""

    def _get_uni_pos(self,fh, AA):
        new_pos = 0
        query_AA= ''
        uni_AA = ''
        pat1 = re.compile(r'Query:\s(\d+)\s+([\w\-\*]+)\s(\d+)')
        pat2= re.compile(r'Sbjct:\s(\d+)\s+([\w\-\*]+)\s(\d+)')
        q = 0
        qline = None
        for line in fh:
            if q == 1 and line.startswith('Sbjct: '):
                sline = pat2.search(line)
                new_pos, query_AA, uni_AA = self._map_uni_pos(qline, sline, AA)
                break
            if line.startswith('Query: '):
                qline = pat1.search(line)
                #print '{0}\t{1}'.format(qline.group(1), qline.group(3))
                if int(qline.group(1)) <= AA and int(qline.group(3)) >= AA:
                    q = 1
        return new_pos, query_AA, uni_AA

    def _map_uni_pos(self,qobj, sobj, aa_pos):
        qoff = 0
        soff = 0
        for i in range(len(qobj.group(2))):
            if qobj.group(2)[i] == '-':
                qoff -= 1
            qpos = i + int(qobj.group(1)) + qoff

            if sobj.group(2)[i] == '-':
                soff -= 1
            spos = i + int(sobj.group(1)) + soff

            if qpos == aa_pos:
                return spos, qobj.group(2)[i] , sobj.group(2)[i]

    def annotate_mutation(self, mutation):
        requiredAnnotations = [self.inputAnnotationName, 'transcript_id']

        # Check that all required annotations are present.
        if not all(field in mutation for field in requiredAnnotations):
            missingList = []
            for r in requiredAnnotations:
                if r not in mutation:
                    missingList.append(r)
            raise MissingAnnotationException("Missing required annotation: " + str(missingList))

        val = mutation[self.inputAnnotationName]
        transcript_id = mutation['transcript_id']
        positionOnly = self._parsePosition(val)
        newPos = ""
        if (positionOnly is not None) and (positionOnly != "") and (transcript_id in self.dbKeys):
            adata = self.db[transcript_id]
            newPos = str(self._get_uni_pos(adata, int(positionOnly))[0])

        mutation.createAnnotation(self.outputAnnotationName, newPos, annotationSource=self.title)
        return mutation


class EnsemblTranscriptDatasource(TranscriptProvider, Datasource):
    """ Similar to a GAF datasource, but uses ensembl transcripts.
        Also, supports gencode

        Though all transcripts for GENCODE can be loaded, it is currently set to ignore any transcripts that are not
            "basic"
    """
    """This is the list of annotations that get populated by this datasource"""
    POPULATED_ANNOTATION_NAMES = set(['variant_type', 'variant_classification', 'other_transcripts', 'gene', 'gene_id', 'annotation_transcript', 'genome_change', 'strand', 'transcript_id', 'secondary_variant_classification', 'protein_change', 'codon_change', 'transcript_change', 'transcript_strand', 'gene', 'gene_type', 'gencode_transcript_tags', 'gencode_transcript_status', 'havana_transcript', 'ccds_id', 'gencode_transcript_type', 'gencode_transcript_name'])
    def __init__(self,  src_file, title='ENSEMBL', version='', tx_mode=TranscriptProvider.TX_MODE_CANONICAL, protocol="file"):
        super(EnsemblTranscriptDatasource, self).__init__(src_file=src_file, title=title, version=version)

        ensembl_index_fname = src_file + ".transcript.idx"
        ensembl_gene_to_transcript_index_fname = src_file + ".transcript_by_gene.idx"
        ensembl_genomic_position_bins_to_transcript_index_fname = src_file + ".transcript_by_gp_bin.idx"

        # Contains a key of transcript id and value of a Transcript class, with sequence data where possible.
        self.transcript_db = shove.Shove(protocol + '://%s' % ensembl_index_fname, "memory://")
        self.transcript_dbkeys = self.transcript_db.keys()

        self.gene_db = shove.Shove(protocol + '://%s' % ensembl_gene_to_transcript_index_fname, "memory://")
        self.gene_dbkeys = self.gene_db.keys()

        self.gp_bin_db = shove.Shove(protocol + '://%s' % ensembl_genomic_position_bins_to_transcript_index_fname, "memory://")
        self.gp_bin_db_dbkeys = self.gp_bin_db.keys()

        self.set_tx_mode(tx_mode)

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

    def annotate_mutation(self, mutation):
        chr = mutation.chr
        start = int(mutation.start)
        end = int(mutation.end)
        txs_unfiltered = self.get_overlapping_transcripts(chr, start, end)
        txs = self._filter_transcripts(txs_unfiltered)
        final_annotation_dict = self._create_blank_set_of_annotations()
        final_annotation_dict['variant_type'] = Annotation(value=TranscriptProviderUtils.infer_variant_type(mutation.ref_allele, mutation.alt_allele), datasourceName=self.title)

        # We have hit IGR if no transcripts come back.  Most annotations can just use the blank set.
        if len(txs) == 0:
            final_annotation_dict['variant_classification'] = self._create_basic_annotation('IGR')
            nearest_genes = self._get_nearest_genes(chr, int(start), int(end))
            final_annotation_dict['other_transcripts'] = self._create_basic_annotation(value='%s (%s upstream) : %s (%s downstream)' % (nearest_genes[0][0], nearest_genes[0][1], nearest_genes[1][0], nearest_genes[1][1]))
            final_annotation_dict['gene'] = self._create_basic_annotation('Unknown')
            final_annotation_dict['gene_id'] = self._create_basic_annotation('0')
        else:
            # Choose the best effect transcript
            chosen_tx = self._choose_transcript(txs, self.get_tx_mode(), final_annotation_dict['variant_type'].value, mutation.ref_allele, mutation.alt_allele, start, end)
            vcer = VariantClassifier()

            final_annotation_dict['annotation_transcript'] = self._create_basic_annotation(chosen_tx.get_transcript_id())
            final_annotation_dict['genome_change'] = self._create_basic_annotation(TranscriptProviderUtils.determine_genome_change(mutation.chr, mutation.start, mutation.end, mutation.ref_allele, mutation.alt_allele, final_annotation_dict['variant_type'].value))
            final_annotation_dict['strand'] = self._create_basic_annotation(chosen_tx.get_strand())

            # final_annotation_dict['transcript_exon'] = self._create_basic_annotation('')
            # final_annotation_dict['transcript_position'] = self._create_basic_annotation('')

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
        return mutation

    def _filter_transcripts(self, txs):
        """ GENCODE transcripts contain tags that are useful for QC.  If tags are present, this method will remove
        transcripts with poor QC.

        For now, just accept "basic"

        :param txs: transcripts to possibly prune
        :return: a list of same or shorter
        """
        return [tx for tx in txs if (not 'tag' in tx.get_other_attributes().keys()) or ('basic' in tx.get_other_attributes()['tag'])]

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
        lvl = tx.get_other_attributes().get('level', None)[0]
        if lvl is None:
            lvl_score = 0
        else:
            lvl_score = 4 - int(lvl)

        type_rank = 2
        type_score = 0
        if tx.get_gene_type() == "protein_coding":
            type_score = 1

        return (lvl_score << lvl_rank) + (type_score << type_rank)

    def get_overlapping_transcripts(self, chr, start, end):
        records = self._get_binned_transcripts(chr, start, end)
        return self._get_overlapping_transcript_records(records, start, end)

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
            records = self.get_overlapping_transcripts(chr, new_start, end)
            nearest_gene_border = 0
            for tx in records:
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
            records = self.get_overlapping_transcripts(chr, start, new_end)
            nearest_gene_border = int(1e9)
            for tx in records:
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

        txs -- a list of transcripts to render.
        transcriptIndicesToSkip -- a list of transcripts that are being used (i.e. not an "other transcript").  This will usually be the canonical or any transcript chosen by tx_mode.
        """
        vcer = VariantClassifier()
        other_transcripts = list()
        for i, ot in enumerate(txs):
            if i not in transcriptIndicesToSkip:
                vc = vcer.variant_classify(tx=ot, variant_type=variant_type, ref_allele=ref_allele, alt_allele=alt_allele, start=start, end=end)
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

