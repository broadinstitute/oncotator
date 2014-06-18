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

from operator import itemgetter
from oncotator.utils import Gaf
from Bio import SeqIO
import cPickle
import shove


def extractGafBins(gaf_txs):
    adjust_neg_strand_coordinates(gaf_txs)
    add_bins_to_gaf(gaf_txs)
    split_gene_field_in_gaf(gaf_txs)
    gaf_genes = create_Genes_dict(gaf_txs)
    add_bins_to_Genes(gaf_genes)
    infer_canonical_transcript_for_genes(gaf_genes, gaf_txs)
    gaf_txs_binned = create_gaf_dict(gaf_txs)
    gaf_genes_binned = create_binned_Genes_dict(gaf_genes)
    fix_mitochondrial_chrs_gaf_binned(gaf_txs_binned)
    remove_chrstr_from_chr_values(gaf_genes_binned)
    remove_redundant_alt_loc_from_binned_genes(gaf_genes_binned)
    return gaf_genes_binned, gaf_txs_binned


def create_gaf_dicts(gaf_fname):
    gaf_fh = open(gaf_fname)
    gaf_txs = parse_gaf(gaf_fh)
    gaf_genes_binned, gaf_txs_binned = extractGafBins(gaf_txs)
    return [gaf_txs_binned, gaf_genes_binned]


def index_gaf(gaf_fname, output_fname):
    cPickle.dump(create_gaf_dicts(gaf_fname), open(output_fname, 'wb'))


def index_gaf_fastas(gaf_transcript_seqs_fname, output_fname, protocol="sqlite"):
    fh_transcripts = SeqIO.parse(gaf_transcript_seqs_fname, 'fasta')
    transcripts_shv = shove.Shove(protocol + ':///%s' % output_fname)

    j = 0
    for transcript in fh_transcripts:
        if j % 1000 == 0:
            print j
        j += 1
        raw_seq = str(transcript.seq)
        transcripts_shv[transcript.name] = raw_seq

    transcripts_shv.close()


def parse_gaf(gaf_file):
    '''
    Returns a list of dicts with the following keys:

    '''
    ##Filter usable transcripts from GAF
    #
    gaf_transcripts = list()
    i = 0

    n_total_transcripts = 0
    n_transcripts_with_invalid_CDS = 0
    n_transcripts_with_no_gene_symbol = 0
    n_useable_transcripts = 0
    n_useable_miRNAs = 0

    print "parsing gaf file..."
    for g in gaf_file:
        g = Gaf.Gaf(g)

        if i % 1000000 == 0 and i != 0: print i
        i += 1
        if g.featureType == 'transcript' and g.compositeType == 'genome':
            n_total_transcripts += 1

            ####remove transcripts that have a CDSstart but no CDSstop
            record = dict()
            #transcript_coords replaces transcript_start and transcript_end lists
            record['transcript_coords'] = [(int(c.split('-')[0]), int(c.split('-')[1])) for c in
                                           g.featureCoordinates.split(',')]
            record['chr'], genomic_coords, record['strand'] = g.compositeCoordinates.split(':')
            record['chr'] = record['chr'].replace('chr', '')
            #genomic_coords replaces genomic_starts and genomic_ends
            record['genomic_coords'] = [(int(c.split('-')[0]), int(c.split('-')[1]),) for c in
                                        genomic_coords.split(',')]
            record['tx_start'], record['tx_end'] = record['genomic_coords'][0][0], record['genomic_coords'][-1][1]
            record['footprint_start'], record['footprint_end'] = record['tx_start'], record['tx_end']
            if record['strand'] == '-':
                record['footprint_start'], record['footprint_end'] = record['footprint_end'], record['footprint_start']
            record['cds_start'], record['cds_stop'] = None, None
            if g.featureInfo != '':
                if g.featureInfo[0] == ';':
                    n_transcripts_with_invalid_CDS += 1
                    continue

                feature_info_lst = g.featureInfo.split(';')
                if len(feature_info_lst) == 1:
                    record['confidence'] = g.featureInfo.split('=')[1]
                elif len(feature_info_lst) == 3:
                    record['confidence'], record['cds_start'], record['cds_stop'] = [int(c.split('=')[1]) for c in
                                                                                     g.featureInfo.split(';')]
                    record['code_len'] = record['cds_stop'] - record['cds_start'] + 1
                else:
                    raise Exception("unexpected number of fields in fieldinfo: %s" % g.featureInfo)

            record['tx_len'] = record['transcript_coords'][-1][1] - record['transcript_coords'][0][0] + 1
            record['n_exons'] = len(record['transcript_coords'])
            record['transcript_id'] = g.featureId
            record['gene'] = g.gene
            record['transcript_coords'] = tuple(record['transcript_coords'])
            record['genomic_coords'] = tuple(record['genomic_coords'])
            record['class'] = 'mRNA'
            record['GeneLocus'] = g.geneLocus
            record['FeatureAlias'] = g.featureAliases
            record['EntryNumber'] = g.entryNumber
            n_useable_transcripts += 1
            gaf_transcripts.append(record)

        elif g.featureType == 'pre-miRNA' and g.compositeCoordinates != 'UNKNOWN':
            record = dict()
            record['transcript_coords'] = [(int(c.split('-')[0]), int(c.split('-')[1])) for c in
                                           g.featureCoordinates.split(',')]
            record['chr'], genomic_coords, record['strand'] = g.compositeCoordinates.split(':')
            record['chr'] = record['chr'].replace('chr', '')
            record['genomic_coords'] = [(int(c.split('-')[0]), int(c.split('-')[1]),) for c in
                                        genomic_coords.split(',')]
            record['transcript_id'] = g.featureId
            record['gene'] = g.gene
            record['transcript_coords'] = tuple(record['transcript_coords'])
            record['genomic_coords'] = tuple(record['genomic_coords'])
            record['n_exons'] = len(record['transcript_coords'])
            record['tx_start'], record['tx_end'] = record['genomic_coords'][0][0], record['genomic_coords'][-1][1]
            record['footprint_start'], record['footprint_end'] = record['tx_start'], record['tx_end']
            if record['strand'] == '-':
                record['footprint_start'], record['footprint_end'] = record['footprint_end'], record['footprint_start']
            record['class'] = 'miRNA'
            record['GeneLocus'] = g.geneLocus
            record['FeatureAlias'] = g.featureAliases
            record['EntryNumber'] = g.entryNumber
            record['confidence'] = 0 # use 0 for miRNAs
            n_useable_miRNAs += 1
            gaf_transcripts.append(record)

    n_total_in_gaf = i
    print "Total lines in GAF: %d" % (n_total_in_gaf)
    print "Total transcripts in GAF: %d" % (n_total_transcripts)
    print "Total useable transcripts in GAF: %d" % (n_useable_transcripts)
    print "Total useable miRNAs in GAF: %d" % (n_useable_miRNAs)
    print "Total unuseable transcripts due to lack of gene: %d" % (n_transcripts_with_no_gene_symbol)
    print "Total unuseable transcripts due to invalid CDS coordinates: %d" % (n_transcripts_with_invalid_CDS)
    return gaf_transcripts


def adjust_neg_strand_coordinates(data):
    for d in data:
        if d['strand'] == '-':
            d['footprint_start'], d['footprint_end'] = d['footprint_end'], d['footprint_start']
            d['genomic_coords'] = tuple((c[1], c[0]) for c in d['genomic_coords'][::-1])
            d['tx_start'], d['tx_end'] = d['tx_end'], d['tx_start']


def add_bins_to_gaf(data):
    for g in data:
        g['bin'] = region2bin(g['footprint_start'], g['footprint_end'])


def split_gene_field_in_gaf(data):
    for g in data:
        g['gaf_gene_id'] = g['gene'] #keep unsplit gene id from GAF
        if g['gaf_gene_id'] == '': ##usually for some miRNAs
            g['gaf_gene_id'] = g['transcript_id']
            g['gene'] = g['transcript_id']
            g['id'] = ''
        else:
            gene_split_lst = g['gene'].split('|')
            if len(gene_split_lst) == 3:
                g['gene'], g['id'], g['MofN'] = gene_split_lst
            elif len(gene_split_lst) == 2:
                g['gene'], g['id'] = gene_split_lst


def create_Genes_dict(gaf):
    genes = dict()
    for t in gaf:
        gaf_gene_id = t['gaf_gene_id']
        try:
            gene_name, id = t['gene'], t['id']
        except KeyError:
            return t

        if gaf_gene_id not in genes:
            genes[gaf_gene_id] = dict()
            genes[gaf_gene_id]['transcripts'] = [t['transcript_id']]
            genes[gaf_gene_id]['class'] = t['class']
            genes[gaf_gene_id]['gene'] = gene_name
            if 'MofN' in t:
                genes[gaf_gene_id]['MofN'] = t['MofN']
            genes[gaf_gene_id]['id'] = id
            genes[gaf_gene_id]['loci'] = list()
            genes[gaf_gene_id]['loci'].extend(t['GeneLocus'].split(';'))
            genes[gaf_gene_id]['confidence'] = t['confidence']
        else:
            genes[gaf_gene_id]['transcripts'].append(t['transcript_id'])
            genes[gaf_gene_id]['loci'].extend(t['GeneLocus'].split(';'))
            if int(t['confidence']) > int(genes[gaf_gene_id]['confidence']):
                genes[gaf_gene_id]['confidence'] = t['confidence']
        del t['GeneLocus']
    for gene in genes.values():
        ##loop through loci and set coordinate values
        alt_loci = list()
        gene['loci'] = sorted(gene['loci'])
        try:
            gene['chr'], temp_coords, gene['strand'] = gene['loci'][0].split(':')
            alt_loci_start_idx = 1
        except ValueError:
            #if first loci is broken, try from second
            gene['chr'], temp_coords, gene['strand'] = gene['loci'][1].split(':')
            alt_loci_start_idx = 2
        gene['start'], gene['end'] = temp_coords.split('-')
        gene['start'], gene['end'] = int(gene['start']), int(gene['end'])
        for locus in gene['loci'][alt_loci_start_idx:]:
            temp_dict = dict()
            try:
                temp_dict['chr'], temp_coords, temp_dict['strand'] = locus.split(':')
            except ValueError:
                continue #skip loci that cannot be parsed (e.g. "chrX:4")
            temp_dict['start'], temp_dict['end'] = temp_coords.split('-')
            temp_dict['start'], temp_dict['end'] = int(temp_dict['start']), int(temp_dict['end'])
            alt_loci.append(temp_dict)
        del gene['loci']
        gene['alt_loci'] = alt_loci
    return genes


def add_bins_to_Genes(genes):
    for g in genes:
        genes[g]['bin'] = region2bin(genes[g]['start'], genes[g]['end'])


def infer_canonical_transcript_for_genes(genes, gaf):
    """
    Canonical transcript for each gene is determined by take top scoring transcript by confidence
    score and using coding length as a tie-breaker.
    """
    idx = dict((t['transcript_id'], i) for i, t in enumerate(gaf))
    for g in genes:
        txs = [gaf[idx[t]] for t in genes[g]['transcripts']]

        canonical_tx = ''
        cds_len = 0
        confidence_score = None
        for t in txs:
            if not canonical_tx:
                canonical_tx = t['transcript_id']
                cds_len = t.get('code_len', 0)
                confidence_score = int(t['confidence'])
                continue
            if int(t['confidence']) > confidence_score:
                canonical_tx = t['transcript_id']
                cds_len = t.get('code_len', 0)
                confidence_score = int(t['confidence'])
            elif int(t['confidence']) == confidence_score and t.get('code_len', 0) > cds_len:
                canonical_tx = t['transcript_id']
                cds_len = t.get('code_len', 0)

        genes[g]['canonical_tx'] = canonical_tx


def create_gaf_dict(gaf):
    gaf = sorted(gaf, key=itemgetter('footprint_end'))
    gaf = sorted(gaf, key=itemgetter('footprint_start'))
    gaf = sorted(gaf, key=itemgetter('chr'))
    gaf_dict = dict()
    for g in gaf:
        if g['chr'] in gaf_dict:
            if g['bin'] in gaf_dict[g['chr']]:
                gaf_dict[g['chr']][g['bin']].append(g)
            else:
                gaf_dict[g['chr']][g['bin']] = [g]
        else:
            gaf_dict[g['chr']] = {g['bin']: [g]}
    return gaf_dict


def create_binned_Genes_dict(genes):
    gene_lst = [genes[g] for g in genes]
    gene_lst = sorted(gene_lst, key=itemgetter('end'))
    gene_lst = sorted(gene_lst, key=itemgetter('start'))
    gene_lst = sorted(gene_lst, key=itemgetter('chr'))
    genes_binned = dict()
    for g in gene_lst:
        if g['chr'] in genes_binned:
            if g['bin'] in genes_binned[g['chr']]:
                genes_binned[g['chr']][g['bin']].append(g)
            else:
                genes_binned[g['chr']][g['bin']] = [g]
        else:
            genes_binned[g['chr']] = {g['bin']: [g]}
    return genes_binned


def fix_mitochondrial_chrs_gaf_binned(gaf_binned):
    if 'M' not in gaf_binned:
        gaf_binned['M'] = dict()

    if 'M_rCRS' in gaf_binned:
        for bin in gaf_binned['M_rCRS']:
            if bin in gaf_binned['M']:
                gaf_binned['M'][bin].extend(gaf_binned['M_rCRS'][bin])
            else:
                gaf_binned['M'][bin] = gaf_binned['M_rCRS'][bin]
        del gaf_binned['M_rCRS']

    if 'MT' in gaf_binned:
        for bin in gaf_binned['MT']:
            if bin in gaf_binned['M']:
                gaf_binned['M'][bin].extend(gaf_binned['MT'][bin])
            else:
                gaf_binned['M'][bin] = gaf_binned['MT'][bin]
        del gaf_binned['MT']


def remove_chrstr_from_chr_values(binned_genes):
    for chr_key in binned_genes:
        for bin in binned_genes[chr_key]:
            for genedict in binned_genes[chr_key][bin]:
                genedict['chr'] = genedict['chr'].replace('chr', '')
                if 'alt_loci' in genedict:
                    for alt_locus in genedict['alt_loci']:
                        alt_locus['chr'] = alt_locus['chr'].replace('chr', '')
        if chr_key.startswith('chr'):
            new_key = chr_key.replace('chr', '')
            binned_genes[new_key] = binned_genes[chr_key]
            del binned_genes[chr_key]


def remove_redundant_alt_loc_from_binned_genes(binned_genes):
    for chr in binned_genes:
        for bin in binned_genes[chr]:
            for genedict in binned_genes[chr][bin]:
                if not genedict['alt_loci']:
                    continue
                to_delete = list()
                for i, locus in enumerate(genedict['alt_loci']):
                    if locus['chr'] == genedict['chr'] and locus['start'] == genedict['start'] and locus['end'] == \
                            genedict['end']:
                        to_delete.append(i)
                for i in to_delete[::-1]:
                    genedict['alt_loci'].pop(i)


def region2bin(beg, end):
    end = end - 1
    if beg >> 17 == end >> 17: return ((1 << 12) - 1) / 7 + (beg >> 17)
    if beg >> 20 == end >> 20: return ((1 << 9) - 1) / 7 + (beg >> 20)
    if beg >> 23 == end >> 23: return ((1 << 6) - 1) / 7 + (beg >> 23)
    if beg >> 26 == end >> 26: return ((1 << 3) - 1) / 7 + (beg >> 26)
    return 0


def region2bins(beg, end):
    bins = [0]
    bins.extend(range(1 + (beg >> 26), 1 + (end >> 26) + 1))
    bins.extend(range(9 + (beg >> 23), 9 + (end >> 23) + 1))
    bins.extend(range(73 + (beg >> 20), 73 + (end >> 20) + 1))
    bins.extend(range(585 + (beg >> 17), 585 + (end >> 17) + 1))
    return bins