#!/usr/bin/env python
# LICENSE_GOES_HERE

import sys
from oncotator.index.TabixIndexer import TabixIndexer
from oncotator.index.gaf import index_gaf, index_gaf_fastas
import pysam

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Usage: oncotator-index <gaf|gaf-seqs|cosmic> <input_file>"
        sys.exit(1)

    datasource_type, input_fname = sys.argv[1:]
    output_fname = input_fname + '.idx'
    if datasource_type == 'cosmic':
        output_fname = input_fname + '.tbi'

    if datasource_type not in ['gaf', 'gaf-seqs', 'cosmic']:
        print "%s is not a valid datasource type." % datasource_type
        sys.exit(1)

    print "Input datasource: %s" % input_fname
    print "Output indexed file: %s" % output_fname

    if datasource_type == 'gaf':
        index_gaf(input_fname, output_fname)
    elif datasource_type == 'gaf-seqs':
        index_gaf_fastas(input_fname, output_fname, protocol="file")
    elif datasource_type == 'cosmic':

        ###This code needs to go in a separate module
        # TODO: Leverage TabixIndexer class
        in_fh = open(input_fname)
        line = in_fh.next()
        line = in_fh.next() if line == '\n' else line
        headers = line.strip('\n').split('\t')
        parse_cosmic_line = lambda input_line: dict(zip(headers, input_line.strip('\n').split('\t')))

        source = 'cosmic_v62' #should not be hardcoded
        feature = 'variation'
        score = '.'
        strand = '.'
        phase = '.'

        new_cosmic_lines = list()
        for line in in_fh:
            if line.strip() == '':
                continue

            tsv_data = parse_cosmic_line(line)
            try:
                chromosome, coords = tsv_data['Mutation GRCh37 genome position'].split(':')
            except ValueError:
                #skip cosmic entries with no position data
                continue
            else:
                start, end = coords.split('-')
                tsv_data['chromosome'] = chromosome
                tsv_data['start'] = start
                tsv_data['end'] = end

                new_cosmic_lines.append(tsv_data)

        #tabix needs file to be sorted
        new_cosmic_lines.sort(key=lambda x: int(x['end']))
        new_cosmic_lines.sort(key=lambda x: int(x['start']))
        new_cosmic_lines.sort(key=lambda x: x['chromosome'])

        headers.extend(['chromosome', 'start', 'end'])
        out = open(output_fname, 'w')
        out.write('#' + '\t'.join(headers) + '\n')
        for tsv_data in new_cosmic_lines:
            out.write('\t'.join([tsv_data[h] for h in headers]) + '\n')
        out.close()

        pysam.tabix_index(filename=output_fname, seq_col=15, start_col=16, end_col=17)

        print("Creating second index for AA position...")
        # Create a second index file by gene and start AA, endAA
        TabixIndexer.indexGeneProteinPosition("Gene name", "Mutation AA", input_fname, output_fname + ".byAA")