#!/usr/bin/env python

import oncotator
import sys

#input_file = '/data/reference/oncotator_0.6/testing/egfr_akt.compendia.maflite'
#output_file = '/data/reference/oncotator_0.6/testing/egfr_akt.0.6.maf'
#errfile = '/data/reference/oncotator_0.6/testing/egfr_akt.0.6.maflite.err'
#gaf_fname = '/data/reference/oncotator_0.6/transcript.genome.v3_0.gaf'
#gaf_transcript_sequences_fname = '/data/reference/oncotator_0.6/UCSCgene.Jan2012.v3_0.fa'
#dbsnp_fname = '/data/reference/oncotator_0.6/dbsnp_135.hg19.vcf.gz'

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print "Usage: oncotator-annotate.py <dbsnp_fname> <gaf_fname> <gaf_seqs_fname> <input_maflite_file> <output_maf_file>"
        sys.exit(1)

    dbsnp_fname = sys.argv[1]
    gaf_fname = sys.argv[2]
    gaf_transcript_sequences_fname = sys.argv[3]
    input_maflite_file = sys.argv[4]
    output_maf_file = sys.argv[5]
    errfile = output_maf_file + '.err'

    class AnnotationSet(oncotator.AnnotationSet):
        set_name = 'Default annotation set'
        gaf = oncotator.Gaf(gaf_fname=gaf_fname, gaf_transcript_sequences_fname=gaf_transcript_sequences_fname,
                            title='Gaf', version='3.0')
        dbsnp = oncotator.dbSNP(src_file = dbsnp_fname, version='dbSNP build 135')

    data = oncotator.read_maflite_file(input_maflite_file)
    data = oncotator.filter_errors_and_write_err_file(data, errfile)
    data = AnnotationSet.annotate_data(data)
    oncotator.write_maf_file(output_maf_file, data, build='hg19')
