#####################################################################################################################
#
# This script will parse a VCF file and create three files in
# oncotator input format: one for SNPs, one for insertions and one for
# deletions. https://www.broadinstitute.org/oncotator/help/#inputformat
# It will use the name of the VCF input as prefix and
# output files are created at the current working directory.
#
#####################################################################################################################

import vcf
import argparse
import os.path
import csv


parser = argparse.ArgumentParser(description='Convert VCF to oncotator input format.')
parser.add_argument('vcf', type=argparse.FileType('r'))
args   = parser.parse_args()

prefix = os.path.basename(args.vcf.name).replace('.vcf','')
vcfr   = vcf.Reader( args.vcf )


with open("%s_snps" % prefix, 'w') as snps, open("%s_ins" % prefix, 'w') as ins, open("%s_dels" % prefix, 'w') as dels:

    snp_wrtr = csv.writer(snps, delimiter=' ')
    ins_wrtr = csv.writer(ins,  delimiter=' ')
    del_wrtr = csv.writer(dels, delimiter=' ')
    
    for v in vcfr:
        if len(v.FILTER)==0:
            for alt in v.ALT:
                if v.is_snp:
                    snp_wrtr.writerow([v.CHROM, v.POS, v.POS, v.REF, alt])
                elif not v.is_deletion and v.is_indel: # insertions
                    ins_wrtr.writerow([v.CHROM, v.POS, v.POS + 1, '-', alt])                
                elif v.is_deletion:
                    del_wrtr.writerow([v.CHROM, v.POS, v.POS+(len(v.REF)-1), v.REF, '-']) 
