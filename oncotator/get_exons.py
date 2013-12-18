from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os
from oncotator.DatasourceFactory import DatasourceFactory


def parseOptions():
    # Setup argument parser
    epilog= '''

    Error file will be [outputFilename].err.

    Transcript datasource is the full path to the directory housing the datasource.

    This script is experimental and has minimal error checking.

    To get the output of this file into a BED format, for example, with a given gene list and gencode datasource :
        get_exons gene_list.txt /my_db_dir/gencode_out2/hg19/gencode_out2.config gencode_basic_exons.txt
        awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' gencode_basic_exons.txt | sort -k1,1 -k2,2n  > gencode_basic_exons.bed.txt

    Then to merge the exons using bedtools, if you would like:
        mergeBed -nms -i gencode_basic_exons.bed.txt > gencode_basic_exons.merged.txt

    '''

    desc = ''' List all exons (non-coding can be included) in a transcript datasource. '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument('input_gene_list_file', type=str, help="A simple file with a gene name (usually Hugo Symbol) on each line")
    parser.add_argument('transcript_ds_config', type=str, help="Path to transcript datasource config file.")
    parser.add_argument('outputFilename', type=str, help="output path.  Must be writable.  Will overwrite existing files.")
    parser.add_argument('--includeNonCoding', action="store_true", help="Whether non-coding regions should be included.")
    parser.add_argument('--padding', type=str, default="0", help="Pad each exon by a number of bases on both sides.  Default: 0")
    # Process arguments
    args = parser.parse_args()

    return args


def main():
    args = parseOptions()
    input_gene_list_file = args.input_gene_list_file
    transcript_ds = args.transcript_ds_config
    outputFilename = args.outputFilename
    isNonCoding = args.includeNonCoding
    padding = args.padding
    ds = DatasourceFactory.createDatasource(transcript_ds, os.path.dirname(transcript_ds))
    input_gene_list_file_fp = file(input_gene_list_file, 'r')
    outputFileFP = file(outputFilename, 'w')
    errorFileFP = file(outputFilename + ".err", 'w')
    for line in input_gene_list_file_fp:
        gene = line.strip()
        exons = ds.retrieveExons(gene, isCodingOnly=(not isNonCoding), padding=int(padding))
        if len(exons) == 0:
            errorFileFP.write("Could not locate " + gene + "\n")
        for e in exons:
            outputFileFP.write('%s\t%s\t%s\t%s\n' % (e[0], e[1], e[2], e[3]))
    print("Done ... " + outputFilename)

if __name__ == "__main__":
    main()