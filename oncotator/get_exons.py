from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os
from oncotator.DatasourceCreator import DatasourceCreator


def parseOptions():
    # Setup argument parser
    epilog= '''

    Error file will be [outputFilename].err.

    Transcript datasource is the full path to the directory housing the datasource.  For example, dbDir/gaf/hg19/

    This script is experimental and has minimal error checking.
    '''

    desc = ''' List all exons (non-coding can be included) in a transcript datasource. '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument('input_gene_list_file', type=str, help="A simple file with a gene name (usually Hugo Symbol) on each line")
    parser.add_argument('transcript_ds', type=str, help="Path to transcript datasource config file.")
    parser.add_argument('outputFilename', type=str, help="output path.  Must be writable.  Will overwrite existing files.")
    parser.add_argument('--includeNonCoding', action="store_true", help="Whether non-coding regions should be included.")
    # Process arguments
    args = parser.parse_args()

    return args



if __name__ == "__main__":
    args = parseOptions()
    input_gene_list_file = args.input_gene_list_file
    transcript_ds = args.transcript_ds
    outputFilename = args.outputFilename
    isNonCoding = args.includeNonCoding
    ds = DatasourceCreator.createDatasource(os.path.dirname(transcript_ds), os.path.basename(transcript_ds))
    input_gene_list_file_fp = file(input_gene_list_file, 'r')
    outputFileFP = file(outputFilename, 'w')
    errorFileFP = file(outputFilename + ".err", 'w')
    for line in input_gene_list_file_fp:
        gene = line.strip()
        exons = ds.retrieveExons(gene, isCodingOnly=(not isNonCoding))
        if len(exons) == 0:
            errorFileFP.write("Could not locate " + gene + "\n")
        for e in exons:
            outputFileFP.write('%s\t%s\t%s\t%s\n' % (e[0], e[1], e[2], e[3]))
    print("Done ... " + outputFilename)