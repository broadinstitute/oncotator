from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os
from os.path import expanduser
from oncotator.DatasourceFactory import DatasourceFactory

__author__ = 'lichtens'

def parseOptions():
    # Process arguments
    desc = '''Output each transcript ID with its protein sequence (tsv) to the specified file.'''
    epilog = '''

    Currently, only works for GENCODE/ENSEMBL datasources.

    NOTE:  Filter option is ignored in the datasource, so all transcripts are rendered.

    This script is not supported and is not guaranteed to function properly at any point.

    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("gencode_ds_loc", type=str, help="Location of the GENCODE datasource config file -- E.g. /bulk/dbDir/gencode_ds/hg19/gencode_ds.config")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")


    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parseOptions()
    output_file = expanduser(args.output_file)
    gencode_ds_loc = expanduser(args.gencode_ds_loc)

    # Instantiate a gencode datasource
    gencode_ds = DatasourceFactory.createDatasource(configFilename=gencode_ds_loc, leafDir=os.path.dirname(gencode_ds_loc))

    # Get all transcript IDs in the datasource
    txs = gencode_ds.getTranscriptDict()
    tx_ids = txs.keys()
    num_tx_ids = len(tx_ids)

    # Initialize output file
    fp = file(output_file, 'w')
    ctr = 0
    for tx_id in tx_ids:
        ctr += 1
        if (ctr % 2000) == 0:
            print(str(ctr) + "/" + str(num_tx_ids))
            fp.flush()

        tx_protein_seq = txs[tx_id].get_protein_seq()
        fp.write("%s\t%s\n" % (tx_id, tx_protein_seq))

    fp.close()

