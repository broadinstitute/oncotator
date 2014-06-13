from argparse import ArgumentParser, RawDescriptionHelpFormatter
from oncotator.utils.GenericTsvReader import GenericTsvReader

__author__ = 'lichtens'

def parseOptions():
    desc = "Collapse transcript xref file, so that each transcript is only listed once and display on command line."
    epilog = ''' This script will effectively load the entire input file into RAM, so be warned.

    python createGencodeXrefTsv.py -k 0 /home/lichtens/broad_oncotator_configs/gencode/gencode.v18.metadata.RefSeq transcript_id,mRNA_id,prot_acc >~/broad_oncotator_configs/gencode_xrefseq.tsv

    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("xref_file", type=str, help="tsv file generated assumes no header.")
    parser.add_argument("column_names", type=str, help="comma separated list of column names, in order, for the tsv.  No spaces")
    parser.add_argument("-k", "--key_column", type=str, help="column number (0-index) to use as a key.  The output will have a unique value for each.")


    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parseOptions()

    xref_file = args.xref_file
    column_names = args.column_names
    key_column = int(args.key_column)

    fp = file(xref_file, 'r')
    key_dict = {}

    for line in fp:
        line_list = line.split("\t")
        key_val = line_list[key_column]
        line_list.remove(key_val)
        try:
            key_dict[key_val].append(line_list)
        except KeyError:
            key_dict[key_val] = []
            key_dict[key_val].append(line_list)

    # output
    # add header
    col_list = column_names.split(",")
    print("\t".join(col_list))
    ks = key_dict.keys()
    for k in ks:
        output_set = key_dict[k]
        if len(output_set) == 1:
            line_list = output_set[0]
            line_list.insert(key_column, k)
            line_list = [x.strip() for x in line_list]
            print("\t".join(line_list))
        else:
            final_line_sets = []
            for line_list in output_set:
                line_list.insert(key_column, k)
                if len(final_line_sets) == 0:
                    final_line_sets = [[] for i in range(0,len(line_list))]
                for i in range(0,len(line_list)):
                    if line_list[i] is None or line_list[i].strip() == "":
                        continue
                    final_line_sets[i].append(line_list[i].strip())
            final_line_str_list = ["|".join(final_line_sets[i]) for i in range(0, len(final_line_sets))]
            final_line_str_list[key_column] = k
            final_line_str = "\t".join(final_line_str_list)
            print(final_line_str)










