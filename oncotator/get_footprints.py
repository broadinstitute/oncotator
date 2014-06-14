# LICENSE_GOES_HERE

from argparse import ArgumentParser, RawDescriptionHelpFormatter


def parseOptions():
    epilog = '''

    This program can use a lot of RAM, as it stores all files into RAM before writing to disk.

    Then to merge the exons using bedtools, if you would like:
        mergeBed -nms -i gencode_basic_exons.bed.txt > gencode_basic_exons.merged.txt
        '''
    desc = ''' Calculate the footprint for each gene in a file.  If more than one file is specified, will output the intersection genes, with the footprints for each.'''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument('-i', '--intersection', action="store_true", help="Only include the intersection of genes in the given files.")
    parser.add_argument('-o', '--output_file', default="footprint.out", type=str, help="Output filename.  Default: footprint.out")
    parser.add_argument('bed_files', nargs='+', type=str, help="Files of genes and merged exons for determining footprints.")
    args = parser.parse_args()

    return args

def main():
    args = parseOptions()
    is_intersection = args.intersection
    bed_files = args.bed_files
    output_file = args.output_file

    gene_dicts = []
    for f in bed_files:

        # dict is gene (str) to footprint (int)
        bed_file_dict = {}
        fp = file(f, 'r')
        for line in fp:
            line_list = line.rstrip().split('\t')
            gene_set = set(line_list[3].split(";"))
            if len(gene_set) > 1:
                print("More than one gene found for exon... in " + f + "  " + str(gene_set))
            for gene in gene_set:
                start = int(line_list[1])
                end = int(line_list[2])
                try:
                    bed_file_dict[gene] += abs(end-start)
                except KeyError:
                    bed_file_dict[gene] = 0
                    bed_file_dict[gene] += abs(end-start)

        gene_dicts.append(bed_file_dict)
        fp.close()

    # Output
    output_fp = file(output_file, 'w')
    output_fp.write("gene\t" + "\t".join(bed_files) + "\n")
    # Create the keys
    gene_keys = set(gene_dicts[0].keys())
    for i in range(1,len(gene_dicts)):
        ks = set(gene_dicts[i].keys())
        if is_intersection:
            gene_keys = gene_keys.intersection(ks)
        else:
            gene_keys = gene_keys.union(ks)

    for g in gene_keys:
        output_list = [str(d[g]) for d in gene_dicts]
        output_fp.write(g + "\t" + "\t".join(output_list) + "\n")
    output_fp.close()

if __name__ == "__main__":
    main()
