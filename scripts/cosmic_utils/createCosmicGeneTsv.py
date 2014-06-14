# LICENSE_GOES_HERE


from argparse import ArgumentParser, RawDescriptionHelpFormatter
import csv
from oncotator.utils.GenericTsvReader import GenericTsvReader


def parseOptions():

    # Process arguments
    desc = '''Create gene table that handles the total_alterations_in_gene and tissue_types_affected'''
    epilog = '''
    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("ds_file", type=str, help="COSMIC datasource filename. For example, 'CosmicCompleteExport_v62_261112.tsv' ")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parseOptions()
    inputFilename = args.ds_file
    outputFilename = args.output_file

    outputHeaders = ['gene', 'total_alterations_in_gene', 'tissue_types_affected']

    tsvReader = GenericTsvReader(inputFilename)
    headers = tsvReader.getFieldNames()
    print('Found headers (input): ' + str(headers))
    if "Gene name" not in headers:
        raise NotImplementedError("Could not find Gene name column in the input file.")

    if 'Primary site' not in headers:
        raise NotImplementedError("Could not find Primary site column in the input file.")

    # Construct dictionary that is [gene][histology/tissue type] = count, where count is the total for that histology
    #   and that gene
    geneDictionary = dict()
    for line in tsvReader:
        gene = line['Gene name']
        # Skip blank genes
        if gene is None or gene.strip() == "":
            continue

        # Skip fusion genes
        if gene.find("/") != -1:
            continue

        if gene not in geneDictionary.keys():
            geneDictionary[gene] = dict()

        site = line['Primary site']
        if site not in geneDictionary[gene].keys():
            geneDictionary[gene][site] = 0
        geneDictionary[gene][site] += 1

    # Write tsv output file.
    tsvWriter = csv.DictWriter(file(outputFilename,'w'), outputHeaders, delimiter='\t', lineterminator="\n")
    tsvWriter.fieldnames = outputHeaders
    tsvWriter.writeheader()
    sortedGenes = sorted(geneDictionary.keys())
    for g in sortedGenes:
        row = dict()
        # Generate
        row['gene'] = g

        tissues = []
        total = 0
        for h in sorted(geneDictionary[g].keys()):
            tissues.append(h + "(" + str(geneDictionary[g][h]) + ")")
            total += geneDictionary[g][h]
        row['total_alterations_in_gene'] = str(total)
        row['tissue_types_affected'] = "|".join(tissues)
        tsvWriter.writerow(row)
    print("Done")
