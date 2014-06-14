# LICENSE_GOES_HERE


'''
Created on Jan 22, 2013

@author: lichtens

Simple utility script to create a fusion gene tsv from the Cosmic tsv.  

The resulting tsv can be used as a generic gene datasource.

NOTE: Use this script at your own risk.
 
'''

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import csv
from oncotator.utils.GenericTsvReader import GenericTsvReader 
from collections import OrderedDict
def parseOptions():
    
    # Process arguments
    desc = ''' '''
    epilog = ''' NOTE: This script will load large portions of the COSMIC input file into RAM.
    egrep "^[A-Z0-9]+/" CosmicMutantExportIncFus_v62_291112.tsv
    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("ds_file", type=str, help="COSMIC datasource filename. For example, 'CosmicCompleteExport_v62_261112.tsv' ")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")
    
    args = parser.parse_args()
    return args

def renderFusionGeneDictEntry(geneKey, fusionGeneDict):
    fusionGeneSubDict = fusionGeneDict[geneKey]
    resultList = []
    for k in fusionGeneSubDict.keys():
        print k
        print str(fusionGeneSubDict[k])
        summaryString = "%(k)s(%(fgene)s)" % {'k':k,'fgene':str(fusionGeneSubDict[k])}
        resultList.append(summaryString)
    return '|'.join(resultList)

if __name__ == '__main__':
    args = parseOptions()
    inputFilename = args.ds_file
    outputFilename = args.output_file
    
    tsvReader = GenericTsvReader(inputFilename)
    headers = tsvReader.getFieldNames()
    print('Found headers (input): ' + str(headers))
    if "Gene name" not in headers:
        raise NotImplementedError("Could not find Gene name column in the input file.")
    
    outputHeaders = ['gene', 'fusion_genes']
    
    # Create a dictionary where key is the gene and value is another dict: {fusion_gene:count}
    fusionGeneDict = OrderedDict()
    
    for line in tsvReader:
        fusionGene = line['Gene name']
        
        if fusionGene.find('/') == -1:
            # Not a fusion gene
            continue
        
        geneListKeys = fusionGene.split('/')
        
        for k in geneListKeys:
            if k not in fusionGeneDict.keys():
                fusionGeneDict[k] = dict()
            
            # Look for the fusion gene
            if fusionGene not in fusionGeneDict[k].keys():
                fusionGeneDict[k][fusionGene] = 0
            
            fusionGeneDict[k][fusionGene] = fusionGeneDict[k][fusionGene] + 1
        
    # Render the fusionGeneDict
    tsvWriter = csv.DictWriter(file(outputFilename,'w'), outputHeaders, delimiter='\t', lineterminator="\n")
    tsvWriter.fieldnames = outputHeaders
    tsvWriter.writeheader()
    for k in fusionGeneDict.keys():
        
        row = dict()
        row['gene'] = k
        row['fusion_genes'] = renderFusionGeneDictEntry(k, fusionGeneDict) 
        tsvWriter.writerow(row)
    
    pass
    
    