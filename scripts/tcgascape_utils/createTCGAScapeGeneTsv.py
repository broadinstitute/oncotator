# LICENSE_GOES_HERE


'''
Created on Jan 23, 2013

@author: lichtens

This script will take a TCGAScape directory and create a generic gene TSV file.
'''
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import os
import csv

def parseOptions():
    
    # Process arguments
    desc = ''' '''
    epilog = ''' 
    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("ds_dir", type=str, help="Directory that contains the gene text files. ")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")
    
    args = parser.parse_args()
    return args

def fetch_tumorscape_columns(data, dbs, title):
    for m in data:
        amp_peaks = []
        del_peaks = []
        tumorscape_dir = dbs['tumorscape_dir']
        tumorscape_data_file = os.path.join(tumorscape_dir, ''.join([m['Hugo_Symbol'], '.txt']))
        if os.path.isfile(tumorscape_data_file):
            lines = open(tumorscape_data_file,'r').readlines()
            lines = [line.strip() for line in lines]
            amp_peaks_section = False
            del_peaks_section = False
            for l in  lines:
                if l == '':
                    amp_peaks_section, del_peaks_section = False, False
                elif l == 'Frequency Amplified':
                    amp_peaks_section, del_peaks_section = True, False
                elif l == 'Frequency Deleted':
                    amp_peaks_section, del_peaks_section = False, True
                elif amp_peaks_section:
                    l = l.split('\t')
                    if l[2] == 'Yes' and float(l[5]) < 0.25:
                        amp_peaks.append('%s(%s;%s)' % (l[0], l[4].strip(), l[5].strip()))
                elif del_peaks_section:
                    l = l.split('\t')
                    if l[2] == 'Yes' and float(l[5]) < 0.25:
                        del_peaks.append('%s(%s;%s)' % (l[0], l[4].strip(), l[5].strip()))
        m['{0}_Amplification_Peaks'.format(title)] = '|'.join(amp_peaks)
        m['{0}_Deletion_Peaks'.format(title)] = '|'.join(del_peaks)
        yield m

def determineGenes(tcgascape_dir):
    ''' Get all of the available text files, which are already named by the Hugo Symbol'''
    dirContents = os.listdir(tcgascape_dir)
    result = []
    for f in dirContents:
        if f.endswith(".txt"):
            geneName = f.replace(".txt", "")
            result.append(geneName)
    
    return result

def add_TCGAscape_to_genes(tcgascape_dir):
    genes = determineGenes(tcgascape_dir)
    result = dict()
    data = fetch_tumorscape_columns([{'Hugo_Symbol':g} for g in genes], {'tumorscape_dir': tcgascape_dir}, 'TCGAscape')
    data = list(data)
    for d in data:
        result[d['Hugo_Symbol']] = dict()
        if any((d.get('TCGAscape_Amplification_Peaks',''), d.get('TCGAscape_Deletion_Peaks',''))):
            result[d['Hugo_Symbol']]['TCGAscape'] = dict()
            for t in ['TCGAscape_Amplification_Peaks', 'TCGAscape_Deletion_Peaks']:
                if d.get(t,''):
                    result[d['Hugo_Symbol']][t[10:]] = d[t]
    return result

def main():
    args = parseOptions() 
    
    data = add_TCGAscape_to_genes(args.ds_dir)
    
    # TODO: Write out a tsv
    fp = file(args.output_file, 'w')
    tsvWriter = csv.DictWriter(fp, ['gene','Amplification_Peaks', 'Deletion_Peaks'], delimiter='\t', lineterminator="\n")
    tsvWriter.writeheader()
    
    for k in data.keys():
        geneDict = data[k]
        
        row = dict()
        row['gene'] = k
        row['Amplification_Peaks'] = geneDict.get('Amplification_Peaks','')
        row['Deletion_Peaks'] = geneDict.get('Deletion_Peaks','')
         
        tsvWriter.writerow(row)
    
    
if __name__ == '__main__':
    main()