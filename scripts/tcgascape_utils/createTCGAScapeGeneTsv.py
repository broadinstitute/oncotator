"""
By downloading the PROGRAM you agree to the following terms of use:

BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY

This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").

WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and

WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.

NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:

1. DEFINITIONS
1.1 "PROGRAM" shall mean the object code and source code known as Oncotator 1.0 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.  BROAD acknowledges that the PROGRAM employs one or more public domain code(s) that are freely available for public use.

2. LICENSE
2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.  LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.

3. OWNERSHIP OF INTELLECTUAL PROPERTY
LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.

Copyright 2014 Broad Institute, Inc.
Notice of attribution:  The Oncotator 1.0 program was made available through the generosity of the Broad Institute, Inc.

LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.

4. INDEMNIFICATION
LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorney fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.

5. NO REPRESENTATIONS OR WARRANTIES
THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.

6. ASSIGNMENT
This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.

7. MISCELLANEOUS
7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
"""


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