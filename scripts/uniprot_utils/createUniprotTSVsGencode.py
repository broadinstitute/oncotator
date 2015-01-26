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
from os.path import expanduser
from oncotator.DatasourceFactory import DatasourceFactory


'''
Created on Jan 28, 2013

@author: lichtens

This script will process UniProt data into files that can be easily loaded into 
    Oncotator.  This script should not be used by end users and is meant for
    developers or those who maintain datasources.
    
    This is a preprocess script.

NOTE: This script takes advantage of the fact that duplicate records contain 
    the same information in all relevant fields.

This script creates the simple_uniprot datasource.
'''
import collections
from Bio import SwissProt
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import os
import csv
import cPickle
from Bio import SeqIO

def parseOptions():
    
    # Process arguments
    desc = ''' Create one tsv files with basic uniprot data.
    Creates the simple uniprot tsv indexed by gene (HUGO symbol).
    This script requires 3GB RAM.'''
    epilog = ''' 
    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("swiss_file", type=str, help="SwissProt file. ")
    parser.add_argument("trembl_file", type=str, help="TREMBL file. ")
    parser.add_argument("gencode_ds", type=str, help="GENCODE datasource config file. ")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")
    
    args = parser.parse_args()
    return args


def get_feature_type(feature):
    if feature[0] in ['INIT_MET', 'SIGNAL', 'TRANSIT', 'PROPEP', 'CHAIN', 'PEPTIDE']:
        return 'molecule_processing'
    elif feature[0] in ['TOPO_DOM', 'TRANSMEM', 'DOMAIN', 'REPEAT', 'CA_BIND', 'ZN_FING', 'DNA_BIND', 'NP_BIND', 'REGION', 'COILED', 'MOTIF', 'COMPBIAS']:
        return 'regions'
    elif feature[0] in ['ACT_SITE', 'METAL', 'BINDING', 'SITE']:
        return 'sites'
    elif feature[0] in ['NON_STD', 'MOD_RES', 'LIPID', 'CARBOHYD', 'DISULFID', 'CROSSLNK']:
        return 'amino_acid_modifications'
    elif feature[0] in ['VAR_SEQ', 'VARIANT']:
        return 'natural_variations'
    elif feature[0] in ['MUTAGEN', 'UNSURE', 'CONFLICT', 'NON_CONS', 'NON_TER']:
        return 'experimental_info'
    elif feature[0] in ['HELIX', 'TURN', 'STRAND']:
        return 'secondary_structure'

def parse_uniprot_data(data):
    return dict((s.entry_name, s) for s in SwissProt.parse(file(uniprot_swiss_fname, 'r')))
    
def get_uniprot_features_table(record):
    features = record.features
    features_dict = collections.defaultdict(list)
    
    for i,f in enumerate(features):  ##Removes FTIds from feature table
        if f[3].find('. /') != -1:
            features[i] = (f[0],f[1],f[2],f[3][:(f[3].find('. /')+1)],f[4])

    for feature in features:
        features_dict[get_feature_type(feature)].append(feature)

    for k in features_dict.keys():  ##Remove nulls
        if '' in features_dict[k]:
            features_dict[k].discard('')
            
    features_dict = dict(features_dict)
    
    if 'natural_variations' in features_dict:   ##Removes dbSNPs and isoform differences
        to_remove = []
        for i,v in enumerate(features_dict['natural_variations']):
            if v[3].find('isoform') != -1 or v[3].find('dbSNP') != -1:
                to_remove.append(i)
        for t in to_remove[::-1]:
            del(features_dict['natural_variations'][t])
        if not features_dict['natural_variations']:
            del(features_dict['natural_variations'])
            
    return features_dict

def get_GO_terms(record):
    goterms = collections.defaultdict(list)
    change = {'P':'biological_process', 'C':'cellular_component', 'F':'molecular_function'}
    cross_references = record.cross_references
    for ref in cross_references:
        if ref[0] == 'GO':
            goterms[change[ref[2][0]]].append((ref[1],ref[2][2:]))
    
    goterms = dict(goterms)
    #for k in goterms.keys():
    #    goterms[k] = '|'.join(goterms[k])
    
    return goterms
    
def get_drugbank_data(record):
    drug_bank_data = []
    for ref in record.cross_references:
        if ref[0] == 'DrugBank':
            drugBankString = ref[2] + "(" + ref[1]+ ")"
            drug_bank_data.append(drugBankString)
    
    return drug_bank_data

def add_uniprot_data_to_Genes(genes, swiss_data, trembl_data):
    gene_list = genes.keys()
    num_genes = gene_list
    for i,gene_symbol in enumerate(gene_list):
        g = genes[gene_symbol]
        try:
            entry_name = g['uniprot_entry_name']
        except KeyError:
            continue
        try:
            uniprot_record = swiss_data[entry_name]
        except KeyError:
            uniprot_record = trembl_data[entry_name]
        go_terms = get_GO_terms(uniprot_record)
        drugbank_list = get_drugbank_data(uniprot_record)
        drugbank_list.sort()
        feature_table = get_uniprot_features_table(uniprot_record)
        if go_terms: g['GO'] = go_terms
        if drugbank_list: g['DrugBank'] = drugbank_list
        if feature_table:
            g['uniprot_features'] = feature_table
            if isinstance(feature_table,str):
                print("Uniprot features is a string: " + g['uniprot_features'])

    return genes


def create_gene_to_uniprot_dict(data, data_id, gene, gene_ids_2_entrynames):
    g = {}
    entry_names = list(gene_ids_2_entrynames[data_id][gene])
    g['uniprot_entry_name'] = entry_names.pop(0)
    if entry_names:
        g['alt_uniprot_entry_names'] = entry_names
    accessions = data[g['uniprot_entry_name']].accessions
    g['uniprot_accession'] = accessions[0]
    if len(accessions) > 1:
        g['alt_uniprot_accessions'] = accessions[1:]
    return g


def add_uniprot_ids_to_Genes(genes, swiss_data, trembl_data):
    """

    :param genes: gene names used in the transcript datasource
    :param swiss_data:
    :param trembl_data:
    :return:
    """

    ### Map GeneIDs to UniProt entry_names
    gene_ids_2_entrynames = {'swiss': collections.defaultdict(list), 'trembl': collections.defaultdict(list)}
    for db in gene_ids_2_entrynames:
        for s in eval(db + '_data').values():

            # Parse out the gene names (incl. synonyms) and create a dict entry.
            #  At least one of these should map to something in GENCODE
            genes_names = []
            gene_name_uniprot_field = s.gene_name
            if gene_name_uniprot_field == "":
                continue
            genes_names_str = gene_name_uniprot_field.split(";")
            for gene_name_str in genes_names_str:
                gene_name_list = gene_name_str.split("=")
                if len(gene_name_list) > 1:
                    gene_name_list = gene_name_str.split("=")[1].split(",")
                    genes_names.extend(gene_name_list)

            for g in genes_names:
                gene_ids_2_entrynames[db][g].append(s.entry_name)

    found_genes = 0
    unfound_genes = 0
    gene_dict = {}
    for gene in genes:
        gene_dict[gene] = {}
        g = gene_dict[gene]
        if gene in gene_ids_2_entrynames['swiss']:
            data = swiss_data
            data_id = 'swiss'
            ##add swiss-prot data to dict
            gene_dict[gene] = create_gene_to_uniprot_dict(data, data_id, gene, gene_ids_2_entrynames)
            gene_dict[gene]['uniprot_db_source'] = data_id
            found_genes += 1

        elif gene in gene_ids_2_entrynames['trembl']:
            ##add trembl data to dict if none found in swiss-prot first
            data = trembl_data
            data_id = 'trembl'
            gene_dict[gene] = create_gene_to_uniprot_dict(data, data_id, gene, gene_ids_2_entrynames)
            gene_dict[gene]['uniprot_db_source'] = data_id
            g['uniprot_db_source'] = data_id
            found_genes += 1
        else:
            print(gene + " not found.")
            unfound_genes += 1

    print(str(found_genes) + " genes were found.")

    return gene_dict

def correct_uniprot_ids_hack(gene_dict, swiss_data):
    #replacement_dict = {
        #'TTN': ('TITIN_HUMAN', 'Q8WZ42'),
        #'MGMT': ('MGMT_HUMAN', 'P16455'),
        #'ETV2': ('ETV2_HUMAN', 'O00321'),
    #}
    genes_with_trembl_id = list()
    for k1,v in gene_dict.items():
        if v.get('uniprot_db_source','') == 'trembl':
            genes_with_trembl_id.append(k1)
    genes_easily_fixed = list()
    for g in genes_with_trembl_id:
        if g + '_HUMAN' in swiss_data:
            genes_easily_fixed.append(g)
    replacement_dict = dict()
    for g in genes_easily_fixed:
        replacement_dict[g] = (g + '_HUMAN', swiss_data[g + '_HUMAN'].accessions[0])
    replacement_dict['TTN'] = ('TITIN_HUMAN', 'Q8WZ42')
    for replacement in replacement_dict.keys():
        if replacement in gene_dict.keys():
            g = gene_dict[replacement]
            if 'alt_uniprot_entry_names' in g:
                g['alt_uniprot_entry_names'].append(g['uniprot_entry_name'])
            else:
                g['alt_uniprot_entry_names'] = [g['uniprot_entry_name']]
            if 'alt_uniprot_accessions' in g.keys():
                g['alt_uniprot_accessions'].append(g['uniprot_accession'])
            else:
                g['alt_uniprot_accessions'] = [g['uniprot_accession']]
            g['uniprot_entry_name'] = replacement_dict[replacement][0]
            g['uniprot_accession'] = replacement_dict[replacement][1]
            g['uniprot_db_source'] = 'swiss'
    return gene_dict

def parseWithPickle(fname, callableParsingFunction, pickleDir=""):
    ''' Pickle dir MUST include appended "/" '''
    pickleFilename = pickleDir + os.path.basename(fname) + ".pkl"
    if os.path.exists(pickleFilename):
        print("Loading pickled structure: " + str(pickleFilename))
        g = cPickle.load(file(pickleFilename, 'r'))
    else:
        print("Parsing...")
        g = callableParsingFunction(file(fname, 'r'))
        print("Writing pickle file: " + str(pickleFilename))
        cPickle.dump(g, file(pickleFilename, 'w'), protocol=0)
    return g



def renderUniprotNaturalVariationsTSV(genesDict, outputFilename):
    headers = ['gene', 'start_AA', 'end_AA','natural_variations']
    tsvWriter = csv.DictWriter(open(outputFilename, 'w'), headers, extrasaction='ignore', delimiter="\t", lineterminator="\n")
    tsvWriter.writeheader()

    # Used to prune duplicates
    genesAlreadyWritten = set()

    for b in genesDict.values():
        for k in b.keys():
            for g in b[k]:

                if g['gene'] in genesAlreadyWritten:
                    continue
                else:
                    genesAlreadyWritten.add(g['gene'])

                outrow = dict()
                outrow['gene'] = g['gene']

                if 'uniprot_features' in g.keys():
                    natural_variations = g['uniprot_features'].get('natural_variations',[])
                    for natvar in natural_variations:
                        outrow['natural_variations'] = natvar[3]
                        outrow['start_AA'] = natvar[1]
                        outrow['end_AA'] = natvar[2]
                        tsvWriter.writerow(outrow)


def renderSimpleUniprotTSV(gene_dict, outputFilename):
    
    headers = ['gene', 'uniprot_entry_name', 'DrugBank','alt_uniprot_accessions','uniprot_accession',
                      'GO_Biological_Process','GO_Cellular_Component','GO_Molecular_Function']
    tsvWriter = csv.DictWriter(open(outputFilename, 'w'), headers, extrasaction='ignore', delimiter="\t", lineterminator="\n")
    tsvWriter.writeheader()
    
    # Used to prune duplicates
    genesAlreadyWritten = set()
    genes = sorted(gene_dict.keys())

    for gene in genes:
        if gene in genesAlreadyWritten:
            print(gene + " was already written.")
            continue
        else:
            genesAlreadyWritten.add(gene)

        g = gene_dict[gene]
        if len(g.keys()) == 0:
            # no data for this gene
            continue

        g['gene'] = gene
        if 'GO' in g.keys():
            g['GO_Biological_Process'] = sorted([t[1] + " (" + t[0] + ")" for t in g['GO'].get('biological_process',[])])
            g['GO_Cellular_Component'] = sorted([t[1] + " (" + t[0] + ")" for t in g['GO'].get('cellular_component',[])])
            g['GO_Molecular_Function'] = sorted([t[1] + " (" + t[0] + ")" for t in g['GO'].get('molecular_function',[])])

        for col in g.keys():
            if isinstance(g[col], list):
                try:
                    tmp = set(g[col])
                    g[col] = "|".join(sorted(list(tmp),key=str.lower))
                except:
                    g[col] = str(g[col])
            # else:
            #     g[col] = str(g[col])

        tsvWriter.writerow(g)


if __name__ == '__main__':


    args = parseOptions()
    uniprot_swiss_fname = args.swiss_file
    uniprot_trembl_fname = args.trembl_file
    output_file = args.output_file

    swiss_data = parseWithPickle(uniprot_swiss_fname, parse_uniprot_data, "pickles/")
    trembl_data = parseWithPickle(uniprot_trembl_fname, parse_uniprot_data, "pickles/")
    
    gencode_ds_loc = expanduser(args.gencode_ds)
    gencode_ds = DatasourceFactory.createDatasource(configFilename=gencode_ds_loc, leafDir=os.path.dirname(gencode_ds_loc))
    gene_ids = gencode_ds.get_gene_symbols()
    print(str(len(gene_ids)) + " genes in datasource.")

    print "Adding uniprot IDs to genes..."
    genesDict = add_uniprot_ids_to_Genes(gene_ids, swiss_data, trembl_data)
    print "Adding uniprot data to genes..."
    genesDict = add_uniprot_data_to_Genes(genesDict, swiss_data, trembl_data)
    
    del trembl_data
    
    print "Correcting uniprot ids..."
    genesDict = correct_uniprot_ids_hack(genesDict, swiss_data)
    
    print("Writing simple_uniprot file... " + output_file)
    renderSimpleUniprotTSV(genesDict, output_file)

    # print("Writing natural variations file... " + output_file + ".natvar.tsv")
    # renderUniprotNaturalVariationsTSV(genesDict, output_file + ".natvar.tsv")
