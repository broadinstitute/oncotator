# LICENSE_GOES_HERE


from argparse import ArgumentParser, RawDescriptionHelpFormatter
import csv
import os
from reportlab.platypus.doctemplate import onDrawStr
import shutil
from createUniprotTSVs import parse_uniprot_data
from createUniprotProteinSeqsAlignments import parseWithShove,generateTranscriptMuts
from oncotator.datasources import Generic_Gene_DataSource, Gaf
from oncotator.utils.TsvFileSorter import TsvFileSorter

__author__ = 'lichtens'


def parseOptions():

    # Process arguments
    desc = ''' Create a tsv file with uniprot data (Natural Variations, Regions, Sites, and Experimental Info).

    This script requires 3GB RAM.'''
    epilog = '''

    Uniprot TSV file will have columns:
        gene
        uniprot_entry_name
        DrugBank
        alt_uniprot_accessions
        uniprot_accession
        GO_Biological_Process
        GO_Cellular_Component
        GO_Molecular_Function

    This can be generated with the createUniprotTSVs.py.  The only columns needed for this script are:
     gene
     uniprot_entry_name

    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("swiss_file", type=str, help="SwissProt file. ")
    parser.add_argument("trembl_file", type=str, help="TREMBL file. ")
    parser.add_argument("gaf_file", type=str, help="GAF file. ")
    parser.add_argument("gaf_transcript_file", type=str, help="GAF transcript file. E.g. UCSCgene.Jan2012.fa")
    parser.add_argument("uniprot_tsv", type=str, help="Uniprot TSV file (used in a simple_uniprot datasource) file. ")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")

    args = parser.parse_args()
    return args

if __name__ == '__main__':

    args = parseOptions()
    uniprot_swiss_fname = args.swiss_file
    uniprot_trembl_fname = args.trembl_file
    output_file = args.output_file
    gaf_file = args.gaf_file
    uniprot_tsv = args.uniprot_tsv
    gaf_transcript_file = args.gaf_transcript_file

    # Go through every record and create an entry for the
    outputHeaders = ["gene", "startAA", "endAA", "region", "site", "natural_variation","experimental_info"]
    tsvWriter = csv.DictWriter(open(output_file, 'w'), outputHeaders, extrasaction='ignore', delimiter="\t", lineterminator="\n")
    tsvWriter.writeheader()

    # TODO: Remove hardcoded paths
    # TODO: Reduce code duplication

    swiss_data = parseWithShove(uniprot_swiss_fname, parse_uniprot_data, "/bulk/pickles/")
    trembl_data = parseWithShove(uniprot_trembl_fname, parse_uniprot_data, "/bulk/pickles/")

    # Use GAF datasource to get list of all possible genes
    gafDS = Gaf(gaf_file, gaf_transcript_file)

    # Use simple_uniprot TSV to get the uniprot_entry_names
    # Create the gene to uniprot info mappings.  But take less RAM.  Given a gene, get the uniprot record.
    uniprotDS = Generic_Gene_DataSource(src_file=uniprot_tsv, title="UniProt", version="2011_09", geneColumnName="gene")

    # key is the uniprot_entry_name from the uniprotDS
    muts = generateTranscriptMuts(gafDS, uniprotDS)

    swissKeys = swiss_data.keys()
    tremblKeys = trembl_data.keys()

    featureTypeToAnnotation = {"SITE":"site", "VARIANT":"natural_variation", "COMPBIAS":"region" , "REGION":"region", "DOMAIN":"region", "CONFLICT":"experimental_info"}
    featureTypes = featureTypeToAnnotation.keys()
    ctr = 0
    numTranscriptsNotInUniprot = 0
    uniprotEntryNameKey = 'UniProt_uniprot_entry_name'
    for m in muts:
        ctr += 1
        if (ctr % 1000) == 0:
            print(str(ctr))

        if m[uniprotEntryNameKey] in swissKeys:
            uniprot_record = swiss_data[m[uniprotEntryNameKey]]

        elif m[uniprotEntryNameKey] in tremblKeys:
            uniprot_record = trembl_data[m[uniprotEntryNameKey]]
        else:
            numTranscriptsNotInUniprot += 1
            continue

        features = uniprot_record.features

        for feature in features:
            # If there is a blank value, then no need to continue.
            if (feature[3] is None) or (feature[3].strip() == ""):
                continue

            # Each feature is a tuple (type, start,end,description)
            if feature[0] not in featureTypes:
                continue
            annotation = featureTypeToAnnotation[feature[0]]

            row = dict()
            for hdr in outputHeaders:
                row[hdr] = ""

            # Test to make sure that both start and end can be cast as integers.  If not, silently skip.
            try:
                tmp1 = int(feature[1])
                tmp2 = int(feature[2])
            except:
                continue

            row['startAA'] = feature[1]
            row['endAA'] = feature[2]
            row['gene'] = m['gene']

            row[annotation] = feature[3]
            tsvWriter.writerow(row)

    print("Could not get uniprot seq for " + str(numTranscriptsNotInUniprot) + " transcripts.")
    print("Attempted " + str(ctr) + " muts")

    print("Creating tabix index")
    print("Creating copy of tsv file (" + output_file + ") ...")
    tabixBasedFilename = output_file + ".copy.tsv"
    shutil.copyfile(output_file, tabixBasedFilename)

    print("Sorting ...")
    tsvFileSorter = TsvFileSorter(fieldNames=['gene','startAA', 'endAA'])
    tsvFileSorter.sortFile(tabixBasedFilename, tabixBasedFilename + ".sorted")
    print("Creating actual index ...")

    # swiss_data[key].features
    #  For each feature, position 0 is name.
    #  Look for "SITE" (site), "VARIANT" (natural_variation),
    # "COMPBIAS" or "REGION" or "DOMAIN"? (region)
    #   create a line for each entry
    # Then add trembl, but only if swiss_prot has not already covered it
    #
    #   Verify with old oncotator code?
    pass