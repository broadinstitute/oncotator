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


from argparse import ArgumentParser, RawDescriptionHelpFormatter
import csv
import logging
import os
from os.path import expanduser
from shove.core import Shove
from createUniprotTSVsGencode import parse_uniprot_data
from createUniprotProteinSeqsAlignments import generateTranscriptMuts
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.datasources.GenericTranscriptDatasource import GenericTranscriptDatasource
from createUniprotTSVsGencode import setup_logging

__author__ = 'lichtens'

def parse_with_shove(fname, callableParsingFunction, pickleDir=""):
    ''' Pickle dir MUST include appended "/" '''
    shoveFilename = pickleDir + os.path.basename(fname) + ".shv"
    if os.path.exists(shoveFilename):
        logging.getLogger(__name__).info("Loading shove structure: " + str(shoveFilename))
        g = Shove("file://" + shoveFilename, "simple://")
    else:
        logging.getLogger(__name__).info("Parsing...")
        tmpStruct = callableParsingFunction(file(fname, 'r'))
        logging.getLogger(__name__).info("Writing shove db: " + str(shoveFilename))
        ks = tmpStruct.keys()
        g = Shove("file://" + shoveFilename)
        for k in ks:
            del tmpStruct[k].references # May be causing an error later down the road.
            g[k] = tmpStruct[k]
    return g

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
    parser.add_argument("gencode_ds", type=str, help="GENCODE datasource config file. ")
    parser.add_argument("uniprot_tsv", type=str, help="Uniprot TSV file (used in a simple_uniprot datasource) file. ")
    parser.add_argument("-p", "--pickles", type=str, default="pickles/", help="Where to place temporary cached files. Default: %(default)s")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")

    args = parser.parse_args()
    return args

if __name__ == '__main__':

    setup_logging()
    args = parseOptions()
    uniprot_swiss_fname = args.swiss_file
    uniprot_trembl_fname = args.trembl_file
    output_file = args.output_file
    uniprot_tsv = args.uniprot_tsv

    pickles_dir = os.path.abspath(os.path.expanduser(args.pickles)) + "/"

    # Go through every record and create an entry for the
    outputHeaders = ["gene", "startAA", "endAA", "region", "site", "natural_variation","experimental_info"]
    tsvWriter = csv.DictWriter(open(output_file, 'w'), outputHeaders, extrasaction='ignore', delimiter="\t", lineterminator="\n")
    tsvWriter.writeheader()

    # TODO: Reduce code duplication

    swiss_data = parse_with_shove(uniprot_swiss_fname, parse_uniprot_data, pickles_dir)
    trembl_data = parse_with_shove(uniprot_trembl_fname, parse_uniprot_data, pickles_dir)

    # Use GENCODE datasource to get list of all possible genes
    gencode_ds_loc = expanduser(args.gencode_ds)
    gencode_ds = DatasourceFactory.createDatasource(configFilename=gencode_ds_loc, leafDir=os.path.dirname(gencode_ds_loc))

    # Use simple_uniprot TSV to get the uniprot_entry_names
    # Create the transcript to uniprot info mappings.  But take less RAM.  Given a gene, get the uniprot record.
    uniprotDS = GenericTranscriptDatasource(src_file=uniprot_tsv, title="UniProt", version="2014_12", geneColumnName="gene")

    # key is the uniprot_entry_name from the uniprotDS
    muts = generateTranscriptMuts(gencode_ds, uniprotDS)

    swissKeys = swiss_data.keys()
    tremblKeys = trembl_data.keys()

    featureTypeToAnnotation = {"SITE":"site", "VARIANT":"natural_variation", "COMPBIAS":"region" , "REGION":"region", "DOMAIN":"region", "CONFLICT":"experimental_info"}
    featureTypes = featureTypeToAnnotation.keys()
    ctr = 0
    numTranscriptsNotInUniprot = 0
    uniprotEntryNameKey = 'UniProt_uniprot_entry_name'
    txs_already_processed = set()
    records_already_processed = set()
    num_txs_with_same_data = 0
    for m in muts:
        if m['transcript_id'] in txs_already_processed:
            continue
        txs_already_processed.add(m['transcript_id'])
        ctr += 1
        if (ctr % 1000) == 0:
            logging.getLogger(__name__).info(str(ctr))

        uniprot_entry_name_key = m[uniprotEntryNameKey]
        if uniprot_entry_name_key in swissKeys:
            uniprot_record = swiss_data[uniprot_entry_name_key]

        elif uniprot_entry_name_key in tremblKeys:
            uniprot_record = trembl_data[uniprot_entry_name_key]
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

            if (feature, m['gene']) in records_already_processed:
                num_txs_with_same_data += 1
                continue

            records_already_processed.add((feature, m['gene']))

            row['startAA'] = feature[1]
            row['endAA'] = feature[2]
            row['gene'] = m['gene']
            row[annotation] = feature[3]
            tsvWriter.writerow(row)

    logging.getLogger(__name__).info("Could not get uniprot seq for " + str(numTranscriptsNotInUniprot) + " unique transcript IDs.")
    logging.getLogger(__name__).info("Attempted " + str(ctr) + " unique transcripts.")
    logging.getLogger(__name__).info("Identical feature data for txs with same gene: " + str(num_txs_with_same_data))

