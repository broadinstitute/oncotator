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
from tempfile import mkdtemp
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.datasources.GenericGeneDatasource import GenericGeneDatasource


'''
Preprocessing script to create UniProt protein sequence alignment to the GENCODE Txs.
'''
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import os
import re
from oncotator.MutationData import MutationData
from createUniprotTSVsGencode import parse_uniprot_data
import tempfile
import subprocess

from shove import Shove
import logging

class UniprotDBError(Exception):
    def __init__(self, err):
        self.err = err
    def __str__(self):
        return repr(self.err)
        
class UniprotPositionError(Exception):
    def __init__(self, err):
        self.err = err
    def __str__(self):
        return repr(self.err)

class ResidueMismatchError(Exception):
    def __init__(self, err):
        self.value = err
    def __str__(self):
        return repr(self.err)

class UnableToMapToSeq2Error(Exception):
    def __init__(self, err):
        self.value = err
    def __str__(self):
        return repr(self.err)

class InputResidueMismatchError(Exception):
    def __init__(self, err):
        self.value = err
    def __str__(self):
        return repr(self.err)        

class OutOfBoundsError(Exception):
    def __init__(self, err):
        self.value = err
    def __str__(self):
        return repr(self.err)


def parseOptions():
    
    # Process arguments
    desc = ''' Create the data for a uniprot aa xform datasource.  This script is meant for developers, not end users.'''
    epilog = '''
     What is being done here?

     Getting a list of all transcripts from a GENCODE datasource.
     For each transcript, getting the gene from GENCODE datasource and uniprot entry name from a uniprot datasource (created with createUniprotTSVs.py)

     For each transcript get the unpiprot protein sequence (using gene) and blast it against the transcript protein_seq.
     The output of blast is written to a database.

     This script still needs testing.
    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("swiss_file", type=str, help="SwissProt file. ")
    parser.add_argument("trembl_file", type=str, help="TREMBL file. ")
    parser.add_argument("gencode_ds", type=str, help="Location of the GENCODE datasource config file-- E.g. /bulk/dbDir/gencode_ds/hg19/gencode_ds.config")
    parser.add_argument("uniprot_tsv", type=str, help="Location of the simple uniprot datasource tsv-- E.g. /bulk/dbDir/simple_uniprot/hg19/simple_uniprot.tsv")
    parser.add_argument("blast_exe", type=str, help="Location of blast executable")
    parser.add_argument("out_tx_matches", type=str, help="Output file for perfect match transcripts -- where uniprot and tx match exactly on protein seq.")
    parser.add_argument("--temp_pickle_store", type=str, help="Store uniprot pickles in the given directory.  Good if you want to run this utility multiple times -- effectively caches some intermediate files")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")

    args = parser.parse_args()
    return args

def parseWithShove(fname, callableParsingFunction, pickleDir=""):
    ''' Pickle dir MUST include appended "/" '''
    shoveFilename = pickleDir + "/" + os.path.basename(fname) + ".shv"
    if os.path.exists(shoveFilename):
        print("Loading shove structure: " + str(shoveFilename))
        g = Shove("file://" + shoveFilename, "simple://")
    else:
        print("Parsing...")
        tmpStruct = callableParsingFunction(file(fname, 'r'))
        print("Writing shove db: " + str(shoveFilename))
        ks = tmpStruct.keys()
        g = Shove("file://" + shoveFilename)
        for k in ks:
            del tmpStruct[k].references # May be causing an error later down the road.
            g[k] = tmpStruct[k]
    return g


def map_uni_pos(qobj, sobj, aa_pos):
    qoff = 0
    soff = 0
    for i in range(len(qobj.group(2))):
        if qobj.group(2)[i] == '-':
            qoff -= 1
        qpos = i + int(qobj.group(1)) + qoff

        if sobj.group(2)[i] == '-':
            soff -= 1
        spos = i + int(sobj.group(1)) + soff

        if qpos == aa_pos:
            return spos, qobj.group(2)[i] , sobj.group(2)[i]


def get_uni_pos(fh, AA):
    new_pos = 0
    query_AA = ''
    uni_AA = ''
    pat1 = re.compile(r'Query:\s(\d+)\s+([\w\-\*]+)\s(\d+)')
    pat2= re.compile(r'Sbjct:\s(\d+)\s+([\w\-\*]+)\s(\d+)')
    q = 0
    qline = None
    for line in fh:
        if q == 1 and line.startswith('Sbjct: '):
            sline = pat2.search(line)
            new_pos, query_AA, uni_AA = map_uni_pos(qline, sline, AA)
            break
        if line.startswith('Query: '):
            qline = pat1.search(line)
            #print '{0}\t{1}'.format(qline.group(1), qline.group(3))
            if int(qline.group(1)) <= AA and int(qline.group(3)) >= AA:
                q = 1
    return new_pos, query_AA, uni_AA

def runAlignment(seq1_name, seq2_name,NP_seq,uni_seq,tmp_dir,bl2seq_path, db):
    alignment_key = '%s' % (seq1_name)

    temp_gencode_fasta_fh,temp_gencode_fasta_fname = tempfile.mkstemp(suffix='.tmp',prefix='annotation.',dir=tmp_dir,text=True)
    temp_uniprot_fasta_fh,temp_uniprot_fasta_fname = tempfile.mkstemp(suffix='.tmp',prefix='annotation.',dir=tmp_dir,text=True)
    temp_align_results_fname = os.path.join(tmp_dir, '%s_%s.alignment' % (seq1_name, seq2_name))
    write_fasta(temp_gencode_fasta_fname, 't1', NP_seq)
    write_fasta(temp_uniprot_fasta_fname, 't2', uni_seq)
    if bl2seq_path.endswith('bl2seq'):
        cmd = [bl2seq_path, '-p', 'blastp', '-F', 'F', '-i', temp_gencode_fasta_fname, '-j',
               temp_uniprot_fasta_fname, '-o', temp_align_results_fname]
    else:
        cmd = [bl2seq_path, '-query', temp_gencode_fasta_fname, '-subject', temp_uniprot_fasta_fname,
               '-out', temp_align_results_fname, '-seg', 'no']
    subprocess.check_output(cmd, close_fds=True)
    ff = open(temp_align_results_fname)
    alignment_data = ff.readlines()
    ff.close()
    db[alignment_key] = alignment_data
    ##cleanup
    os.close(temp_gencode_fasta_fh)
    os.close(temp_uniprot_fasta_fh)
    for t in [temp_gencode_fasta_fname, temp_uniprot_fasta_fname]: #, temp_align_results_fname]:
        os.remove(t)

def write_fasta(fname, header, seq):
    OUT = open(fname, 'w')
    OUT.write('>' + header + '\n')
    OUT.write(seq + '\n')
    OUT.close()

def generateTranscriptMuts(gafDS,uniprotDS):
    tDict = gafDS.getTranscriptDict()
    for transcriptID in tDict.keys():
        m = MutationData()
        m.createAnnotation('gene', tDict[transcriptID]['gene'])
        m.createAnnotation('transcript_id', transcriptID)
        m = uniprotDS.annotate_mutation(m)
        yield m

if __name__ == '__main__':
    args = parseOptions()
    uniprot_swiss_fname = expanduser(args.swiss_file)
    uniprot_trembl_fname = expanduser(args.trembl_file)
    output_file = args.output_file
    gencode_ds_loc = expanduser(args.gencode_ds)

    out_tx_matches = args.out_tx_matches
    out_tx_matches_fp = file(out_tx_matches, 'w')

    uniprot_tsv = expanduser(args.uniprot_tsv)
    blast_exe = args.blast_exe

    gencode_ds = DatasourceFactory.createDatasource(configFilename=gencode_ds_loc, leafDir=os.path.dirname(gencode_ds_loc))
    #uniprotDS = DatasourceFactory.createDatasource(configFilename=uniprot_ds_loc, leafDir=os.path.dirname(uniprot_ds_loc))
    uniprotDS = GenericGeneDatasource(src_file=uniprot_tsv, title="UniProt", version="2014_12", geneColumnName="gene")

    tmp_dir = args.temp_pickle_store
    if tmp_dir is None:
        tmp_dir = mkdtemp(prefix="onco_unipickles_")
    swiss_data = parseWithShove(uniprot_swiss_fname, parse_uniprot_data, tmp_dir)
    trembl_data = parseWithShove(uniprot_trembl_fname, parse_uniprot_data, tmp_dir)
    alignmentDB = Shove("file://" + output_file, "simple://")

    # Go through each transcript
    txs = gencode_ds.getTranscriptDict()
    tx_ids = txs.keys()
    num_tx_ids = len(tx_ids)
    swissKeys = swiss_data.keys()
    tremblKeys = trembl_data.keys()

    uniprotEntryNameKey = 'UniProt_uniprot_entry_name'

    numNotInProteinSeqs = 0
    numTranscriptsNotInUniprot = 0
    ctr = 0
    for tx_id in tx_ids:
        ctr += 1
        if (ctr % 2000) == 0:
            print(str(ctr) + "/" + str(num_tx_ids))

        tx_protein_seq = txs[tx_id].get_protein_seq()
        if tx_protein_seq is None or tx_protein_seq.strip() == "" or tx_protein_seq.strip() == "*":
            numNotInProteinSeqs += 1
            continue


        # Create a fake dummy mutation and annotate the gene and the simple_uniprot info
        m = MutationData()
        m.createAnnotation('gene', txs[tx_id].get_gene())
        m = uniprotDS.annotate_mutation(m)
        uniprot_entry_key = m[uniprotEntryNameKey]
        if uniprot_entry_key in swissKeys:
            uniprot_record = swiss_data[uniprot_entry_key]

        elif uniprot_entry_key in tremblKeys:
            uniprot_record = trembl_data[uniprot_entry_key]
        else:
            numTranscriptsNotInUniprot += 1
            continue
        uniprot_seq = uniprot_record.sequence

        # print(m['transcript_id'] + " " + m[uniprotEntryNameKey])
        # "/bulk/blast-2.2.26/bin/bl2seq" is blast_exe for Lee's laptop VM

        # When doing the comparison, tx protein includes stop codon at the end, uniprot does not.
        if tx_protein_seq[0:-1] == uniprot_seq:
            out_tx_matches_fp.write(tx_id + "\n")

        runAlignment(tx_id, uniprot_entry_key, tx_protein_seq, uniprot_seq, tmp_dir, blast_exe, alignmentDB)

    print("Could not get protein seq for " + str(numNotInProteinSeqs) + " transcripts.")
    print("Could not get uniprot seq for " + str(numTranscriptsNotInUniprot) + " transcripts.")
    print("Attempted " + str(ctr) + " muts")

