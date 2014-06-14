# LICENSE_GOES_HERE
from os.path import expanduser
from tempfile import mkdtemp
from oncotator.DatasourceFactory import DatasourceFactory


'''
Preprocessing script to create UniProt protein sequence alignment to the GENCODE Txs.
'''
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import os
import re
from oncotator.MutationData import MutationData
from createUniprotTSVs import parse_uniprot_data
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
    parser.add_argument("simple_uniprot_ds", type=str, help="Location of the simple uniprot datasource -- E.g. /bulk/dbDir/simple_uniprot/hg19/simple_uniprot.config")
    parser.add_argument("blast_exe", type=str, help="Location of blast executable")
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

    temp_refseq_fasta_fh,temp_refseq_fasta_fname = tempfile.mkstemp(suffix='.tmp',prefix='annotation.',dir=tmp_dir,text=True)
    temp_uniprot_fasta_fh,temp_uniprot_fasta_fname = tempfile.mkstemp(suffix='.tmp',prefix='annotation.',dir=tmp_dir,text=True)
    temp_align_results_fname = os.path.join(tmp_dir, '%s_%s.alignment' % (seq1_name, seq2_name))
    write_fasta(temp_refseq_fasta_fname, 't1', NP_seq)
    write_fasta(temp_uniprot_fasta_fname, 't2', uni_seq)
    if bl2seq_path.endswith('bl2seq'):
        cmd = [bl2seq_path, '-p', 'blastp', '-F', 'F', '-i', temp_refseq_fasta_fname, '-j',
               temp_uniprot_fasta_fname, '-o', temp_align_results_fname]
    else:
        cmd = [bl2seq_path, '-query', temp_refseq_fasta_fname, '-subject', temp_uniprot_fasta_fname,
               '-out', temp_align_results_fname, '-seg', 'no']
    subprocess.check_output(cmd, close_fds=True)
    ff = open(temp_align_results_fname)
    alignment_data = ff.readlines()
    ff.close()
    db[alignment_key] = alignment_data
    ##cleanup
    os.close(temp_refseq_fasta_fh)
    os.close(temp_uniprot_fasta_fh)
    for t in [temp_refseq_fasta_fname, temp_uniprot_fasta_fname]: #, temp_align_results_fname]:
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
    uniprot_ds_loc = expanduser(args.simple_uniprot_ds)
    blast_exe = args.blast_exe

    gencode_ds = DatasourceFactory.createDatasource(configFilename=gencode_ds_loc, leafDir=os.path.dirname(gencode_ds_loc))
    uniprotDS = DatasourceFactory.createDatasource(configFilename=uniprot_ds_loc, leafDir=os.path.dirname(uniprot_ds_loc))

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
        runAlignment(tx_id, uniprot_entry_key, tx_protein_seq, uniprot_seq, tmp_dir, blast_exe, alignmentDB)

    print("Could not get protein seq for " + str(numNotInProteinSeqs) + " transcripts.")
    print("Could not get uniprot seq for " + str(numTranscriptsNotInUniprot) + " transcripts.")
    print("Attempted " + str(ctr) + " muts")

