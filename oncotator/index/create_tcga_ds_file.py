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
import sys
import pandas
import numpy


def add_full_mut_key_columns(input_df, chr_col_name='Chromosome', start_col_name='Start_position',
                             end_col_name='End_position', ref_col_name='Reference_Allele',
                             alt_col_name='Tumor_Seq_Allele1'):
    chrom_series, start_series = input_df[chr_col_name].apply(str), input_df[start_col_name].apply(str)
    end_series, ref_series = input_df[end_col_name].apply(str), input_df[ref_col_name]
    tum_series = input_df[alt_col_name]
    input_df[
        'full_mut_key'] = chrom_series + '_' + start_series + '_' + end_series + '_' + ref_series + '_' + tum_series


def main(tcga_mutations_fname, output_fname):
    print "Loading tcga data..."
    tcga_df = pandas.read_table(tcga_mutations_fname, skiprows=1, skipfooter=1)
    get_alt_allele = lambda x: x['tumor_seq_allele_1'] if x['tumor_seq_allele_1'] != x['reference_allele'] else x[
        'tumor_seq_allele_2']
    tcga_df['alt_allele'] = tcga_df.apply(get_alt_allele, axis=1)

    add_full_mut_key_columns(tcga_df, chr_col_name='chromosome', start_col_name='start_position',
                             end_col_name='end_position',
                             ref_col_name='reference_allele', alt_col_name='alt_allele')

    tcga_mut_counts = dict(tcga_df.full_mut_key.value_counts())
    all_tcga_diseases = set(tcga_df.tcga_disease)
    total_disease_sample_counts = dict()
    for tcga_disease in all_tcga_diseases:
        total_disease_sample_counts[tcga_disease] = len(
            set(tcga_df[tcga_df.tcga_disease == tcga_disease].tumor_sample_barcode))

    print "Indexing tcga data..."
    indexed_tcga_data = dict()
    for idx, row in tcga_df.iterrows():
        disease = row['tcga_disease']
        key = row['full_mut_key']
        if key not in indexed_tcga_data:
            indexed_tcga_data[key] = dict()

        if disease not in indexed_tcga_data[key]:
            indexed_tcga_data[key][disease] = 0

        indexed_tcga_data[key][disease] += 1

    print "Parsing tcga data..."
    tcga_mut_count = dict()
    tcga_muts_disease_names = dict()
    tcga_muts_disease_counts = dict()
    tcga_muts_disease_percents = dict()
    all_mut_keys = set()
    for mut_key in indexed_tcga_data:
        all_mut_keys.add(mut_key)
        disease_data = sorted(indexed_tcga_data[mut_key].items(), key=lambda x: x[1], reverse=True)
        tcga_muts_disease_names[mut_key] = ':'.join([d[0] for d in disease_data])
        tcga_muts_disease_counts[mut_key] = ':'.join([str(d[1]) for d in disease_data])
        get_mut_percent = lambda d: '%.2f' % (float(d[1]) / total_disease_sample_counts[d[0]] * 100)
        tcga_muts_disease_percents[mut_key] = ':'.join([get_mut_percent(d) for d in disease_data])
        tcga_mut_count[mut_key] = sum(d[1] for d in disease_data)

    print "Writing parsed tcga data..."
    output_headers = ['chr', 'start', 'end', 'ref_allele', 'alt_alele', 'tcga_mut_count', 'tcga_muts_disease_names',
                      'tcga_muts_disease_counts', 'tcga_muts_disease_percents']
    out_fh = open(output_fname, 'w')
    out_fh.write('\t'.join(output_headers) + '\n')

    all_mut_keys = sorted(all_mut_keys)
    for mut_key in all_mut_keys:
        output_fields = mut_key.split('_')
        output_fields.extend(
            [str(tcga_mut_count[mut_key]), tcga_muts_disease_names[mut_key], tcga_muts_disease_counts[mut_key],
             tcga_muts_disease_percents[mut_key]])
        out_fh.write('\t'.join(output_fields) + '\n')

    out_fh.close()


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Usage: python create_tcga_ds_file.py input_tcga_file output_oncotator_ds_file"
        print "Example: create_tcga_ds_file.py mutations.2013042800.txt tcga_2013042800.tsv"
        sys.exit(1)

    tcga_mutations_fname, output_fname = sys.argv[1:]
    main(tcga_mutations_fname, output_fname)



