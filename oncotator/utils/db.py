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


from oncotator.utils.io import read_delim
import os
import operator
import cPickle
import logging

def get_db_data(src_file, title, use_binary, index_mode, indexColumnNames='gene'):
    '''TODO: Finish documentation
	
	indexColumnNames is provided so that the column name can vary from datasource to datasource.  For multiple columns, use a comma
		separated string with no spaces (e.g. chr,start,end)
    '''
    # TODO: This is clunky and a case has to be added for each index_mode.
    if index_mode == 'gene' or index_mode == 'transcript':
        create_db_obj_from_txt_function = create_db_obj_from_txt_GENEINDEXED
    elif index_mode == 'genomic_pos' or index_mode=='gene_protein_pos':
        create_db_obj_from_txt_function = create_db_obj_from_txt_POSINDEXED

    if use_binary:
        bin_file = os.path.splitext(src_file)[0] + '.bin'

        if os.path.exists(bin_file):
            logging.info(("Loading %s..." % bin_file))
            db_obj, output_headers = cPickle.load(open(bin_file, 'rb'))
        else:
            db_obj, output_headers = create_db_obj_from_txt_function(src_file, title, indexColumnNames=indexColumnNames)
            logging.info("Saving indexed table as binary file for future use.")
            logging.info(bin_file)
            cPickle.dump((db_obj, output_headers), open(bin_file, 'wb'))

    else:
        db_obj, output_headers = create_db_obj_from_txt_function(src_file, title)

    return db_obj, output_headers

def read_data_and_add_title_str_to_each_header(data_0, title, indexCols='gene'):
	##Add Title string to beginning of each header
	indexColsList = indexCols.split(',')
	excludeList = indexColsList
	data = list()
	for d in data_0:
		t_dict = dict()
		for k in d:
			if k not in excludeList:
				new_key = '_'.join([title,k])
			else:
				new_key = k
			t_dict[new_key] = d[k]
		data.append(t_dict)
	return data

def create_db_obj_from_txt_GENEINDEXED(txt_fname, title, indexColumnNames='gene'):
	data_0, colnames = read_delim(txt_fname, return_list=False)
	data = read_data_and_add_title_str_to_each_header(data_0, title, indexColumnNames)
		
	new_colnames = ['%s_%s' % (title, c.strip()) for c in colnames if c not in [indexColumnNames]]

	annotation_table = dict()
	
	for d in data:
		if d[indexColumnNames] in annotation_table:
			print("ERROR: " + indexColumnNames + " indices must be unique: %s found multiple times in %s!" % (d[indexColumnNames], txt_fname))
			sys.exit(1)

		annotation_table[d[indexColumnNames]] = dict()
		for field in new_colnames:
			annotation_table[d[indexColumnNames]][field] = d[field].strip()
			
	return annotation_table, new_colnames
	
def create_db_obj_from_txt_POSINDEXED(txt_fname, title, indexColumnNames='chr,start,end'):
	data_0, colnames = read_delim(txt_fname, return_list=False)
	data = read_data_and_add_title_str_to_each_header(data_0, title, indexColumnNames)
	indexColumnList = indexColumnNames.split(',')
	
	if len(set(colnames).intersection( set(indexColumnList) )) <> len(set(indexColumnList)):
	
		print("ERROR: Invalid header values!")
		print("If indexing by genomic position, headers must include '" + indexColumnNames + "'.")
		raise Exception("If indexing by genomic position, headers must include '" + indexColumnNames + "'.")
	
	excludeList = indexColumnList
	excludeList.append('build')
		
	new_colnames = ['%s_%s' % (title, c) for c in colnames if c not in excludeList]

	annotation_table = dict()

	chr = indexColumnList[0]
	start = indexColumnList[1]
	end = indexColumnList[2]
	
	for d in data:
		if d[chr].startswith('chr'): d[chr] = d[chr][3:]

		d['bin'] = region2bin(int(d[start]), int(d[end]))
		d['start'] = d[start]
		d['end'] = d[end] 	
	data = sorted(data, key=operator.itemgetter(end))
	data = sorted(data, key=operator.itemgetter(start))
	data = sorted(data, key=operator.itemgetter(chr))

	for d in data:
		if d[chr] in annotation_table:
			if d['bin'] in annotation_table[d[chr]]:
				annotation_table[d[chr]][d['bin']].append(d)
			else:
				annotation_table[d[chr]][d['bin']] = [d]
		else:
			annotation_table[d[chr]] = {d['bin']:[d]}
			
	return annotation_table, new_colnames
	
def region2bin(beg, end):
	end = end-1
	if beg>>17 == end>>17: return ((1<<12)-1)/7 + (beg>>17)
	if beg>>20 == end>>20: return ((1<<9)-1)/7 + (beg>>20)
	if beg>>23 == end>>23: return ((1<<6)-1)/7 + (beg>>23)
	if beg>>26 == end>>26: return ((1<<3)-1)/7 + (beg>>26)
	return 0
	
def region2bins(beg, end):
	#end = end -1
	bins = [0]
	bins.extend(range(1+(beg>>26), 1+(end>>26)+1))
	bins.extend(range(9+(beg>>23), 9+(end>>23)+1))
	bins.extend(range(73+(beg>>20), 73+(end>>20)+1))
	bins.extend(range(585+(beg>>17), 585+(end>>17)+1))
	return bins
	
def get_overlapping_records(records, start, end, startColumnName='start', endColumnName='end'):
	out_records = list()
	for r in records:
		if test_overlap(start, end, int(r[startColumnName]), int(r[endColumnName])):
			out_records.append(r)
	return out_records
	
def get_binned_data(db_dict, chr, start, end):
	if chr not in db_dict:
		return []
	
	bins = region2bins(start, end)
	records = list()
	for b in bins:
		records.extend(db_dict[chr].get(b, []))
		
	return records

def test_overlap(a_st, a_en, b_st, b_en):
	if (a_st >= b_st and a_st <= b_en) or (a_en >= b_st and a_en <= b_en) or \
		(a_st <= b_st and a_en >= b_en) or (a_st >= b_st and a_en <= b_en):
		return True
	else:
		return False

def get_summary_output_string(lst_of_strs):
	#takes in list of strings and outputs formatted string for maf report
    out_string = ''
    summary_dict = dict()
    for s in lst_of_strs:
        if (s is None) or (s.strip() == ""):
            continue
        if s in summary_dict: summary_dict[s] += 1
        else: summary_dict[s] = 1
    summary_list = [(k, v) for k,v in summary_dict.items()]
	# summary_list.sort(key=operator.itemgetter(1), reverse=True)
    summary_list.sort()
    out_string = '|'.join(['%s' % (s[0]) for s in summary_list])
    return out_string