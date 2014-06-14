# LICENSE_GOES_HERE


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