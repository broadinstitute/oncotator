# LICENSE_GOES_HERE



def gen_cat(sources):
	""" Generic cat generator for file input.
	"""
	for s in sources:
		for item in s:
			yield item.strip('\n')
			
def read_delim(infile, return_list=False, headers=None, delim='\t'):
	""" Utility function for reading txt files.
		Headers are inferred if not provided.
		Returns data generator and list of headers.
		Use return_list=True to return list of data instead of generator.
	"""
	
	IN = open(infile, 'r')
	return read_delim_from_fileobject(IN, return_list, headers=headers, delim=delim)
	
def read_delim_from_fileobject(file, return_list=False, headers=None, delim='\t'):
	if headers:
		colnames = headers
	else:
		headers = file.readline()
		while headers.startswith('#'):
			headers = file.readline()
		colnames = headers.strip('\n').strip('\r').split(delim)
	datalines = gen_cat([file])
	fields = (m.split('\t') for m in datalines)
	data = (dict(zip(colnames,f)) for f in fields)
	if return_list:
		data = list(data)
	return data, colnames

