# LICENSE_GOES_HERE

class Options:
	"""
	A container for options that control how a annotation database set should be handled
	when used for mutation annotation.
	"""
	
	def __init__(self, **kwargs):
		self.datasources = []
		
	def add_datasource(self, datasource):
		self.datasources.append(datasource)
