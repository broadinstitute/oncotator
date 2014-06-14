# LICENSE_GOES_HERE


'''
Created on Mar 4, 2013

@author: lichtens
'''
from oncotator.datasources.Datasource import Datasource


class MockExceptionThrowingDatasource(Datasource):
    '''
    Throws an exception whenever you initialize it or try to annotate
    '''


    def __init__(self,src_file="", title='', version=""):
        '''
        Constructor
        '''
        raise NotImplementedError("This will never be implemented.  This class throws exceptions.")
    
    def annotate_mutations(self, muts):
        raise NotImplementedError("Annotation from this datasource will never be implemented.  This class throws exceptions.")