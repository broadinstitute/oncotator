# LICENSE_GOES_HERE



'''
Created on Nov 15, 2012

@author: lichtens
'''

class GafDatasourceException(Exception):
    '''
    TODO: Documentation
    '''


    def __init__(self, value):
        '''
        
        '''
        self.value = value

    def __str__(self):
         return repr(self.value)