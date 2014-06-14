# LICENSE_GOES_HERE


'''
Created on Nov 13, 2012

@author: lichtens
'''


class MutationValidationFailureException(Exception):
    '''
    Exception class arising when a MutationData object fails validation
    '''


    def __init__(self, value):
        '''
        
        '''
        self.value = value

    def __str__(self):
        return repr(self.value)
