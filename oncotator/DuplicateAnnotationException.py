# LICENSE_GOES_HERE


'''
Created on Oct 22, 2012

@author: lichtens
'''

class DuplicateAnnotationException(Exception):
    '''
    Exception class arising when an annotation with the same name is being written.
    '''


    def __init__(self, value):
        '''
        
        '''
        self.value = value

    def __str__(self):
         return repr(self.value)
