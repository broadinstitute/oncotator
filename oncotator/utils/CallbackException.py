# LICENSE_GOES_HERE


"""
Created on Jan 9, 2013

@author: mgupta
"""

class CallbackException(Exception):
    """
    Exception class arising when input VCF has syntax errors.
    """

    def __init__(self, value):
        """
        
        """
        self.value = value

    def __str__(self):
        return repr(self.value)