# LICENSE_GOES_HERE


"""
Created on Oct 22, 2012

@author: lichtens
"""

class MafliteMissingRequiredHeaderException(Exception):
    """
    Exception class arising when a maflite file does not have all of the required headers.
    """


    def __init__(self, value):
        """
        
        """
        self.value = value

    def __str__(self):
        return repr(self.value)
