# LICENSE_GOES_HERE

class MissingRequiredAnnotationException(Exception):
    """
    Exception class arising when a required annotation is missing.
    """

    def __init__(self, value):
        """
        
        """
        self.value = value

    def __str__(self):
        return repr(self.value)