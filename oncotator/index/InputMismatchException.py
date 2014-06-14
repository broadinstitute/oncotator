# LICENSE_GOES_HERE


class InputMismatchException(Exception):
    """
    Exception class arising when a tsv file does not have all of the required headers.
    """

    def __init__(self, value):
        """
        
        """
        self.value = value

    def __str__(self):
        return repr(self.value)
