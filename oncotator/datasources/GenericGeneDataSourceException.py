# LICENSE_GOES_HERE

class GenericGeneDataSourceException(Exception):
    def __init__(self, str):
        """

        """
        self.value = str

    def __str__(self):
        return repr(self.value)