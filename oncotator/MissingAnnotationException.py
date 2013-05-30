
class MissingAnnotationException(Exception):
    '''
    Thrown when a datasource needs an annotation (e.g. Generic_GeneDatasource needs "gene") that is not present
        in a mutation.
    '''


    def __init__(self, value):
        '''

        '''
        self.value = value

    def __str__(self):
        return repr(self.value)