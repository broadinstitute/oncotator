# LICENSE_GOES_HERE

from oncotator.datasources.Datasource import Datasource


class PositionTransformingDatasource(Datasource):
    """ Given a coordinate from system A, translates it to system B.

    Useful for wrapping mappings, such as GAF transcript protein position to uniprot protein position.

    """

    def __init__(self, src_file, title='', version=None, use_binary=True):
        raise NotImplementedError("Not implemented yet")

    def annotate_mutation(self, mutation):
        raise NotImplementedError("Not implemented yet")