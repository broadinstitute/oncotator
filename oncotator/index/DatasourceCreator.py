
from abc import abstractmethod


class DatasourceCreator(object):
    """
    This is the base class for creating databases and congif files.
    """

    @abstractmethod
    def createDatasource(self, destDir, ds_file, index_column_names, column_names, configFilename, ds_type, ds_name,
                         ds_version, ds_match_mode, annotation_column_names, indexCols):

        raise NotImplementedError