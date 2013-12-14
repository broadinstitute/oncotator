from DatasourceCreator import DatasourceCreator
from ConfigParser import ConfigParser
import os
import shutil


class GenericTsvDatasourceCreator(DatasourceCreator):

    def __init__(self):
        pass

    def createDatasource(self, destDir, ds_file, index_column_names=None, column_names=None):
        baseDSFile = os.path.basename(ds_file)
        shutil.copy(ds_file, os.path.join(destDir, baseDSFile))

        return baseDSFile

    def createConfigFile(self, configFilename, baseDSFile, ds_type, ds_name, ds_version, indexCols,
                         column_names=None, annotation_column_names=None):
        """


        :param configFilename: configuration filename
        :param baseDSFile: data source base filename
        :param ds_type: data source type
        :param ds_name: data source title
        :param ds_version: data source version
        :param column_names: column names in the input data source file
        :param annotation_column_names: column names whose values are used for annotation
        :param indexCols: named tuple consisting of index column type and corresponding column names
        """
        config = ConfigParser()
        filePtr = open(configFilename, 'w')
        config.add_section("general")
        config.set("general", "version", ds_version)
        config.set("general", "title", ds_name)
        config.set("general", "type", ds_type)
        config.set("general", "src_file", baseDSFile)
        config.set("general", indexCols.type, indexCols.names)
        config.write(filePtr)
        filePtr.close()