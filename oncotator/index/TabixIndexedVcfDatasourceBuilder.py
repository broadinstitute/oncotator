from DatasourceBuilder import DatasourceBuilder
from oncotator.index.TabixIndexer import TabixIndexer
import ConfigParser
import os


class TabixIndexedVcfDatasourceBuilder(DatasourceBuilder):

    def __init__(self):
        pass

    def createDatasource(self, destDir, ds_file, index_column_names=None, column_names=None):
        tabixIndexedFile = TabixIndexer.index(destDir=destDir, inputFilename=ds_file, preset="vcf")
        baseDSFile = os.path.basename(tabixIndexedFile)

        return baseDSFile

    def createConfigFile(self, configFilename, baseDSFile, ds_type, ds_name, ds_version, column_names=None,
                         annotation_column_names=None, indexCols=None):
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
        config = ConfigParser.ConfigParser()
        filePtr = open(configFilename, 'w')
        config.add_section("general")
        config.set("general", "version", ds_version)
        config.set("general", "title", ds_name)
        config.set("general", "type", ds_type)
        config.set("general", "src_file", baseDSFile)
        config.write(filePtr)
        filePtr.close()