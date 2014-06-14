# LICENSE_GOES_HERE
from DatasourceCreator import DatasourceCreator
from TabixIndexer import TabixIndexer
import ConfigParser
import os


class TabixIndexedVcfDatasourceCreator(DatasourceCreator):

    def __init__(self):
        pass

    def createDatasource(self, destDir, ds_file, configFilename, ds_type, ds_name, ds_version,
                         ds_match_mode, index_column_names=None, indexCols=None, annotation_column_names=None):
        # Create database
        baseDSFile = self._createDatabase(destDir, ds_file)

        # Create config file
        self._createConfigFile(configFilename, baseDSFile, ds_type, ds_name, ds_version, ds_match_mode)

    def _createDatabase(self, destDir, ds_file):
        tabixIndexedFile = TabixIndexer.index(destDir=destDir, inputFilename=ds_file, preset="vcf")
        baseDSFile = os.path.basename(tabixIndexedFile)

        return baseDSFile

    def _createConfigFile(self, configFilename, baseDSFile, ds_type, ds_name, ds_version, ds_match_mode):
        """


        :param configFilename: configuration filename
        :param baseDSFile: data source base filename
        :param ds_type: data source type
        :param ds_name: data source title
        :param ds_version: data source version
        :param ds_match_mode: describes how to annotate mutations from an indexed tsv or indexed vcf datasources
        """
        config = ConfigParser.ConfigParser()
        filePtr = open(configFilename, 'w')
        config.add_section("general")
        config.set("general", "version", ds_version)
        config.set("general", "title", ds_name)
        config.set("general", "type", ds_type)
        config.set("general", "src_file", baseDSFile)
        config.set("general", "match_mode", ds_match_mode)
        config.write(filePtr)
        filePtr.close()