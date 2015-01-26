from DatasourceCreator import DatasourceCreator
from ConfigParser import ConfigParser
import os
import shutil

class BigWigDatasourceCreator(DatasourceCreator):

    def createDatasource(self, destDir, ds_file, configFilename, ds_type, ds_name, ds_version, **kwargs):
        # Create database
        baseDSFile = self._createDatabase(destDir, ds_file)

        # Create config file
        self._createConfigFile(configFilename, baseDSFile, ds_type, ds_name, ds_version)

    def _createDatabase(self, destDir, ds_file):
        baseDSFile = os.path.basename(ds_file)
        shutil.copy(ds_file, os.path.join(destDir, baseDSFile))

        return baseDSFile

    def _createConfigFile(self, configFilename, baseDSFile, ds_type, ds_name, ds_version):
        """


        :param configFilename: configuration filename
        :param baseDSFile: data source base filename
        :param ds_type: data source type
        :param ds_name: data source title
        :param ds_version: data source version
        """
        config = ConfigParser()
        filePtr = open(configFilename, 'w')
        config.add_section("general")
        config.set("general", "version", ds_version)
        config.set("general", "title", ds_name)
        config.set("general", "type", ds_type)
        config.set("general", "src_file", baseDSFile)
        config.write(filePtr)
        filePtr.close()