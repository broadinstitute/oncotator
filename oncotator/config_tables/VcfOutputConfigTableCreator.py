# LICENSE_GOES_HERE

import string

from VcfOutputConfigTable import VcfOutputConfigTable
from oncotator.utils.ConfigUtils import ConfigUtils
from ConfigTableCreator import ConfigTableCreator


class VcfOutputConfigTableCreator(ConfigTableCreator):

    def __init__(self):
        pass

    def createConfigTableKeys(self, configParser, configTable):
        # Parse fields from FORMAT section of the config file
        """


        :param configParser:
        :param configTable:
        """
        # Parse fields from INFO section of the config file
        table = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(configParser, "INFO")
        for name, ID in table.items():
            configTable.addInfoFieldName(name, ID)

        table = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(configParser, "FORMAT")
        for name, ID in table.items():
            configTable.addFormatFieldName(name, ID)

        table = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(configParser, "OTHER")
        for name, ID in table.items():
            configTable.addOtherFieldName(name, ID)

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser, "INFO_DESCRIPTION")
        for name, description in table.items():
            configTable.addInfoFieldNameDescription(name, string.join(description, ","))

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser, "FORMAT_DESCRIPTION")
        for name, description in table.items():
            configTable.addFormatFieldNameDescription(name, string.join(description, ","))

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser, "FILTER_DESCRIPTION")
        for name, description in table.items():
            configTable.addFilterFieldNameDescription(name, string.join(description, ","))

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser, "SPLIT_TAGS")
        for fieldType, names in table.items():
            configTable.addFieldNamesToSplitSet(fieldType, names)

        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser, "NOT_SPLIT_TAGS")
        for fieldType, names in table.items():
            configTable.addFieldNamesToNotSplitSet(fieldType, names)

    def getConfigTable(self, configFilename, filename=None):
        """



        :param configFilename:
        :param filename:
        :return:
        """
        configParser = ConfigUtils.createConfigParser(configFilename, ignoreCase=False)
        configTable = VcfOutputConfigTable()
        self.createConfigTableKeys(configParser=configParser, configTable=configTable)

        return configTable