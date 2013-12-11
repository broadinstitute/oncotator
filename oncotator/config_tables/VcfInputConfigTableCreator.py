import vcf
from VcfInputConfigTable import VcfInputConfigTable
from ConfigTableCreator import ConfigTableCreator
from oncotator.utils.ConfigUtils import ConfigUtils


class VcfInputConfigTableCreator(ConfigTableCreator):

    def __init__(self):
        pass

    def createConfigTableKeys(self, configParser, configTable):
        # Parse fields from FORMAT section of the config file
        """


        :param configParser:
        :param configTable:
        """
        table = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(configParser, "INFO")
        for ID, name in table.items():
            configTable.addInfoFieldID(ID, name)

        # Parse fields from FORMAT section of the config file
        table = ConfigUtils.buildReverseAlternativeDictionaryFromConfig(configParser, "FORMAT")
        for ID, name in table.items():
            configTable.addFormatFieldID(ID, name)

        # Parse fields from NOT_SPLIT_TAGS section of the config file
        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser, "NOT_SPLIT_TAGS")
        for fieldType, IDs in table.items():
            configTable.addFieldIDsToNotSplitSet(fieldType, IDs)

        # Parse fields from SPLIT_TAGS section of the config file
        table = ConfigUtils.buildAlternateKeyDictionaryFromConfig(configParser, "SPLIT_TAGS")
        for fieldType, IDs in table.items():
            configTable.addFieldIDsToSplitSet(fieldType, IDs)

    def createConfigTable(self, vcfReader, configTable):
        """


        :param vcfReader:
        :param configTable:
        """
        # Remove INFO field IDs that are present in the config file but not in the input file
        for ID in configTable.getInfoFieldIDs():
            if ID not in vcfReader.infos.keys():
                configTable.removeInfoFieldID(ID)

        # Insert INFO fields IDs that are missing in the config file but are present in the input file
        for ID in vcfReader.infos.keys():
            if ID not in configTable.getInfoFieldIDs():
                name = ID
                configTable.addInfoFieldID(ID, name)

        # Remove FORMAT field IDs that are present in the config file but not in the input file
        for ID in configTable.getFormatFieldIDs():
            if ID not in vcfReader.formats.keys():
                configTable.removeFormatFieldID(ID)

        # Insert FORMAT fields IDs that are missing in the config file but are present in the input file
        for ID in vcfReader.formats.keys():
            if ID not in configTable.getFormatFieldIDs():
                name = ID
                configTable.addFormatFieldID(ID, name)

    def getConfigTable(self, configFilename, filename):
        """


        :return:
        """
        configParser = ConfigUtils.createConfigParser(configFilename, ignoreCase=False)
        configTable = VcfInputConfigTable()
        vcfReader = vcf.Reader(filename=filename, strict_whitespace=True)

        self.createConfigTableKeys(configParser=configParser, configTable=configTable)
        self.createConfigTable(vcfReader=vcfReader, configTable=configTable)

        return configTable