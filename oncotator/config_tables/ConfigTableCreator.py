
from abc import abstractmethod


class ConfigTableCreator(object):

    @abstractmethod
    def getConfigTable(self, configFilename, filename):

        """

        :param configFilename:
        :param filename: depending on the data type, this is
        :raise:
        """
        raise NotImplementedError
