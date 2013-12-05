import abc

__author__ = 'lichtens'

class TranscriptFilter(object):
    """Base class for transcript filters """

    @abc.abstractmethod
    def filter(self, txs):
        """

        :param txs: list of Transcripts
        :return:
        """
        raise NotImplementedError("This is an abstract version")