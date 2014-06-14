# LICENSE_GOES_HERE

class DummyTranscriptFilter(object):
    """ Class for transcript filter that does not actually filter anything """

    def filter(self, txs):
        """
        Does nothing

        :param txs: list of Transcripts
        :return: txs
        """
        return txs