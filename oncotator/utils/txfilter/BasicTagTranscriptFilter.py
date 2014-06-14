# LICENSE_GOES_HERE

class BasicTagTranscriptFilter(object):
    """ Class for transcript filter that only passes transcripts that have a tag of basic.

     This is used by Basic GENCODE transcript providers.
     """

    def filter(self, txs):
        """ GENCODE transcripts contain tags that are useful for QC.  If tags are present, this method will remove
        transcripts with poor QC.

        For now, just accept "basic"  Prune any transcripts that have no tags.

        :param txs: transcripts to possibly prune
        :return: a list of same or shorter
        """
        return [tx for tx in txs if ('tag' in tx.get_other_attributes().keys()) and ('basic' in tx.get_other_attributes()['tag'])]
