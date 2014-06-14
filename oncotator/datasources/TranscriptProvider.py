# LICENSE_GOES_HERE

import abc


class TranscriptProvider(object):

    TX_MODE_CANONICAL = "CANONICAL"
    TX_MODE_BEST_EFFECT = "EFFECT"
    TX_MODE_CHOICES = [TX_MODE_CANONICAL, TX_MODE_BEST_EFFECT] # Simply list the ones above here.

    @abc.abstractmethod
    def getTranscriptDict(self):
        """ Return a dict containing all transcripts where key is the transcript ID.
        """
        return

    @abc.abstractmethod
    def get_transcript(self, tx_id):
        """ Return a Transcript instance for the given transcript ID.
        """
        return

    @abc.abstractmethod
    def get_transcripts_by_pos(self, chr, start, end):
        """ Return Transcript instances for the given chromosome start and end.
        """
        return

    @abc.abstractmethod
    def retrieveExons(self, gene, padding=10, isCodingOnly=False):
        """Return a list of (chr, start, end) tuples for each exon in each transcript"""
        return

    @abc.abstractmethod
    def set_tx_mode(self, tx_mode):
        # TODO: Throw exception if not in TX_MODE_CHOICES
        return

    @abc.abstractmethod
    def get_tx_mode(self):
        return