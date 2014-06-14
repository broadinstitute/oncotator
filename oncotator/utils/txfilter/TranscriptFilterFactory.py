# LICENSE_GOES_HERE

import logging
from oncotator.utils.txfilter.BasicTagTranscriptFilter import BasicTagTranscriptFilter
from oncotator.utils.txfilter.DummyTranscriptFilter import DummyTranscriptFilter


class TranscriptFilterFactory(object):

    TRANSCRIPT_FILTER_DICT = {"basic": BasicTagTranscriptFilter, "dummy": DummyTranscriptFilter}

    @staticmethod
    def create_instance(transcript_filter_type):
        """

        :param transcript_filter_type: (str)
        :return: a transcript filter instance of the proper type.
        """
        if transcript_filter_type not in TranscriptFilterFactory.TRANSCRIPT_FILTER_DICT.keys():
            logging.getLogger(__name__).warn(transcript_filter_type + " is not a valid type: " + str(TranscriptFilterFactory.TRANSCRIPT_FILTER_DICT.keys()) + " .... using a dummy filter that does nothing.")
            transcript_filter_type = "dummy"

        return TranscriptFilterFactory.TRANSCRIPT_FILTER_DICT[transcript_filter_type]()


