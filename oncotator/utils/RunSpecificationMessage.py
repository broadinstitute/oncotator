import logging


class RunSpecificationMessage(object):
    """A simple entity class for relaying messages regarding a run specification.

    Though this class uses the logging levels for its own level field, a caller is able to override.
    """
    def __init__(self, level, message):
        """
        Constructor using logging level and a message.

        :param logging_level level: the level of this message (e.g. logging.WARN)  Should use the logging enumeration, but this is not enforced.
        :param str message:
        """
        self._level = level
        self._message = message

    @property
    def level(self):
        return self._level

    @level.setter
    def level(self, value):
        self._level = value

    @property
    def message(self):
        return self._message

    @message.setter
    def message(self, value):
        self._message = value