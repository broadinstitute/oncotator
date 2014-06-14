# LICENSE_GOES_HERE
import abc


class Cache(object):
    """ Base class for Caches
    """

    @abc.abstractmethod
    def retrieve_from_cache(self, key):
        """Given a key, retrieve from the cache."""
        pass

    @abc.abstractmethod
    def store_into_cache(self, key, value):
        pass

    @abc.abstractmethod
    def close_cache(self):
        pass