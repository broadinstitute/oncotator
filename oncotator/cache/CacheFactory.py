import logging
from oncotator.cache.DummyCache import DummyCache
from oncotator.cache.MemcacheCache import MemcacheCache
from oncotator.cache.ShoveCache import ShoveCache


class CacheFactory(object):
    """Currently, only supports Shove caches. """

    @staticmethod
    def createCache(url):
        """
        :param url:
        :return: instance of Cache
        """

        if url is None:
            return DummyCache()

        cache_url = None
        if url.startswith("file"):

            # We should set up in memory cache as well.
            cache_url = "memory://"

        if url.startswith("memcache"):
            logging.getLogger(__name__).info("Using memcache: " + str(url) )
            return MemcacheCache(url)

        logging.getLogger(__name__).info("Creating shove cache: " + str(url) + ", " + str(cache_url))
        return ShoveCache(url, cache_url)