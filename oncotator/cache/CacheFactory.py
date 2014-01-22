import logging
from oncotator.cache.DummyCache import DummyCache
from oncotator.cache.MemcacheCache import MemcacheCache
from oncotator.cache.ShoveCache import ShoveCache


class CacheFactory(object):
    """Currently, only supports Shove cache and a custom interface to memcache. """

    @staticmethod
    def createCache(url, is_thread_safe=False):
        """

        :param is_thread_safe: If this is True, uses a slower, but thread safe version of the Cache, if available.
        :param url:
        :return: instance of Cache
        """

        if url is None or url == "":
            return DummyCache()

        cache_url = None
        if url.startswith("file"):
            logging.getLogger(__name__).info("Creating file-based cache: " + str(url) )
            if is_thread_safe:
                cache_url = "memory://"
            else:
                # We should set up in-memory cache as well.
                logging.getLogger(__name__).info("Cache being set up as faster, but NOT thread safe." )
                cache_url = "simple://"

        if url.startswith("memcache"):
            logging.getLogger(__name__).info("Creating memcache cache: " + str(url) )
            return MemcacheCache(url)

        logging.getLogger(__name__).info("Creating shove cache: " + str(url) + ", " + str(cache_url))
        return ShoveCache(url, cache_url)