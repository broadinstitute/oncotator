from oncotator.cache.Cache import Cache
import memcache

class MemcacheCache(Cache):
    """
    """
    def __init__(self, url):
        ip_address = url.split("://")[1]
        self.mc = memcache.Client([ip_address], debug=0)

    def retrieve_from_cache(self, key):
        """Given a key, retrieve from the cache."""
        obj = self.mc.get(key)
        if not obj:
            return None
        else:
            return obj

    def store_into_cache(self, key, value):
        self.mc.set(key, value)


    def close_cache(self):
        """Do nothing """
        pass