from oncotator.cache.Cache import Cache
import memcache

class MemcacheCache(Cache):
    """

    Please note that each mutation stored in the memcache uses about 12k due to pickling, if set to default pickling protocol.
    With the pickling protocol set here, each entry uses about 8k.

    """
    def __init__(self, url):
        ip_address = url.split("://")[1]
        self.mc = memcache.Client([ip_address], debug=0, pickleProtocol=2)

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