import logging
from shove import Shove
from oncotator.cache.Cache import Cache


class ShoveCache(Cache):
    """Expects a url in the form that shove requires.

    Maintains a cache of keys to speed performance."""

    def __init__(self, url_db, url_cache):
        self.db = Shove(url_db, url_cache, optimize=False)

        logging.getLogger(__name__).info("Loading Oncotator cache keys...")
        self.dbKeys = set(self.db.keys())

    def retrieve_from_cache(self, key):
        if key not in self.dbKeys:
            return None
        return self.db[key]

    def store_into_cache(self, key, value):
        self.db[key] = value
        self.dbKeys.add(key)

    def close_cache(self):
        self.db.close()