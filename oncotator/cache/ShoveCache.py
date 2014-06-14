# LICENSE_GOES_HERE

import logging
from shove import Shove
from oncotator.cache.Cache import Cache


class ShoveCache(Cache):
    """Expects a url in the form that shove requires.

    Maintains a cache of keys to speed performance."""

    def __init__(self, url_db, url_cache):
        self.db = Shove(url_db, url_cache, optimize=False, max_entries=2000)

    def retrieve_from_cache(self, key):
        try:
            val = self.db[key]
            return val
        except KeyError:
            return None

    def store_into_cache(self, key, value):
        self.db[key] = value

    def close_cache(self):
        self.db.close()