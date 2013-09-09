from oncotator.cache.ShoveCache import ShoveCache

__author__ = 'lichtens'

import unittest


class ShoveCacheTest(unittest.TestCase):
    def test_store(self):
        cache_file = "out/shove.cache"
        cache = ShoveCache("file://" + cache_file, "memory://")
        key1 = "dummy1"
        value1 = ["1", "5"]
        key2 = "dummy2"
        value2 = {'boo':'1', 'boo2':'2'}
        cache.store_into_cache(key1, value1)
        cache.store_into_cache(key2, value2)

        retrieved_val1 = cache.retrieve_from_cache("dummy1")
        self.assertTrue(retrieved_val1[1] == "5")
        retrieved_val2 = cache.retrieve_from_cache("dummy2")
        self.assertTrue(retrieved_val2["boo2"] == "2")


if __name__ == '__main__':
    unittest.main()
