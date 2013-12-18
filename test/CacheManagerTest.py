import unittest
from oncotator.MutationData import MutationData
from oncotator.cache.CacheManager import CacheManager

class CacheManagerTest(unittest.TestCase):
    _multiprocess_can_split_ = True
    def test_cached_annots(self):
        """Test to make sure that we are not storing annotations that should not be cached.  Also, tests a simple store and retrieve."""
        cache_file = "out/shove.managertest.annots.cache"
        cm = CacheManager()
        fake_db_dir_key = "blah"
        cm.initialize("file://" + cache_file, fake_db_dir_key, is_read_only=False)
        m = MutationData()
        m.createAnnotation("blah1", "val1", annotationSource="INPUT")
        m.createAnnotation("blah2", "val5", annotationSource="some_datasource")
        cm.store_annotations_in_cache(m)
        annots = cm.retrieve_cached_annotations(m)

        self.assertTrue(len(annots.keys()) == 1)
        self.assertTrue(annots["blah2"].getValue() == "val5")

    def test_cached_annots_dummy_cache(self):
        """Test dummy cache.  Also, tests a simple store and retrieve, which should be None."""
        cm = CacheManager()
        fake_db_dir_key = "blah"
        cm.initialize(None, fake_db_dir_key, is_read_only=False)
        m = MutationData()
        m.createAnnotation("blah1", "val1", annotationSource="INPUT")
        m.createAnnotation("blah2", "val5", annotationSource="some_datasource")
        cm.store_annotations_in_cache(m)
        annots = cm.retrieve_cached_annotations(m)
        self.assertTrue(annots is None)



if __name__ == '__main__':
    unittest.main()
