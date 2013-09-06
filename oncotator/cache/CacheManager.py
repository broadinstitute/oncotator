from oncotator.cache.CacheFactory import CacheFactory
from oncotator.utils.MutUtils import MutUtils


class CacheManager(object):
    def __init__(self):
        self.__cache = None
        self.__is_read_only = True
        self.__db_dir_key = None

    def initialize(self, url, db_dir_key, is_read_only=True):
        self.set_read_only(is_read_only)
        self.set_db_dir_key(db_dir_key)

        #TODO: Double check that multiple calls to initialize will not break or leave dangling resources.
        self.__cache = CacheFactory.createCache(url)

    def retrieve_cached_annotations(self, m):
        """
        :param m: mutation
        :return: list of Annotations, or None, if cache miss.
        """
        cache_key = MutUtils.create_variant_key_by_mutation(m, self.get_db_dir_key())
        return self.get_cache().retrieve_from_cache(cache_key)

    def _determine_annotations_to_cache(self, m):
        """Returns a list of the raw annotations that should be stored, if we are caching """
        return {annot:m.getAnnotation(annot) for annot in m.keys() if m.getAnnotation(annot).datasourceName not in ["INPUT", "__ATTR__", "OUTPUT", "MANUAL", "DEFAULT"]}

    def store_annotations_in_cache(self, m):

        if self.is_read_only():
            return

        cache_key = MutUtils.create_variant_key_by_mutation(m, self.get_db_dir_key())
        annotations = self._determine_annotations_to_cache(m)
        self._store_basic_annotations_in_cache(cache_key, annotations)

    def close_cache(self):
        self.__cache.close_cache()

    def _store_basic_annotations_in_cache(self, key, annotations):
        self.get_cache().store_into_cache(key, annotations)

    def get_cache(self):
        return self.__cache

    def set_cache(self, value):
        self.__cache = value

    def get_db_dir_key(self):
        return self.__db_dir_key

    def set_db_dir_key(self, value):
        self.__db_dir_key = value

    def set_read_only(self, value):
        self.__is_read_only = value

    def is_read_only(self):
        return self.__is_read_only