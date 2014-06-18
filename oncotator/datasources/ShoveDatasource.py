from shove.core import Shove
from oncotator.datasources.Datasource import Datasource


class ShoveDatasource(Datasource):
    """Datasource backed by Shove.

    Initialization is done using a Shove URI.

    For example:
    file:///absolute/path/to/a/folder/
    leveldb:///absolute/path/to/a/folder/

    Values returned from a

    """

    def __init__(self, src_file, title='', version=None, shove_cache_url="simple://", max_entries=500):
        """

        :param src_file: Shove URI for store
        :param title:
        :param version:
        :param shove_cache_url: Shove URI for the datasource cache.  Defaults to NOT thread safe simple:// with 500 entries.
        One entry will be
        :param max_entries:
        """

        self._db_store = Shove(src_file, shove_cache_url, max_entries=max_entries)


