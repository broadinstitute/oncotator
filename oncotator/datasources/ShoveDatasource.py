from shove.core import Shove
from oncotator.datasources.Datasource import Datasource


class ShoveDatasource(Datasource):
    """Datasource backed by Shove.

    Initialization is done using a Shove URI.

    For example:
    file:///absolute/path/to/a/folder/
    leveldb:///absolute/path/to/a/folder/

    Values returned from a ShoveDatasrource will be

    Index columns must match exactly in order to receive a value from the ShoveDatasource

    """

    def __init__(self, src_file, title, version, index_cols, shove_cache_url="simple://", max_entries=500):
        """

        :param src_file: Shove URI for store
        :param title:
        :param version:
        :param shove_cache_url: Shove URI for the datasource cache.  Defaults to NOT thread safe simple:// with 500 entries.
        One entry will be
        :param max_entries:
        """
        super(ShoveDatasource, self).__init__(src_file, title=title, version=version)
        self._db_store = Shove(src_file, shove_cache_url, max_entries=max_entries)

    def annotate_mutation(self, mutation):
        """ Mutations are annotated only with exact matches of chr, start, end, ref, and alt.
        """

        # create hash for this mutation

        # extract value for this hash from the db

        # Annotate

    @staticmethod
    def generate_hash(m):
        return "%s_%s_%s_%s_%s" % (m.chr, m.start, m.end, m.ref_allele, m.alt_allele)


