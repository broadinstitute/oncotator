from shove.core import Shove
from oncotator.datasources.Datasource import Datasource
import leveldb


class ShoveDatasource(Datasource):
    """Datasource backed by Shove.

    Initialization is done using a Shove URI.

    For example:
    file:///absolute/path/to/a/folder/
    leveldb:///absolute/path/to/a/folder/

    Values returned from a ShoveDatasrource will be

    Index columns must match exactly in order to receive a value from the ShoveDatasource

    """

    def __init__(self, src_file, title, version, index_cols, annotation_columns, shove_cache_url="simple://", max_entries=500):
        """

        :param src_file: Shove URI for store
        :param title:
        :param version:
        :param shove_cache_url: Shove URI for the datasource cache.  Defaults to NOT thread safe simple:// with 500 entries.
        One entry will be
        :param max_entries:
        """
        super(ShoveDatasource, self).__init__(src_file, title=title, version=version)
        # self._db_store = Shove(src_file, shove_cache_url, max_entries=max_entries)
        self._db_store = leveldb.LevelDB(src_file)
        self._annotation_columns = annotation_columns

    def annotate_mutation(self, mutation):
        """ Mutations are annotated only with exact matches of chr, start, end, ref, and alt.
        """

        # create hash for this mutation
        h = ShoveDatasource.generate_hash(mutation)

        # extract value for this hash from the db
        annotations_list = []
        try:
            annotations_list = self._db_store.Get(h)
        except KeyError:
            # do nothing
            pass

        # Annotate
        for i,col in enumerate(self._annotation_columns):
            if len(annotations_list) <= i:
                mutation.createAnnotation(col, "", self.title)
            else:
                mutation.createAnnotation(col, annotations_list[i], self.title)

        return mutation



    @staticmethod
    def generate_hash(m):
        return ShoveDatasource.generate_hash(m.chr, m.start, m.end, m.ref_allele, m.alt_allele)

    @staticmethod
    def generate_hash(chrom, start,end,ref,alt):
        return "%s_%s_%s_%s_%s" % (chrom, start, end, ref, alt)

