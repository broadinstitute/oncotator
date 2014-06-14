# LICENSE_GOES_HERE
import os
import string
from oncotator.utils.Hasher import Hasher
import collections
from oncotator.index.DatasourceCreatorFactory import DatasourceBuilderFactory


class DatasourceInstallUtils(object):

    indexCols = collections.namedtuple("indexCols", ["type", "names"])

    @staticmethod
    def getIndexCols(dsType, index_columns):
        if dsType.endswith("gp_tsv"):
            return DatasourceInstallUtils.indexCols("genomic_position_cols", index_columns)
        elif dsType.endswith("gene_tsv"):
            return DatasourceInstallUtils.indexCols("gene_col", index_columns)
        elif dsType.endswith("transcript_tsv"):
            return DatasourceInstallUtils.indexCols("transcript_col", index_columns)
        elif dsType.endswith("gpp_tsv"):
            return DatasourceInstallUtils.indexCols("gene_protein_position_cols", index_columns)
        elif dsType.endswith("indexed_tsv"):
            return DatasourceInstallUtils.indexCols("index_column_names", index_columns)
        return None

    @staticmethod
    def create_datasource_md5_file(datasource_dir):
        """datasource_dir should be the /db_dir/ds_name/genome_build.
        For example,
        create_datasource_md5_file("/home/user/my_db_dir/gaf/hg19")
        """
        if datasource_dir.endswith('/'):
            datasource_dir = datasource_dir[:-1]
        md5_filename = os.path.abspath(datasource_dir) + ".md5"
        print("md5 being written to: " + os.path.abspath(md5_filename))
        hasher = Hasher()
        hashcode = hasher.create_hashcode_for_dir(datasource_dir)
        fp = file(md5_filename, "w")
        fp.write(hashcode)
        fp.close()

    @staticmethod
    def create_datasource(destDir, ds_file, ds_foldername, ds_name, ds_type, ds_version, index_columns=None,
                          ds_annotation_columns=None, ds_match_mode="exact"):
        """

        :param ds_match_mode: describes how to annotate mutations from an indexed tsv or indexed vcf datasources
        :param destDir: temporary destination directory (tmpdir/ds_foldername/genome_build)
        :param ds_file: data source filename
        :param ds_foldername:
        :param ds_name:
        :param ds_type: data source type (indexed_vcf, indexed_tsv, etc.)
        :param ds_version: data source version
        :param index_columns:
        :param ds_annotation_columns: if data source is of type indexed tsv,
        """
        index_columns = [] if index_columns is None else index_columns
        datasourceBuilder = DatasourceBuilderFactory.getDatasourceCreatorInstance(ds_type)
        configFilename = os.path.join(*[destDir, string.join([ds_foldername, "config"], ".")])
        datasourceBuilder.createDatasource(destDir=destDir, ds_file=ds_file, index_column_names=index_columns,
                                           configFilename=configFilename, ds_type=ds_type, ds_name=ds_name,
                                           ds_version=ds_version, ds_match_mode=ds_match_mode,
                                           annotation_column_names=ds_annotation_columns,
                                           indexCols=DatasourceInstallUtils.getIndexCols(ds_type, index_columns))

        DatasourceInstallUtils.create_datasource_md5_file(destDir)