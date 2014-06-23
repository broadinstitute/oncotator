from ConfigParser import ConfigParser
import logging
import os
from oncotator.datasources.SnpOnlyLevelDbDatasource import SnpOnlyLevelDbDatasource
from oncotator.index.DatasourceCreator import DatasourceCreator
import leveldb
from oncotator.utils.GenericTsvReader import GenericTsvReader
import copy

class SnpOnlyLevelDbDatasourceCreator(DatasourceCreator):

    def __init__(self):
        pass

    def _createConfigFile(self, out_config_filename, baseDSFile, ds_name, ds_version, indexCols, annotation_columns):
        """

        :param out_config_filename: configuration filename to be created
        :param baseDSFile: data source base filename
        :param ds_type: data source type
        :param ds_name: data source title
        :param ds_version: data source version
        :param indexCols: list of string
        :param annotation_columns: List of str.  Should be populated.  empty list indicates that this datasource
            should not actually annotate anything.
        """
        config = ConfigParser()
        filePtr = open(out_config_filename, 'w')
        config.add_section("general")
        config.set("general", "version", ds_version)
        config.set("general", "title", ds_name)
        config.set("general", "type", "snp_leveldb")
        config.set("general", "src_file", baseDSFile)
        config.set("general", "index_column_names", ",".join(indexCols))
        config.set("general", "annotation_column_names", ",".join(annotation_columns))
        config.write(filePtr)
        filePtr.close()

    def createDatasource(self, destDir, ds_file, index_column_names, out_config_filename, ds_type, ds_name, ds_version,
                         ds_match_mode, annotation_column_names, indexCols):
        """


        :param destDir:
        :param ds_file:
        :param index_column_names:
        :param out_config_filename:
        :param ds_type:
        :param ds_name:
        :param ds_version:
        :param ds_match_mode:
        :param annotation_column_names:
        :param indexCols: list of the index columns.  Assumed to be five corresponding to chrom, start, end, ref, and alt.
        """
        output_filename = destDir + "/" + ds_name + ".leveldb"
        src_file = os.path.basename(output_filename)
        db = leveldb.LevelDB(output_filename, create_if_missing=True)

        tsv_file = ds_file
        tsv_reader = GenericTsvReader(tsv_file)

        logging.getLogger(__name__).info("Creating SNP LevelDB for the following index headers: " + str(index_column_names))
        logging.getLogger(__name__).info("Creating SNP LevelDB for the following data headers: " + str(annotation_column_names))

        # Create the config file
        self._createConfigFile(out_config_filename + "/" + ds_name + ".config", src_file, ds_name, ds_version, index_column_names, annotation_columns=annotation_column_names)

        batch = leveldb.WriteBatch()
        for i,line_dict in enumerate(tsv_reader):

            chrom = line_dict[index_column_names[0]]
            start = line_dict[index_column_names[1]]
            end = line_dict[index_column_names[2]]
            ref = line_dict[index_column_names[3]]
            alt = line_dict[index_column_names[4]]

            h = SnpOnlyLevelDbDatasource.generate_hash(chrom, start, end, ref, alt)
            if i % 5000 == 0:
                logging.getLogger(__name__).info("Rendering %d entries" % (i))

            line_list = [line_dict.get(k, "") for k in annotation_column_names]
            db.Put(h, ",".join(line_list))
        db.Write(batch, sync = True)



