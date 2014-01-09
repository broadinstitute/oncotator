from DatasourceCreator import DatasourceCreator
from TabixIndexer import TabixIndexer
from ConfigParser import ConfigParser
import os


class TabixIndexedTsvDatasourceCreator(DatasourceCreator):

    def __init__(self):
        pass

    def createDatasource(self, destDir, ds_file, index_column_names, column_names):
        index_column_names = index_column_names.split(",")
        column_names = column_names.split(",")
        index_columns = []
        for index_column_name in index_column_names:
            index_columns += [column_names.index(index_column_name)]

        tabixIndexedFile = TabixIndexer.index(destDir=destDir, inputFilename=ds_file, fileColumnNumList=index_columns,
                                              preset="tsv")
        baseDSFile = os.path.basename(tabixIndexedFile)

        return baseDSFile

    def createConfigFile(self, configFilename, baseDSFile, ds_type, ds_name, ds_version, ds_match_mode, column_names,
                         annotation_column_names, indexCols):
        """


        :param configFilename: configuration filename
        :param baseDSFile: data source base filename
        :param ds_type: data source type
        :param ds_name: data source title
        :param ds_version: data source version
        :param column_names: column names in the input data source file
        :param annotation_column_names: column names whose values are used for annotation
        :param indexCols: named tuple consisting of index column type and corresponding column names
        :param ds_match_mode: describes how to annotate mutations from an indexed tsv or indexed vcf datasources
        """
        config = ConfigParser()
        filePtr = open(configFilename, 'w')
        config.add_section("general")
        config.set("general", "version", ds_version)
        config.set("general", "title", ds_name)
        config.set("general", "type", ds_type)
        config.set("general", "src_file", baseDSFile)
        config.set("general", "column_names", column_names)
        config.set("general", "annotation_column_names", annotation_column_names)
        config.set("general", indexCols.type, indexCols.names)
        config.set("general", "match_mode", ds_match_mode)
        config.write(filePtr)
        filePtr.close()