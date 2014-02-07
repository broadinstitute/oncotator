from DatasourceCreator import DatasourceCreator
from TabixIndexer import TabixIndexer
from ConfigParser import ConfigParser
import os
import pandas
import string
import numpy
import logging
from oncotator.index.InputMismatchException import InputMismatchException


class TabixIndexedTsvDatasourceCreator(DatasourceCreator):

    def __init__(self):
        self.columnDataTypes = {}

    def createDatasource(self, destDir, ds_file, index_column_names, column_names, configFilename, ds_type, ds_name,
                         ds_version, ds_match_mode, annotation_column_names, indexCols):
        # Create database
        baseDSFile = self._createDatabase(destDir, ds_file, index_column_names, annotation_column_names,
                                          column_names)
        # Create config file
        self._createConfigFile(configFilename, baseDSFile, ds_type, ds_name, ds_version, ds_match_mode, column_names,
                               annotation_column_names, indexCols)

    def _createDatabase(self, destDir, ds_file, index_column_names, annotation_column_names, column_names):
        index_column_names = index_column_names.split(",")
        annotation_column_names = annotation_column_names.split(",")
        column_names = column_names.split(",")
        index_columns = []
        for index_column_name in index_column_names:
            index_columns += [column_names.index(index_column_name)]

        column_names = set(column_names)
        annotation_column_names = set(annotation_column_names)

        # Read the column names and determine whether they exist or not in the file
        data = pandas.read_csv(filepath_or_buffer=ds_file, delimiter="\t", iterator=True, chunksize=1)
        for chunk in data:
            index = chunk.columns
            fieldNames = set(index)
            missingColumns = fieldNames.difference(column_names)
            if len(missingColumns) != 0:
                msg = "The input tsv, %s, is missing the following columns: %s." \
                    % (ds_file, string.join(missingColumns, ", "))
                raise InputMismatchException(msg)

            missingColumns = annotation_column_names.difference(fieldNames)
            if len(missingColumns) != 0:
                msg = "The input tsv, %s, is missing the following annotation columns: %s." \
                      % (ds_file, string.join(missingColumns))
                raise InputMismatchException(msg)
            break

        # Iterate through the file and determine column's data type
        data = pandas.read_csv(filepath_or_buffer=ds_file, delimiter="\t", iterator=True, chunksize=10000,
                               usecols=annotation_column_names, na_values=["", ".", "-"])
        for chunk in data:
            index = chunk.columns

            # Missing values default to float data type
            for idx in index:
                if numpy.issubdtype(chunk[idx].dtype, numpy.inexact):
                    if idx not in self.columnDataTypes or self.columnDataTypes[idx] not in ("String",):
                        self.columnDataTypes[idx] = "Float"
                elif numpy.issubdtype(chunk[idx].dtype, numpy.integer):
                    if idx not in self.columnDataTypes or self.columnDataTypes[idx] not in ("Float", "String",):
                        self.columnDataTypes[idx] = "Integer"
                elif numpy.issubdtype(chunk[idx].dtype, numpy.character):
                    self.columnDataTypes[idx] = "String"
                elif numpy.issubdtype(chunk[idx].dtype, numpy.bool_):
                    self.columnDataTypes[idx] = "Flag"
                elif numpy.issubdtype(chunk[idx].dtype, numpy.object_):
                    self.columnDataTypes[idx] = "String"
                else:
                    self.columnDataTypes[idx] = "String"

        tabixIndexedFile = TabixIndexer.index(destDir=destDir, inputFilename=ds_file, fileColumnNumList=index_columns,
                                              preset="tsv")
        baseDSFile = os.path.basename(tabixIndexedFile)
        logging.getLogger(__name__).info("%s file was created." % tabixIndexedFile)
        return baseDSFile

    def _createConfigFile(self, configFilename, baseDSFile, ds_type, ds_name, ds_version, ds_match_mode, column_names,
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
        config.optionxform = lambda option: option
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

        config.add_section("data_types")
        for column_name in self.columnDataTypes.keys():
            if column_name in annotation_column_names:
                config.set('data_types', column_name, self.columnDataTypes[column_name])

        config.write(filePtr)
        filePtr.close()
        logging.getLogger(__name__).info("%s file was created." % configFilename)
