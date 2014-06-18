"""
By downloading the PROGRAM you agree to the following terms of use:

BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY

This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").

WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and

WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.

NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:

1. DEFINITIONS
1.1 "PROGRAM" shall mean the object code and source code known as Oncotator 1.0 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.  BROAD acknowledges that the PROGRAM employs one or more public domain code(s) that are freely available for public use.

2. LICENSE
2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.  LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.

3. OWNERSHIP OF INTELLECTUAL PROPERTY
LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.

Copyright 2014 Broad Institute, Inc.
Notice of attribution:  The Oncotator 1.0 program was made available through the generosity of the Broad Institute, Inc.

LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.

4. INDEMNIFICATION
LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorney fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.

5. NO REPRESENTATIONS OR WARRANTIES
THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.

6. ASSIGNMENT
This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.

7. MISCELLANEOUS
7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
"""
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

    def _determine_column_names(self, ds_file):
        column_names = []
        data = pandas.read_csv(filepath_or_buffer=ds_file, delimiter="\t", iterator=True, chunksize=1)
        for chunk in data:
            column_names = list(chunk.columns)
            break
        return column_names

    def createDatasource(self, destDir, ds_file, index_column_names, configFilename, ds_type, ds_name,
                         ds_version, ds_match_mode, annotation_column_names, indexCols):
        # Create database
        column_names = self._determine_column_names(ds_file)
        if not index_column_names:
            index_column_names = column_names
        else:
            index_column_names = index_column_names.split(",")

        if not annotation_column_names:  # no annotation column names were specified; default to index columns
            annotation_column_names = [index_column_name for index_column_name in index_column_names
                                       if index_column_name not in index_column_names]
        else:
            annotation_column_names = annotation_column_names.split(",")

        if not all([index_column_name in column_names for index_column_name in index_column_names]):
            raise ValueError("index_column names must be a subset of column names.")
        if not all([annotation_column_name in column_names for annotation_column_name in annotation_column_names]):
            raise ValueError("annotation_columm names must be a subset of column names.")

        baseDSFile = self._createDatabase(destDir, ds_file, ds_match_mode, index_column_names, annotation_column_names,
                                          column_names)
        # Create config file
        column_names = string.join(column_names, ",")
        annotation_column_names = string.join(annotation_column_names, ",")
        self._createConfigFile(configFilename, baseDSFile, ds_type, ds_name, ds_version, ds_match_mode, column_names,
                               annotation_column_names, indexCols)

    def _createDatabase(self, destDir, ds_file, ds_match_mode, index_column_names, annotation_column_names,
                        column_names):
        index_columns = []
        for index_column_name in index_column_names:  # ensures that all index columns are in the column names list
            index_columns += [column_names.index(index_column_name)]

        if len(index_columns) != 3 and len(index_columns) != 5:
            raise ValueError("Wrong number of index columns.  Must be a comma separated list of length 3 or 5.")

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
                if ds_match_mode != "exact" or len(index_column_names) != 3:
                    if numpy.issubdtype(chunk[idx].dtype, numpy.inexact):
                        if idx not in self.columnDataTypes or self.columnDataTypes[idx] not in ("String",):
                            self.columnDataTypes[idx] = "Float"
                    elif numpy.issubdtype(chunk[idx].dtype, numpy.integer):
                        if idx not in self.columnDataTypes or self.columnDataTypes[idx] not in ("Float", "String",):
                            self.columnDataTypes[idx] = "Integer"
                    elif numpy.issubdtype(chunk[idx].dtype, numpy.bool_):
                        self.columnDataTypes[idx] = "Flag"
                    else:
                        self.columnDataTypes[idx] = "String"
                else:
                    self.columnDataTypes[idx] = "String"

            if ds_match_mode != "exact" or len(index_column_names) != 3:
                break

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
