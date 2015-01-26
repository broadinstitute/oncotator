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

import logging
import os
from oncotator.MockExceptionThrowingDatasource import MockExceptionThrowingDatasource
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.utils.RunSpecification import RunSpecification
from utils.ConfigUtils import ConfigUtils
from oncotator.datasources.Gaf import Gaf
from oncotator.datasources.ReferenceDatasource import ReferenceDatasource
from oncotator.datasources.dbSNP import dbSNP
from oncotator.datasources.Cosmic import Cosmic
from oncotator.datasources.GenericGeneDatasource import GenericGeneDatasource
from oncotator.datasources.GenericTranscriptDatasource import GenericTranscriptDatasource
from oncotator.datasources.GenericVariantClassificationDatasource import GenericVariantClassificationDatasource
from oncotator.datasources.GenericGenomicPositionDatasource import GenericGenomicPositionDatasource
from oncotator.datasources.GenericGeneProteinPositionDatasource import GenericGeneProteinPositionDatasource
from oncotator.datasources.PositionTransformingDatasource import PositionTransformingDatasource
from oncotator.datasources.TranscriptToUniProtProteinPositionTransformingDatasource import TranscriptToUniProtProteinPositionTransformingDatasource
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.datasources.GenericGenomicMutationDatasource import GenericGenomicMutationDatasource
from oncotator.datasources.TabixIndexedTsvDatasource import IndexedTsvDatasource
from oncotator.datasources.TabixIndexedVcfDatasource import IndexedVcfDatasource
from oncotator.datasources.ChangeTransformingDatasource import ChangeTransformingDatasource
from oncotator.datasources.BigWigDatasource import BigWigDatasource
from utils.MultiprocessingUtils import LoggingPool
from oncotator import NGSLIB_INSTALLED

#TODO:  futures (python lib -- 2.7 backport exists on pypi) is more flexible and less error prone


def createDatasource(t):
    """ Create a datasource given a tuple (configFilename, leafDir).  This method should not be used and is only for a workaround to enable multiprocessing. """
    return DatasourceFactory.createDatasourceGivenTuple(t)


class DatasourceFactory(object):
    """
    Static class that creates instances of datasources.
    #TODO: Rename as a factory, since DatasourceCreator is a Datasource Factory
    
    """

    def __init__(self,params):
        """
        Constructor
        """
        raise NotImplementedError('This class is static, call methods as DatasourceCreator.<methodName>()')
    
    @staticmethod
    def createDatasource(configFilename, leafDir):
        configParser = ConfigUtils.createConfigParser(configFilename)
        return DatasourceFactory.createDatasourceFromConfigParser(configParser, leafDir)
    
    @staticmethod
    def createDatasourceGivenTuple(configTuple):
        """ Calls createDatasourceFromConfigParser with a two-entry tuple using 
            exact same arguments. """
        return DatasourceFactory.createDatasource(configTuple[0], configTuple[1])

    @staticmethod
    def _retrieve_hash_code(leafDir):
        hashcode = ""
        md5_filename = os.path.dirname(leafDir) + ".md5"
        if os.path.exists(md5_filename):
            logging.getLogger(__name__).debug("md5 found for " + leafDir)
            md5_fp = file(md5_filename, 'r')
            hashcode = md5_fp.read()
            md5_fp.close()
        return hashcode

    @staticmethod
    def _log_missing_column_name_msg(colnames, indexColnames):
        for colname in indexColnames:
            if colname not in colnames:
                msg = "%s is missing from column name list." % colname
                logging.getLogger(__name__).warn(msg)

    @staticmethod
    def createDatasourceFromConfigParser(configParser, leafDir):
        """
        configParser -- config parser instance from the config file in the leafdir. For information on config file format/conventions see (TODO)
        
        leafDir -- contains the file and necessary files (post indexing and install steps) to instantiate a datasource.

        """
        result = None
        # Determine the type
        dsType = configParser.get("general", "type")
        
        # TODO: Replace these if statements with something a bit more robust, such as a proper dependency injection framework
        filePrefix = leafDir + "/"
        if dsType == "gaf":
            gaf_fname = filePrefix + configParser.get('general', 'gaf_fname')
            gaf_transcript_sequences_fname = filePrefix + configParser.get('general', 'gaf_transcript_seqs_fname')
            result = Gaf(gaf_fname, gaf_transcript_sequences_fname, title=configParser.get("general", "title"), version=configParser.get("general", "version"), protocol=configParser.get("general", "protocol"))
        elif dsType == "dbsnp":
            result = dbSNP(filePrefix + configParser.get('general', 'src_file'), title=configParser.get('general', 'title'), version=configParser.get('general', 'version'))
        elif dsType == "ensembl":
            result = EnsemblTranscriptDatasource(filePrefix + configParser.get('general', 'src_file'),
                                                 title=configParser.get('general', 'title'),
                                                 version=configParser.get('general', 'version'),
                                                 tx_filter=configParser.get('general', 'transcript_filter'))
        elif dsType == "cosmic":
            result = Cosmic(src_file=filePrefix + configParser.get('general', 'src_file'), version=configParser.get('general', 'version'), gpp_tabix_file=filePrefix + configParser.get('general', 'gpp_src_file'))
        elif dsType == 'ref':
            if configParser.has_option('general', 'windowSizeRef'):
                window_size = configParser.get('general', 'windowSizeRef')
            else:
                window_size = 10
            result = ReferenceDatasource(filePrefix, title=configParser.get("general", "title"), version=configParser.get('general', 'version'), windowSizeRef=window_size)
        elif dsType == 'gene_tsv':
            result = GenericGeneDatasource(src_file=filePrefix + configParser.get('general', 'src_file'), title=configParser.get("general", "title"), version=configParser.get('general', 'version'), geneColumnName=configParser.get('general', 'gene_col'))
        elif dsType == 'transcript_tsv':
            result = GenericTranscriptDatasource(src_file=filePrefix + configParser.get('general', 'src_file'), title=configParser.get("general", "title"), version=configParser.get('general', 'version'), geneColumnName=configParser.get('general', 'transcript_col'))
        elif dsType == 'vc_tsv':
            result = GenericVariantClassificationDatasource(src_file=filePrefix + configParser.get('general', 'src_file'), title=configParser.get("general", "title"), version=configParser.get('general', 'version'), geneColumnName=configParser.get('general', 'vc_col'))
        elif dsType == 'gp_tsv':
            result = GenericGenomicPositionDatasource(src_file=filePrefix + configParser.get('general', 'src_file'), title=configParser.get("general", "title"), version=configParser.get('general', 'version'), gpColumnNames=configParser.get('general', 'genomic_position_cols'))
        elif dsType == 'gm_tsv':
            result = GenericGenomicMutationDatasource(src_file=filePrefix + configParser.get('general', 'src_file'), title=configParser.get("general", "title"), version=configParser.get('general', 'version'), gpColumnNames=configParser.get('general', 'genomic_position_cols'))
        elif dsType == 'gm_tsv_reverse_complement':
            result = GenericGenomicMutationDatasource(src_file=filePrefix + configParser.get('general', 'src_file'), title=configParser.get("general", "title"), version=configParser.get('general', 'version'), gpColumnNames=configParser.get('general', 'genomic_position_cols'), use_complementary_strand_alleles_for_negative_strand_transcripts=True)
        elif dsType == 'gpp_tsv':
            result = GenericGeneProteinPositionDatasource(src_file=filePrefix + configParser.get('general', 'src_file'),title=configParser.get("general", "title"), version=configParser.get('general', 'version'), gpColumnNames=configParser.get('general', 'gene_protein_position_cols'))
        elif dsType == "transcript_to_uniprot_aa":
            result = TranscriptToUniProtProteinPositionTransformingDatasource(title=configParser.get("general", "title"),
                                                                              version=configParser.get('general', 'version'),
                                                                              src_file="file://" + filePrefix + configParser.get('general', 'src_file'), # three slashes for sqlite
                                                                              inputPositionAnnotationName=configParser.get('general', 'inputPositionAnnotationName'),
                                                                              outputPositionAnnotationName=configParser.get('general','outputPositionAnnotationName'))
        
        elif dsType == "mock_exception":
            result = MockExceptionThrowingDatasource(title=configParser.get("general", "title"), version=configParser.get('general', 'version'))
        elif dsType == "indexed_vcf":
            result = IndexedVcfDatasource(src_file=filePrefix + configParser.get('general', 'src_file'),
                                           title=configParser.get("general", "title"),
                                           version=configParser.get('general', 'version'),
                                           match_mode=configParser.get('general', 'match_mode'))
        elif dsType == "indexed_tsv":
            columnNames = configParser.get("general", "column_names")
            columnNames = columnNames.split(",")

            annotationColumnNames = configParser.get("general", "annotation_column_names")
            annotationColumnNames = annotationColumnNames.split(",")

            indexColumnNames = configParser.get("general", "index_column_names")
            indexColumnNames = indexColumnNames.split(",")

            DatasourceFactory._log_missing_column_name_msg(columnNames, annotationColumnNames)

            columnDataTypes = dict()
            for columnName in annotationColumnNames:
                if columnName.strip() == "":
                    continue
                columnDataTypes[columnName] = configParser.get("data_types", columnName)

            result = IndexedTsvDatasource(src_file=filePrefix + configParser.get("general", "src_file"),
                                           title=configParser.get("general", "title"),
                                           version=configParser.get("general", "version"),
                                           colNames=columnNames,
                                           annotationColNames=annotationColumnNames,
                                           indexColNames=indexColumnNames,
                                           match_mode=configParser.get("general", "match_mode"),
                                           colDataTypes=columnDataTypes)

        
        elif dsType == 'bigwig':
            if not NGSLIB_INSTALLED:
                raise RuntimeError("Bigwig datasource found in db-dir but ngslib library not installed.")
            result = BigWigDatasource(src_file=filePrefix + configParser.get('general', 'src_file'), title=configParser.get("general", "title"), version=configParser.get('general', 'version'))
        else:
            raise RuntimeError('Unknown datasource type: %s' % dsType)


        hashcode = DatasourceFactory._retrieve_hash_code(leafDir)
        result.set_hashcode(hashcode)
        return result
    
    @staticmethod
    def createDatasources(datasourceDir, genomeBuild="hg19", isMulticore=False, numCores=4, tx_mode="CANONICAL"):
        """
        Scrapes a directory and creates a list of datasource instances.
        
        Example datasource dir:
            dataSourceDir/gaf/hg19/
        
        Datasources are sorted as follows:
            Gaf
            Any genomic position indexed
            Any gene indexed
            
            Note: presumably any datasource that provides 'gene' and/or 'transcript' annotation(s) should be the first on the list. 
            
        Datasource config file:
        [general]
        version=1.0
        title=myDatasource
        type=tsv
        
        numCores is ignored if isMulticore == False


        """
        # TODO: Note that createDatasources does not honor the tx-mode
        dsQueueList = []

        # Get a list of all of the directories
        dsDirs = []
        if os.path.exists(datasourceDir):
            dirs = os.listdir(datasourceDir)
        else:
            logging.getLogger(__name__).warn("%s does not exist, so there will be no datasources.")
            dirs = []
        for d in dirs:
            tmpD = os.path.join(datasourceDir, d)
            fullD = os.path.join(*[tmpD, genomeBuild, ""])
            if not os.path.exists(fullD):
                logging.getLogger(__name__).warn("Potential datasource directory is missing a genome build subdirectory and will be ignored: " + tmpD)
            elif os.path.isdir(fullD):
                dsDirs.append(tmpD)
                logging.getLogger(__name__).debug("Potential datasource directory: " + fullD)
            else:
                logging.getLogger(__name__).warn("Potential datasource is not a directory: " + fullD)

        # Look for the config file.  This will have the same name as the the directory 
        # Note that the genome version is not in the dsDirs variable.
        for ds in dsDirs:
            tmpDs = ds
            if tmpDs.endswith(os.sep):
                tmpDs = tmpDs[:-1]
            configFilename = os.path.join(*[tmpDs, genomeBuild, os.path.basename(tmpDs) + ".config"])
            if not os.path.exists(configFilename):
                logging.getLogger(__name__).warn("Could not find config file for datasource: " + configFilename)
            else:
                logging.getLogger(__name__).debug("Found config file for datasource: " + configFilename)
                
                # Queue the datasource for instantiation
                logging.getLogger(__name__).info("Queuing datasource creation for " + configFilename)
                dsQueueList.append((configFilename, os.path.join(*[tmpDs, genomeBuild, ""])))

        result = []        
        if not isMulticore:
            for dsTuple in dsQueueList:
                result.append(DatasourceFactory.createDatasourceGivenTuple(dsTuple))
        else:
            result = DatasourceFactory._createDatasourcesMulticore(numCores, dsQueueList)

        return DatasourceFactory.sortDatasources(result)
    
    @staticmethod
    def _createDatasourcesMulticore(numProcesses, datasourceTuples):
        """ Private method to create the datasource list using multiple cores where possible.
        Currently, multicore functionality is only supported for gene, gp, and transcript tsvs."""
        result = []
        if len(datasourceTuples) > 0:
            logging.getLogger(__name__).info("Creating pool")
            # TODO: Create a default pool like Poco in C++.  That way multiple classes can share the same pool.
            p = LoggingPool(processes=numProcesses)
            logging.getLogger(__name__).info("Pool created")

            # Split the datasources into tmpQueue, which holds the datasources that can be initialized in parallel.
            tmpQueue = []
            tmpResult = []
            for dsTuple in datasourceTuples:
                configParser = ConfigUtils.createConfigParser(dsTuple[0]) 
                if configParser.get("general", "type") in ["gene_tsv", "gp_tsv", "gpp_tsv", "transcript_tsv"]:
                    tmpQueue.append(dsTuple)
                else:
                    result.append(DatasourceFactory.createDatasourceGivenTuple(dsTuple))
            if len(tmpQueue) > 0:
                tmpResult = p.map(createDatasource, tmpQueue)
            result.extend(tmpResult)
            logging.getLogger(__name__).info("Mapping complete: " + str(len(tmpResult)) + " datasources created in multiprocess")
            p.close()
            p.join()
        else:
            logging.getLogger(__name__).info("No datasources to initialize")
        return result
    
    @staticmethod
    def sortDatasources(datasources):
        """ 1) Make sure to put the gene-indexed datasources at the end of the list (in any order otherwise).  This is so that a gene annotation can be created ahead of annotating by gene. 
            2) Put position transforming datasources at the front (though see next step)
            3) Make sure that any Transcript datasources are put up front, but still behind ref_hg."""
        newlist = sorted(datasources, key=lambda k: isinstance(k, GenericGeneDatasource))
        newlist = sorted(newlist, key=lambda k: isinstance(k, BigWigDatasource))
        newlist = sorted(newlist, key=lambda k: isinstance(k, ChangeTransformingDatasource))
        newlist = sorted(newlist, key=lambda k: not isinstance(k, PositionTransformingDatasource))
        newlist = sorted(newlist, key=lambda k: not isinstance(k, TranscriptProvider))
        newlist = sorted(newlist, key=lambda k: not isinstance(k, ReferenceDatasource))
        return newlist
