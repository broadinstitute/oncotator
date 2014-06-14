# LICENSE_GOES_HERE


import logging
import os
from oncotator.MockExceptionThrowingDatasource import MockExceptionThrowingDatasource
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
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
from utils.MultiprocessingUtils import LoggingPool

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
                columnDataTypes[columnName] = configParser.get("data_types", columnName)

            result = IndexedTsvDatasource(src_file=filePrefix + configParser.get("general", "src_file"),
                                           title=configParser.get("general", "title"),
                                           version=configParser.get("general", "version"),
                                           colNames=columnNames,
                                           annotationColNames=annotationColumnNames,
                                           indexColNames=indexColumnNames,
                                           match_mode=configParser.get("general", "match_mode"),
                                           colDataTypes=columnDataTypes)

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
        dirs = os.listdir(datasourceDir)
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
            3) Make sure that any Gaf datasources are put up front."""
        newlist = sorted(datasources, key=lambda k: isinstance(k, GenericGeneDatasource))
        newlist = sorted(newlist, key=lambda k: isinstance(k, ChangeTransformingDatasource))
        newlist = sorted(newlist, key=lambda k: not isinstance(k, PositionTransformingDatasource))
        newlist = sorted(newlist, key=lambda k: not isinstance(k, ReferenceDatasource))
        newlist = sorted(newlist, key=lambda k: not isinstance(k, TranscriptProvider))
        return newlist