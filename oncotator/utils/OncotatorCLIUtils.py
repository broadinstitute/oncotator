# LICENSE_GOES_HERE
import os
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.input.InputMutationCreator import InputMutationCreatorOptions
from oncotator.input.VcfInputMutationCreator import VcfInputMutationCreator
from oncotator.output.VcfOutputRenderer import VcfOutputRenderer
import logging


"""
Class file for OncotatorCLIUtils


Created on Dec 19, 2012

@author: lichtens


"""

# Step 1
from oncotator.input.MafliteInputMutationCreator import MafliteInputMutationCreator
from oncotator.output.TcgaMafOutputRenderer import TcgaMafOutputRenderer
from ConfigUtils import ConfigUtils

from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.output.SimpleOutputRenderer import SimpleOutputRenderer
from oncotator.output.SimpleBedOutputRenderer import SimpleBedOutputRenderer
from oncotator.output.TcgaVcfOutputRenderer import TcgaVcfOutputRenderer


class RunSpecification(object):
    """ This class contains a specification for an oncotator run.  
    
    This class is quite simple and simply holds parameters.  Annotator classes will be expected to be able to initialize and annotate
        from an instance of this class.
        
        inputCreator -- instance
        outputCreator -- instance
        manualAnnotations -- A dict of name-value pairs that will be annotated onto all mutations.
            name is the name of the annotation.
            value is the value of the annotation.
            The source will always be listed as MANUAL
            TODO: The Annotator takes care of some of this.  Move the relevant documentation.
        datasources -- A list of datasources (instance of Datasource).
        isMulticore -- use multicore processing, where available (True/False)
        numCores -- number of cores to use if isMulticore is True.  Otherwise, this is ignored.
        cache_url -- if None, implies that there is no cache.
        """
    def __init__(self):
        self.__inputCreator = None
        self.__outputRenderer = None
        self.__inputFilename = None
        self.__outputFilename = None
        self.__manualAnnotations = None
        self.__defaultAnnotations = None
        self.__datasources = None
        self.__isMulticore = False
        self.__numCores = None
        self.__cache_url = None
        self.__is_read_only_cache=True
        self.__is_skip_no_alts = False
        pass

    def get_cache_url(self):
        return self.__cache_url

    def set_cache_url(self, value):
        self.__cache_url = value

    def del_cache_url(self):
        del self.__cache_url

    def set_is_read_only_cache(self, value):
        self.__is_read_only_cache = value

    def get_is_read_only_cache(self):
        return self.__is_read_only_cache

    def del_is_read_only_cache(self):
        del self.__is_read_only_cache

    def get_is_multicore(self):
        return self.__isMulticore


    def get_num_cores(self):
        return self.__numCores


    def set_is_multicore(self, value):
        self.__isMulticore = value


    def set_num_cores(self, value):
        self.__numCores = value


    def del_is_multicore(self):
        del self.__isMulticore


    def del_num_cores(self):
        del self.__numCores


    def get_datasources(self):
        return self.__datasources


    def set_datasources(self, value):
        self.__datasources = value


    def del_datasources(self):
        del self.__datasources


    def get_input_creator(self):
        return self.__inputCreator


    def get_output_renderer(self):
        return self.__outputRenderer


    def get_manual_annotations(self):
        return self.__manualAnnotations


    def set_input_creator(self, value):
        self.__inputCreator = value


    def set_output_renderer(self, value):
        self.__outputRenderer = value

    def set_manual_annotations(self, value):
        self.__manualAnnotations = value

    def del_input_creator(self):
        del self.__inputCreator

    def del_output_renderer(self):
        del self.__outputRenderer

    def del_manual_annotations(self):
        del self.__manualAnnotations

    def get_default_annotations(self):
        return self.__defaultAnnotations

    def set_default_annotations(self, value):
        self.__defaultAnnotations = value

    def del_default_annotations(self):
        del self.__defaultAnnotations

    def get_is_skip_no_alts(self):
        return self.__is_skip_no_alts

    def set_is_skip_no_alts(self, value):
        self.__is_skip_no_alts = value

    def del_is_skip_no_alts(self):
        del self.__is_skip_no_alts

    def initialize(self, inputCreator, outputRenderer, manualAnnotations=None, datasources=None, isMulticore=False, numCores=4, defaultAnnotations=None, cacheUrl=None, read_only_cache=True, is_skip_no_alts=False):
        self.inputCreator = inputCreator
        self.outputRenderer = outputRenderer
        self.manualAnnotations = manualAnnotations if manualAnnotations is not None else dict()
        self.datasources = datasources if datasources is not None else []
        self.isMulticore = isMulticore
        self.numCores = numCores
        self.defaultAnnotations = defaultAnnotations if defaultAnnotations is not None else dict()
        self.cacheUrl = cacheUrl
        self.isReadOnlyCache = read_only_cache
        self.isSkipNoAlts = is_skip_no_alts

    
    inputCreator = property(get_input_creator, set_input_creator, del_input_creator, "inputCreator's docstring")
    outputRenderer = property(get_output_renderer, set_output_renderer, del_output_renderer, "outputRenderer's docstring")
    manualAnnotations = property(get_manual_annotations, set_manual_annotations, del_manual_annotations, "manualAnnotations's docstring")
    defaultAnnotations = property(get_default_annotations, set_default_annotations, del_default_annotations, "Annotations that are populated only when the annotation does not exist or is empty ('' or None)")
    datasources = property(get_datasources, set_datasources, del_datasources, "datasources's docstring")
    isMulticore = property(get_is_multicore, set_is_multicore, del_is_multicore, "isMulticore's docstring")
    numCores = property(get_num_cores, set_num_cores, del_num_cores, "numCores's docstring")
    cacheUrl = property(get_cache_url, set_cache_url, del_cache_url, "cacheUrl's docstring")
    isReadOnlyCache = property(get_is_read_only_cache, set_is_read_only_cache, del_is_read_only_cache, "isReadOnlyCache's docstring")
    isSkipNoAlts = property(get_is_skip_no_alts, set_is_skip_no_alts, del_is_skip_no_alts, "isSkipNoAlts's docstring")


class OncotatorCLIUtils(object):
    """
    Utility methods for implementing a command line interface (or any presentation layer class).
    
    Provides utilities for starting an Oncotator session given a set of parameters.
    
    IMPORTANT: If more input/output formats are needed, simply add to the dicts created in  
        createInputFormatNameToClassDict()
        createOutputFormatNameToClassDict()

        First parameter is the class name and the second is the default config file.  If no config file is needed, then
            just specify empty string.  Constructors of the classes must be able to take a config file, even if it
            is ignored.
        
        That will enable support for a new type.  The CLI drives the available 
            formats from the above methods. key is a unique string (e.g. TCGAMAF)
            and the value is the class name (NOT in quotes).  See the method code.
        
    """
    
    def __init__(self, params):
        """
        Constructor -- Never use this.  All methods should be called from a static context.  Throws an exception.
        """
        raise NotImplementedError('This class should not be instantiated.  All methods are static.')

    @staticmethod
    def createInputFormatNameToClassDict():
        """ Poor man's dependency injection. Change this method to support 
        more input formats."""
        return {'MAFLITE': (MafliteInputMutationCreator, 'maflite_input.config'),
                "VCF": (VcfInputMutationCreator, 'vcf.in.config')}

    @staticmethod
    def createOutputFormatNameToClassDict():
        """ Poor man's dependency injection. Change this method to support 
        more output formats."""
        return {'TCGAMAF': (TcgaMafOutputRenderer, "tcgaMAF2.4_output.config"),
                "SIMPLE_TSV": (SimpleOutputRenderer, ""),
                'SIMPLE_BED': (SimpleBedOutputRenderer, ""),
                'TCGAVCF': (TcgaVcfOutputRenderer, "tcgaVCF1.1_output.config"),
                'VCF': (VcfOutputRenderer, "vcf.out.config")}

    @staticmethod
    def getSupportedOutputFormats():
        """ Lists the supported output formats """
        tmp = OncotatorCLIUtils.createOutputFormatNameToClassDict()
        return tmp.keys()     
    
    @staticmethod
    def getSupportedInputFormats():
        """ Lists the supported input formats """
        tmp = OncotatorCLIUtils.createInputFormatNameToClassDict()
        return tmp.keys()

    @staticmethod
    def create_input_creator(inputFilename, inputFormat, genome_build="hg19", input_creator_options=None):
        inputCreatorDict = OncotatorCLIUtils.createInputFormatNameToClassDict()
        if inputFormat not in inputCreatorDict.keys():
            raise NotImplementedError("The inputFormat specified: " + inputFormat + " is not supported.")
        else:
            inputConfig = inputCreatorDict[inputFormat][1]
            inputCreator = inputCreatorDict[inputFormat][0](inputFilename, inputConfig, genome_build, input_creator_options)
        return inputCreator

    @staticmethod
    def create_output_renderer(outputFilename, outputFormat, otherOptions):
        outputRendererDict = OncotatorCLIUtils.createOutputFormatNameToClassDict()
        if outputFormat not in outputRendererDict.keys():
            raise NotImplementedError("The outputFormat specified: " + outputFormat + " is not supported.")
        else:
            outputConfig = outputRendererDict[outputFormat][1]
            outputRenderer = outputRendererDict[outputFormat][0](outputFilename, outputConfig, otherOptions)
        return outputRenderer

    @staticmethod
    def create_run_spec(inputFormat, outputFormat, inputFilename, outputFilename, globalAnnotations=None,
                        datasourceDir=None, genomeBuild="hg19", isMulticore=False, numCores=4,
                        defaultAnnotations=None, cacheUrl=None, read_only_cache=True,
                        tx_mode=TranscriptProvider.TX_MODE_CANONICAL, is_skip_no_alts=False, other_opts=None):
        """ This is a very simple interface to start an Oncotator session.  As a warning, this interface may notbe supported in future versions.
        
        If datasourceDir is None, then the default location is used.  TODO: Define default location.
        
        IMPORTANT: Current implementation attempts to annotate using a default set of datasources.
        
        TODO: Make sure that this note above is no longer the case.  Current implementation attempts to annotate using a default set of datasources
        TODO: This method may get refactored into a separate class that handles RunConfigutaion objects. 
        """  
        # TODO: Use dependency injection for list of name value pairs?  Otherwise, set it up as an attribute on this class.
        # TODO: Use dependency injection to return instance of the input/output classes
        # TODO: Support more than the default configs.
        # TODO: On error, list the supported formats (both input and output) 
        # TODO: Make sure that we can pass in both a class and a config file, not just a class.

        globalAnnotations = dict() if globalAnnotations is None else globalAnnotations
        defaultAnnotations = dict() if defaultAnnotations is None else defaultAnnotations
        other_opts = dict() if other_opts is None else other_opts

        other_opts[InputMutationCreatorOptions.IS_SKIP_ALTS] = is_skip_no_alts

        # Step 1 Initialize input and output
        inputCreator = OncotatorCLIUtils.create_input_creator(inputFilename, inputFormat, genomeBuild, other_opts)
        outputRenderer = OncotatorCLIUtils.create_output_renderer(outputFilename, outputFormat, other_opts)

        # Step 2 Datasources
        datasourceList = DatasourceFactory.createDatasources(datasourceDir, genomeBuild, isMulticore=isMulticore, numCores=numCores, tx_mode=tx_mode)

        #TODO: Refactoring needed here to specify tx-mode (or any option not in a config file) in a cleaner way.
        for ds in datasourceList:
            if isinstance(ds, TranscriptProvider):
                logging.getLogger(__name__).info("Setting %s %s to tx-mode of %s..." % (ds.title, ds.version, tx_mode))
                ds.set_tx_mode(tx_mode)

        result = RunSpecification()
        result.initialize(inputCreator, outputRenderer, manualAnnotations=globalAnnotations, datasources=datasourceList,
                          isMulticore=isMulticore, numCores=numCores, defaultAnnotations=defaultAnnotations,
                          cacheUrl=cacheUrl, read_only_cache=read_only_cache, is_skip_no_alts=is_skip_no_alts)
        return result
    
    @staticmethod
    def createManualAnnotationsGivenConfigFile(configFile):
        """
        Assumes a config file as:
        [manual_annotations]
        # annotation3 has blank ('') value
        override:annotation1=value1,annotation2=value2,annotation3=,annotation4=value4
        
        Returns a dictionary:
        {annotation1:value1,annotation2:value2,annotation3:'',annotation4=value4}
        """
        if (configFile is None) or (configFile.strip() == ''):
            return dict()

        if not os.path.exists(configFile):
            logging.getLogger(__name__).warn("Could not find annotation config file: " + configFile + "  ... Proceeding without it.")
            return dict()

        config = ConfigUtils.createConfigParser(configFile)
        opts = config.get('manual_annotations','override').split(',')
        result = dict()
        for optTmp in opts:
            opt = optTmp.split('=')
            if (len(opt) == 1) or (opt[1] is None):
                opt[1] = '' 
            result[opt[0]] = opt[1]
        return result

    @staticmethod
    def determineAllAnnotationValues(commandLineManualOverrides, overrideConfigFile):
        manualOverrides = OncotatorCLIUtils.createManualAnnotationsGivenConfigFile(overrideConfigFile)
        for clmo in commandLineManualOverrides:
            if clmo.find(":") == -1:
                logging.getLogger(__name__).warn("Could not parse annotation: " + str(clmo) + "   ... skipping")
                continue
            keyval = clmo.split(':', 1)
            manualOverrides[keyval[0]] = keyval[1]

        return manualOverrides