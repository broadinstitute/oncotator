import logging
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.input.InputMutationCreator import InputMutationCreatorOptions
from oncotator.utils.OncotatorCLIUtils import OncotatorCLIUtils
from oncotator.utils.RunSpecification import RunSpecification

class RunSpecificationFactory(object):
    """Static class for getting instances of RunSpecification"""

    def __init__(self):
        """Throw an exception since this is a static class """
        raise NotImplementedError("RunSpecificationFactory is a static class and should not be instantiated.")

    @staticmethod
    def create_run_spec(inputFormat, outputFormat, inputFilename, outputFilename, globalAnnotations=None,
                        datasourceDir=None, genomeBuild="hg19", isMulticore=False, numCores=4,
                        defaultAnnotations=None, cacheUrl=None, read_only_cache=True,
                        tx_mode=TranscriptProvider.TX_MODE_CANONICAL, is_skip_no_alts=False, other_opts=None, annotating_type=None):
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
                          cacheUrl=cacheUrl, read_only_cache=read_only_cache, is_skip_no_alts=is_skip_no_alts, annotating_type=annotating_type)
        return result