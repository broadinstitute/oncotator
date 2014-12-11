import logging
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.input.InputMutationCreator import InputMutationCreatorOptions
from oncotator.utils.OncotatorCLIUtils import OncotatorCLIUtils
from oncotator.utils.OptionConstants import OptionConstants
from oncotator.utils.RunSpecification import RunSpecification
from oncotator.utils.RunSpecificationException import RunSpecificationException
from oncotator.utils.RunSpecificationMessage import RunSpecificationMessage


class RunSpecificationFactory(object):
    """Static class for getting instances of RunSpecification"""

    def __init__(self):
        """Throw an exception since this is a static class """
        raise NotImplementedError("RunSpecificationFactory is a static class and should not be instantiated.")

    @staticmethod
    def _validate_run_spec_parameters(inputFormat, outputFormat, inputFilename, outputFilename, globalAnnotations,
                        datasourceDir, genomeBuild, isMulticore, numCores,
                        defaultAnnotations, cacheUrl, read_only_cache,
                        tx_mode, is_skip_no_alts, other_opts, annotating_type):
        """


        :param inputFormat:
        :param outputFormat:
        :param inputFilename:
        :param outputFilename:
        :param globalAnnotations:
        :param datasourceDir:
        :param genomeBuild:
        :param isMulticore:
        :param numCores:
        :param defaultAnnotations:
        :param cacheUrl:
        :param read_only_cache:
        :param tx_mode:
        :param is_skip_no_alts:
        :param other_opts:
        :param annotating_type:
        :returns list: List of RunSpecificationMessage
        """
        result = []
        if outputFormat=="TCGAVCF":
            result.append(RunSpecificationMessage(logging.WARN, "TCGA VCF output is not supported and should be considered experimental when used outside of the Broad Institute.  Outside of the Broad Institute, use of -o VCF is more likely to be desired by users."))

        if is_skip_no_alts and (outputFormat == "VCF"):
            result.append(RunSpecificationMessage(logging.WARN, "--skip-no-alt specified when output is a VCF.  This is likely to generate errors."))
        if is_skip_no_alts and (inputFormat != "VCF"):
            result.append(RunSpecificationMessage(logging.INFO, "--skip-no-alt specified when input is not VCF.  skip-no-alt is not going to do anything."))

        is_no_prepend = other_opts.get(OptionConstants.NO_PREPEND, True)
        if is_no_prepend and (outputFormat != "TCGAMAF"):
            result.append(RunSpecificationMessage(logging.INFO, "no prepend specified when output is not TCGAMAF.  Ignoring and proceeding."))

        if outputFormat=="GENE_LIST" and inputFormat != "SEG_FILE":
            result.append(RunSpecificationMessage(logging.ERROR, "Output format of GENE_LIST is only supported when input is SEG_FILE"))

        if outputFormat not in ["GENE_LIST", "SIMPLE_TSV"] and inputFormat == "SEG_FILE":
            result.append(RunSpecificationMessage(logging.WARN, "Input format of SEG_FILE is only supported when output is GENE_LIST or SIMPLE_TSV"))

        if inputFormat == "VCF" and outputFormat == "VCF" and other_opts.get(OptionConstants.VCF_OUT_INFER_GENOTYPES):
            result.append(RunSpecificationMessage(logging.WARN,"Infer genotypes option has been set to true.  "
                        "Because the input is a VCF file, infer genotypes will have no effect on the output."))
        elif inputFormat != "VCF" and outputFormat == "VCF" and not other_opts.get(OptionConstants.VCF_OUT_INFER_GENOTYPES):
            result.append(RunSpecificationMessage(logging.WARN,"Infer genotypes option has been set to false.  "
                        "Because the input is not a VCF file, genotype field may not be rendered properly."))

        if outputFormat == "VCF" and inputFormat == "VCF" and other_opts.get(OptionConstants.INFER_ONPS):
            result.append(RunSpecificationMessage(logging.WARN,("Inferring ONPs may cause issues with VCF->VCF annotation.")))

        if other_opts.get(OptionConstants.INFER_ONPS):
            result.append(RunSpecificationMessage(logging.INFO,("Output file order may be different than input file because we're combining SNPs into ONPs")))

        return result

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

        globalAnnotations = dict() if globalAnnotations is None else globalAnnotations
        defaultAnnotations = dict() if defaultAnnotations is None else defaultAnnotations
        other_opts = dict() if other_opts is None else other_opts

        other_opts[InputMutationCreatorOptions.IS_SKIP_ALTS] = is_skip_no_alts

        # Step 0 Validate given parameters and log messages.  If an error or critical is found, throw an exception.
        validation_messages = RunSpecificationFactory._validate_run_spec_parameters(inputFormat, outputFormat, inputFilename, outputFilename, globalAnnotations,
                        datasourceDir, genomeBuild, isMulticore, numCores,
                        defaultAnnotations, cacheUrl, read_only_cache,
                        tx_mode, is_skip_no_alts, other_opts, annotating_type)
        for msg in validation_messages:
            logging.getLogger(__name__).log(msg.level, msg.message)
            if (msg.level == logging.ERROR) or (msg.level == logging.CRITICAL):
                raise RunSpecificationException(msg.message)

        # Step 1 Initialize input and output
        inputCreator = OncotatorCLIUtils.create_input_creator(inputFilename, inputFormat, genomeBuild, other_opts)
        outputRenderer = OncotatorCLIUtils.create_output_renderer(outputFilename, outputFormat, other_opts)

        # Step 2 Datasources
        if datasourceDir:
            datasource_list = DatasourceFactory.createDatasources(datasourceDir, genomeBuild, isMulticore=isMulticore, numCores=numCores, tx_mode=tx_mode)
        else:
            datasource_list = []

        #TODO: Refactoring needed here to specify tx-mode (or any option not in a config file) in a cleaner way.
        for ds in datasource_list:
            if isinstance(ds, TranscriptProvider):
                logging.getLogger(__name__).info("Setting %s %s to tx-mode of %s..." % (ds.title, ds.version, tx_mode))
                ds.set_tx_mode(tx_mode)

                if other_opts.get(OptionConstants.CUSTOM_CANONICAL_TX_LIST_FILE, None) is not None:
                    cc_txs_fp = file(other_opts[OptionConstants.CUSTOM_CANONICAL_TX_LIST_FILE], 'r')
                    cc_txs = [tx.rsplit(".", 1)[0] for tx in cc_txs_fp]
                    cc_txs_fp.close()
                    ds.set_custom_canonical_txs(cc_txs)
                    logging.getLogger(__name__).info(str(len(cc_txs)) + " custom canonical transcripts specified.")
                else:
                    logging.getLogger(__name__).info("No custom canonical transcripts specified.")

        result = RunSpecification()
        result.initialize(inputCreator, outputRenderer, manualAnnotations=globalAnnotations, datasources=datasource_list,
                          isMulticore=isMulticore, numCores=numCores, defaultAnnotations=defaultAnnotations,
                          cacheUrl=cacheUrl, read_only_cache=read_only_cache, is_skip_no_alts=is_skip_no_alts, annotating_type=annotating_type)
        return result