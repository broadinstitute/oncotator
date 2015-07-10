import logging
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MutationDataFactory import MutationDataFactory
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
                        datasource_list, genomeBuild, isMulticore, numCores,
                        defaultAnnotations, cacheUrl, read_only_cache,
                        tx_mode, is_skip_no_alts, other_opts, annotating_type):
        """


        :param inputFormat:
        :param outputFormat:
        :param inputFilename:
        :param outputFilename:
        :param globalAnnotations:
        :param datasource_list:
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

        collapse_filter_cols = other_opts.get(OptionConstants.COLLAPSE_FILTER_COLS, False)
        if  all([inputFormat == 'VCF', outputFormat == "VCF", collapse_filter_cols]):
            result.append(RunSpecificationMessage(logging.ERROR, "collapse-filter-cols flag can only be used with VCF input and non-VCF (preferably TCGAMAF/SIMPLE_TSV) output formats."))

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

        if other_opts.get(OptionConstants.REANNOTATE_TCGA_MAF_COLS) and all([inputFormat != "TCGAMAF", inputFormat !="MAFLITE"]):
            result.append(RunSpecificationMessage(logging.WARN, "Attempting to prune TCGA MAF columns from an input not specified as TCGA MAF (or MAFLITE).  This is will be ignored."))

        if (not other_opts.get(OptionConstants.REANNOTATE_TCGA_MAF_COLS) or not OptionConstants.ALLOW_ANNOTATION_OVERWRITING) and inputFormat == "TCGAMAF":
            result.append(RunSpecificationMessage(logging.ERROR, "Internal inconsistency found:  TCGAMAF input, but pruning of columns set to False or allow annotation overwriting set to false.  Please contact oncotator@broadinstitute.org"))

        if other_opts.get(OptionConstants.REANNOTATE_TCGA_MAF_COLS) and all([inputFormat != "TCGAMAF", inputFormat !="MAFLITE"]):
            result.append(RunSpecificationMessage(logging.WARN, "Asking to reannotate a tCGA MAF on an input that is not specified as a TCGA MAF.  Proceeding, but errors may result."))

        if other_opts.get(OptionConstants.REANNOTATE_TCGA_MAF_COLS) and all([outputFormat != "TCGAMAF", outputFormat !="SIMPLE_TSV"]):
            result.append(RunSpecificationMessage(logging.WARN, "Asking to reannotate a TCGA MAF for an output that is not specified as a TCGA MAF.  This is currently not supported.  Proceeding, but errors are quite likely."))

        if other_opts.get(OptionConstants.ALLOW_ANNOTATION_OVERWRITING) and all([outputFormat != "TCGAMAF", outputFormat !="SIMPLE_TSV"]):
            result.append(RunSpecificationMessage(logging.WARN, "Asking to overwrite annotations for output format that is not TCGAMAF nor SIMPLE_TSV.  This is currently not supported.  Proceeding, but errors (or inconsistent annotations) are quite likely."))

        return result

    @staticmethod
    def create_run_spec(input_format, output_format, input_filename, output_filename, global_annotations=None,
                        datasource_dir=None, genomeBuild="hg19", is_multicore=False, num_cores=4,
                        default_annotations=None, cache_url=None, read_only_cache=True,
                        tx_mode=TranscriptProvider.TX_MODE_CANONICAL, is_skip_no_alts=False, other_opts=None, annotating_type=None):
        """ This is a very simple interface to start an Oncotator session.  As a warning, this interface may notbe supported in future versions.

        If datasourceDir is None, then no datasources are used

        """
        if datasource_dir:
            datasource_list = DatasourceFactory.createDatasources(datasource_dir, genomeBuild, isMulticore=is_multicore, numCores=num_cores, tx_mode=tx_mode)
        else:
            datasource_list = []

        global_annotations = dict() if global_annotations is None else global_annotations
        default_annotations = dict() if default_annotations is None else default_annotations
        other_opts = dict() if other_opts is None else other_opts

        #TODO: Refactoring needed here to specify tx-mode (or any option not in a config file) in a cleaner way.
        for ds in datasource_list:
            if isinstance(ds, TranscriptProvider):
                logging.getLogger(__name__).info("Setting %s %s to tx-mode of %s..." % (ds.title, ds.version, tx_mode))
                ds.set_tx_mode(tx_mode)

                if other_opts.get(OptionConstants.CUSTOM_CANONICAL_TX_LIST_FILE, None) is not None:
                    cc_txs_filename = other_opts[OptionConstants.CUSTOM_CANONICAL_TX_LIST_FILE]
                    cc_txs_fp = file(cc_txs_filename, 'r')
                    cc_txs = [tx.rsplit(".", 1)[0] for tx in cc_txs_fp]
                    cc_txs_fp.close()
                    ds.set_custom_canonical_txs(cc_txs)
                    logging.getLogger(__name__).info(str(len(cc_txs)) + " custom canonical transcripts specified.")
                else:
                    logging.getLogger(__name__).info("No custom canonical transcripts specified.")

        return RunSpecificationFactory.create_run_spec_given_datasources(input_format, output_format, input_filename, output_filename, global_annotations,
                        datasource_list, genomeBuild, is_multicore, num_cores,
                        default_annotations, cache_url, read_only_cache,
                        tx_mode, is_skip_no_alts, other_opts, annotating_type)

    @staticmethod
    def create_run_spec_given_datasources(input_format, output_format, input_filename, output_filename, global_annotations=None,
                        datasource_list=None, genomeBuild="hg19", is_multicore=False, num_cores=4,
                        default_annotations=None, cache_url=None, read_only_cache=True,
                        tx_mode=TranscriptProvider.TX_MODE_CANONICAL, is_skip_no_alts=False, other_opts=None, annotating_type=None):
        """Same as create_run_spec, but a list of datasource instances can be used.  Typically, this method is only called
        by automated tests."""


        global_annotations = dict() if global_annotations is None else global_annotations
        default_annotations = dict() if default_annotations is None else default_annotations
        datasource_list = [] if datasource_list is None else datasource_list

        other_opts = dict() if other_opts is None else other_opts

        if input_format == "TCGAMAF" and not other_opts[OptionConstants.REANNOTATE_TCGA_MAF_COLS]:
            other_opts[OptionConstants.REANNOTATE_TCGA_MAF_COLS] = True

        other_opts[InputMutationCreatorOptions.IS_SKIP_ALTS] = is_skip_no_alts

        # Step 0 Validate given parameters and log messages.  If an error or critical is found, throw an exception.
        validation_messages = RunSpecificationFactory._validate_run_spec_parameters(input_format, output_format, input_filename, output_filename, global_annotations,
                        datasource_list, genomeBuild, is_multicore, num_cores,
                        default_annotations, cache_url, read_only_cache,
                        tx_mode, is_skip_no_alts, other_opts, annotating_type)
        for msg in validation_messages:
            logging.getLogger(__name__).log(msg.level, msg.message)
            if (msg.level == logging.ERROR) or (msg.level == logging.CRITICAL):
                raise RunSpecificationException(msg.message)

        # Step 1 Initialize input and output
        is_allow_annotation_overwriting = other_opts.get(OptionConstants.ALLOW_ANNOTATION_OVERWRITING, False)
        mutation_data_factory = MutationDataFactory(is_allow_annotation_overwriting)

        inputCreator = OncotatorCLIUtils.create_input_creator(input_filename, input_format, mutation_data_factory, genomeBuild, other_opts)
        outputRenderer = OncotatorCLIUtils.create_output_renderer(output_filename, output_format, other_opts)

        result = RunSpecification()
        result.initialize(inputCreator, outputRenderer, manualAnnotations=global_annotations, datasources=datasource_list,
                          isMulticore=is_multicore, numCores=num_cores, defaultAnnotations=default_annotations,
                          cacheUrl=cache_url, read_only_cache=read_only_cache, is_skip_no_alts=is_skip_no_alts, annotating_type=annotating_type,
                          is_allow_annotation_overwriting=is_allow_annotation_overwriting)
        return result