# LICENSE_GOES_HERE
from ConfigParser import SafeConfigParser
from argparse import RawTextHelpFormatter, ArgumentParser
import logging
import os
import shutil
import tempfile
from oncotator.index.GenericTsvDatasourceCreator import GenericTsvDatasourceCreator
from oncotator.utils.install.DatasourceInstallUtils import DatasourceInstallUtils
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from oncotator.utils.txfilter.TranscriptFilterFactory import TranscriptFilterFactory
from oncotator.utils.version import VERSION


def setup_logging():
    # Add a console logger to the root logger, which means that all loggers generated will have the console dump.
    #    Output on the console will be the same as what is in the log file.
    logging_format = '%(asctime)s %(levelname)s [%(name)s:%(lineno)d] %(message)s'
    logging.basicConfig(level=logging.INFO, format=logging_format)
    logging.getLogger(__name__).info("Version: " + VERSION)


def parseOptions():

    epilog = """    This utility can require a lot of RAM (~4GB for gencode.v18).
    Creation of a gencode datasource can require as much as two hours.

    NOTE about Filter:
    Please see the filter option.  Since this defaults to a GENCODE specific filter, which can be problematic for
        ENSEMBL-only.

    Use "dummy" for ENSEMBL-only datasources
    Use "basic" for GENCODE datasources, unless you want to annotate using every available transcript.

    Note that all transcripts are present in a datasource, so if a filter change is needed to a datasource that has
        already been generated, you can edit the config file, instead of re-creating the entire datasource.

    IF you wish to have HGVS support, you must provide the protein mapping file (--protein-map-file).

    """
    desc = "Create a gencode/ensembl based datasource."
    parser = ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter, epilog=epilog)
    parser.add_argument("gtf_files", type=str, help="Location of the gtf files.  Multiple files can be specified as a comma separated list (e.g. file1,file2) without spaces ")
    parser.add_argument("fasta_files", type=str, help="Location of the fasta file (cDNA) associated with the gtf files.  Multiple files can be specified as a comma separated list (e.g. file1,file2) without spaces")
    parser.add_argument("output_dir", type=str, help="Datasource output location.  This directory should NOT already exist.")
    parser.add_argument("genome_build", type=str, help="Genome build -- this must be specified correctly by the user.  For example, hg19.")
    parser.add_argument("--name", type=str, help="name of the datasource.  For example, ensembl.  Or GENCODE", default="ensembl")
    parser.add_argument("version", type=str, help="version.  For example, v18")
    parser.add_argument("--filter", type=str, help="Filter to use from " + str(TranscriptFilterFactory.TRANSCRIPT_FILTER_DICT.keys()) + ".  For non-GENCODE ENSEMBL, this should be set to dummy. default: basic", default="basic")
    parser.add_argument("-p", "--protein-map-file", type=str, help="Protein mapping file (a tsv with transcript ID to protein ID .... Typically, for ENSEMBL or GENCODE a file with ENST to ENSP mappings).")

    # Process arguments
    args = parser.parse_args()

    logger = logging.getLogger(__name__)
    logger.info("Args: " + str(args))

    return args

def main():
    setup_logging()
    args = parseOptions()
    gtf_files = args.gtf_files.split(",")
    fasta_files = args.fasta_files.split(",")
    output_dir = args.output_dir
    genome_build = args.genome_build
    name = args.name
    ver = args.version
    tx_filter = args.filter
    protein_map_file = args.protein_map_file

    # create temp dir
    tmpDir = tempfile.mkdtemp(prefix="onco_ensembl_ds_")
    try:
        logging.getLogger(__name__).info("Creating tmp dir (" + tmpDir + ") ....")
        ds_build_dir = tmpDir + "/" + genome_build + "/"
        os.mkdir(ds_build_dir)

        if not (args.gtf_files.lower().find("gencode") !=-1) and tx_filter == "basic":
            logging.getLogger(__name__).warn("basic filter requested for (apparently) a non-gencode set of GTFs.  If this is an ENSEMBL run (not GENCODE), please specify dummy, using --filter.")

        logging.getLogger(__name__).info("Creating config file...")
        config_filename = ds_build_dir + "/" + name + ".config"
        logging.getLogger(__name__).info("config file being written to: " + os.path.abspath(config_filename))

        config_file_creator = GenericTsvDatasourceCreator()
        idx_cols = DatasourceInstallUtils.indexCols("dummy_option", "dummy_values")
        config_file_creator._createConfigFile(configFilename=config_filename + ".tmp", baseDSFile=os.path.basename(gtf_files[0]),ds_type="ensembl", ds_version=ver, ds_name=name, indexCols=idx_cols)

        # Append the tx_filter and protein map file
        config_parser = SafeConfigParser()
        fp = file(config_filename + ".tmp", 'r')
        config_parser.readfp(fp)
        fp.close()
        config_parser.set("general", "transcript_filter", tx_filter)

        # Write updated config file
        fp = file(config_filename, 'w')
        config_parser.write(fp)
        fp.close()

        logging.getLogger(__name__).info("Starting index construction (temp location: " + ds_build_dir + ") ...")
        factory = GenomeBuildFactory()
        factory.construct_ensembl_indices(gtf_files, fasta_files, ds_build_dir + os.path.basename(gtf_files[0]), protein_id_mapping_file=protein_map_file)

        logging.getLogger(__name__).info("Creating datasource md5...")
        DatasourceInstallUtils.create_datasource_md5_file(ds_build_dir)


        logging.getLogger(__name__).info("Copying created datasource from temp directory to final location (" + output_dir + ")...")
        shutil.copytree(symlinks=True, src=tmpDir, dst=output_dir)

    except Exception as e:
        import traceback
        logging.getLogger(__name__).fatal((e.__repr__()) + " " + traceback.format_exc())
        logging.getLogger(__name__).info(""""If you are getting and error such as:  KeyError: 'ENST00000474204.1'), then you may be out of disk space in /tmp/.""")

    # Remove the tempdir
    logging.getLogger(__name__).info("Done...")
    logging.getLogger(__name__).info("Removing ..." + tmpDir + '/')
    shutil.rmtree(tmpDir)

if __name__ == '__main__':
    main()