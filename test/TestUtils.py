# LICENSE_GOES_HERE
from functools import wraps
import shutil
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.utils.MultiprocessingUtils import MyManager
from ConfigParser import SafeConfigParser
import os
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.datasources.Gaf import Gaf
from oncotator.datasources.dbSNP import dbSNP
from oncotator.datasources.ReferenceDatasource import ReferenceDatasource
import logging
import traceback
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory


def data_provider_decorator(fn_data_provider):
    """Data provider decorator, allows another callable to provide the data for the test.
    Modified from https://pypi.python.org/pypi/unittest-data-provider/1.0.0
    to work with nose and to accumulate assertion errors."""

    def test_decorator(fn):
        @wraps(fn)
        def repl(self, *args):
            assertion_errors = []
            ctr = 0
            for i in fn_data_provider():
                try:
                    ctr += 1
                    fn(self, *i)
                except AssertionError as ae:
                    stack_trace = traceback.format_exc()
                    assertion_errors.append("\n\n ==== Assertion error on data %s: %s -- %s\n\n%s" % (str(ctr), str(i), ae.message, stack_trace))
            if len(assertion_errors) > 0:
                raise AssertionError("\n"+"\n".join(assertion_errors) + "\n" + str(len(assertion_errors)) + " of " + str(ctr) + " tests failed.")
        return repl
    return test_decorator


class TestUtils(object):
    """
    Stateless class implementing methods that are useful in the unit tests.
    """


    def __init__(self, params):
        """
        Nothing to do when initializing this class.  
        """
        pass
    
    @staticmethod
    def createUnitTestConfig():
        config = SafeConfigParser()
        config.readfp(file('configs/default-test.config', 'r'))
        config.read(['configs/personal-test.config'])
        return config
    
    @staticmethod
    def createTranscriptProviderDatasource(config, tx_mode="CANONICAL", protocol="file"):
        """ Creates a GENCODE or Gaf 3.0 datasource from a config file.  Determines which is available automatically,
            For GAF 3.0, assumes a gaf3.0 section with keys: gaf_fname and gaf_transcript_seqs_fname

            """
        if os.path.exists(config.get("gencode", "gencodeDir")):
            gencode_dir = config.get("gencode", "gencodeDir")
            result_ds = EnsemblTranscriptDatasource(gencode_dir + "/gencode.v19.annotation.gtf", title="GENCODE", version="TEST v19", tx_filter="basic")
        else:
            gaf_fname = config.get("gaf3.0", "gaf_fname")
            gaf_transcripts_fname = config.get("gaf3.0", "gaf_transcript_seqs_fname")
            result_ds = Gaf(gaf_fname, gaf_transcripts_fname, tx_mode=tx_mode, protocol=protocol)
        return result_ds

    @staticmethod
    def createReferenceDatasource(config):
        refFilename = config.get("ref_hg", "refDir")
        return ReferenceDatasource(refFilename)

    @staticmethod
    def createGafDatasourceProxy(config, tx_mode="CANONICAL", protocol="file"):
        """ Creates a Gaf 3.0 datasource from a config file and make a proxy.
            Assumes a gaf3.0 section with keys: gaf_fname and gaf_transcript_seqs_fname
            """
        MyManager.register('Gaf', Gaf)
        manager = MyManager()
        manager.start()
        gaf_fname = config.get("gaf3.0", "gaf_fname")
        gaf_transcripts_fname = config.get("gaf3.0", "gaf_transcript_seqs_fname")
        gafDatasource = manager.Gaf(gaf_fname, gaf_transcripts_fname, tx_mode=tx_mode, protocol="file")
        return gafDatasource

    @staticmethod
    def createDbSnpDatasource(config):
        dbsnpFilename = config.get("dbSNP", "dbSNPFilename")
        return dbSNP(dbsnpFilename)

    @staticmethod
    def createCosmicDatasource(config):
        """ Creates a Cosmic datasource from a config file.
            """
        cosmic_dirname = config.get("COSMIC", "CosmicDir")
        cosmicDatasource = DatasourceFactory.createDatasource(cosmic_dirname + "/cosmic.config", cosmic_dirname)
        return cosmicDatasource

    @staticmethod
    def setupLogging(filename, package_name):
        """Static method to set up logging for a unit test.  Argument is usually __file__ from within the unit
        test file.  Should be called globally in a unit test.
        Usually called as follows::
                ...
                        TestUtils.setupLogging(__file__, __name__)
                        class MyClassTest(unittest.TestCase):
                ...

        """
        #TODO: Low priority: Use decorator instead of TestUtils.setupLogging(...)
        curdir = os.path.dirname(filename) + '/'
        logging.basicConfig(filemode='w', filename=(os.path.join(curdir, 'out/oncotator_unitTest_' + package_name + '.log')),
                            level=logging.DEBUG, format='%(asctime)s %(levelname)s [%(name)s:%(lineno)d]  %(message)s')

    @staticmethod
    def _create_test_gencode_ds(base_output_filename, protein_id_mapping_file="testdata/gencode/ensembl_id_mappingsGRCh37.p13.txt"):
        genes = ["MAPK1", "MUC16", "PIK3CA", "YPEL1", "KRTAP4-7", "MAT2A"]
        gtf_list = []
        fasta_list = []
        for gene in genes:
            gtf_list.append("testdata/gencode/" + gene + ".gencode.v18.annotation.gtf")
            fasta_list.append("testdata/gencode/" + gene + ".gencode.v18.pc_transcripts.fa")
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(gtf_list, fasta_list, base_output_filename, protein_id_mapping_file=protein_id_mapping_file)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, title="GENCODE", version="v18", tx_filter="basic")
        return ensembl_ds