import os
import shutil
import copy
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MutationData import MutationData
from oncotator.index.SnpOnlyLevelDbDatasourceCreator import SnpOnlyLevelDbDatasourceCreator
from oncotator.utils.GenericTsvReader import GenericTsvReader
from test.TestUtils import TestUtils

import unittest

TestUtils.setupLogging(__file__, __name__)
class SnpOnlyLevelDbDatasourceTest(unittest.TestCase):

    _multiprocess_can_split_ = True
    def _create_test_ds(self, input_tsv, dir_name, index_cols):

        base_name = "test_snp_leveldb"

        full_name = dir_name + "/" + base_name

        if os.path.exists(full_name):
            shutil.rmtree(full_name)

        os.makedirs(full_name)

        tsv_reader = GenericTsvReader(input_tsv, commentPrepend="%")
        annotation_cols = copy.copy(tsv_reader.getFieldNames())
        for icol in index_cols:
            if icol in annotation_cols:
                annotation_cols.remove(icol)

        ds_creator = SnpOnlyLevelDbDatasourceCreator()
        ds_creator.createDatasource(full_name, input_tsv, index_cols, dir_name + "/" +  base_name, "snp_leveldb", base_name, "TEST",
                         "exact", annotation_cols, [])

        config_filename = "out/test_simple_annotate_snp_only_leveldb/test_snp_leveldb/test_snp_leveldb.config"
        ds = DatasourceFactory.createDatasource(os.path.abspath(config_filename), os.path.dirname(config_filename))

        return ds

    def test_simple_annotate(self):

        ds = self._create_test_ds("testdata/small_tsv_leveldb/dbNSFP2.4_variant.chr1_cut5000.tsv", os.path.abspath("out/test_simple_annotate_snp_only_leveldb"), ["chr", "pos(1-coor)", "pos(1-coor)", "ref", "alt"])
        m = MutationData()
        # 1	35138	T	A
        m.chr = "1"
        m.start = "35138"
        m.end = "35138"
        m.ref_allele = "T"
        m.alt_allele = "A"
        m = ds.annotate_mutation(m)

        self.assertTrue(m['phyloP100way_vertebrate_rankscore'] == "0.19875")

if __name__ == '__main__':
    unittest.main()
