from shove.core import Shove
from oncotator.MutationData import MutationData
from oncotator.datasources.ShoveDatasource import ShoveDatasource

__author__ = 'lichtens'

import unittest


class ShoveDatasourceTest(unittest.TestCase):
    def test_creation_scrap(self):

        # chr, start, end, ref, alt
        protocol = "sqlite"
        output_url = protocol + ":///dbNSFP_chr1.sqlite"

        shove_tmp = Shove(output_url, "simple://", max_entries=300)

        index_cols = ["#chr", "pos(1-coor)", "pos(1-coor)", "ref", "alt"]
        tsv_file = "/xchip/cga_home/mgupta/mgupta/dbNSFP/dbNSFP2.4_variant.chr1"
        fp = file(tsv_file, 'r')
        for i,line in enumerate(fp):
            line_list = line.split("\t")
            line_list[-1] = line_list[-1].strip()
            m = MutationData()
            m.chr = line_list[0]
            m.start = line_list[1]
            m.end = line_list[1]
            m.ref_allele = line_list[2]
            m.alt_allele = line_list[3]
            h = ShoveDatasource.generate_hash(m)
            if i % 5000 == 0:
                print("Rendering %d" % (i))
            try:
                shove_tmp[h].append(line_list[4:])
            except KeyError:
                shove_tmp[h] = [line_list[4:]]

        shove_ds = ShoveDatasource(title="dbNSFP_chr1_test_" + protocol, version = "TEST", src_file=output_url, index_cols=index_cols)

if __name__ == '__main__':
    unittest.main()
