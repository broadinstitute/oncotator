from oncotator.utils.txfilter.TranscriptFilterFactory import TranscriptFilterFactory
from test.TestUtils import TestUtils

import unittest


class BasicTagTranscriptFilterTest(unittest.TestCase):
    def test_basic_tag_filtering(self):
        """Test several cases for the BasicTagTranscriptFilter"""
        tx_filter = TranscriptFilterFactory.create_instance("basic")

        ensembl_ds = TestUtils._create_test_gencode_ds("out/basic_tag_filter_ensembl_ds")
        tx_dict = ensembl_ds.getTranscriptDict()

        tx = tx_dict["ENST00000215832.6"]
        self.assertTrue(len(tx_filter.filter([tx])) == 1)

        attrib_dict = tx.get_other_attributes()
        attrib_dict.pop('tag', None)

        self.assertTrue(len(tx_filter.filter([tx])) == 0)

if __name__ == '__main__':
    unittest.main()
