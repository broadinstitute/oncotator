# LICENSE_GOES_HERE

import os
from oncotator.utils.Hasher import Hasher
from test.TestUtils import TestUtils

__author__ = 'lichtens'

import unittest

TestUtils.setupLogging(__file__, __name__)
class HasherTest(unittest.TestCase):
    _multiprocess_can_split_ = True
    def test_directory_hash(self):
        """Test that we can read a hashcode for a directory."""
        test_dir = "testdata/thaga_janakari_gene_ds/hg19"
        hasher = Hasher()
        
        self.assertTrue(hasher.create_hashcode_for_dir(test_dir)=="cd04c3299a21ce5fbe1bb9dc0d73d907", "Hashed directory did not match ground truth. (" + hasher.create_hashcode_for_dir(test_dir) + ")  for path:  " + os.path.abspath(test_dir))

    def test_simple_hash(self):
        """Test that the single md5 call (static) functions correctly."""
        guess = Hasher.md5_hash("blah\n")
        self.assertTrue(guess == "0d599f0ec05c3bda8c3b8a68c32a1b47")

if __name__ == '__main__':
    unittest.main()
