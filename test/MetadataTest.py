# LICENSE_GOES_HERE
import unittest
from oncotator.Annotation import Annotation
from oncotator.Metadata import Metadata
from TestUtils import TestUtils

TestUtils.setupLogging(__file__, __name__)
class MetadataTest(unittest.TestCase):

    def test_simple_create(self):
        """Test absurdly simple creation"""
        md = Metadata()
        md['fake1'] = Annotation("1")

        self.assertEqual(md['fake1'].value, "1")


if __name__ == '__main__':
    unittest.main()
