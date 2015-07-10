from oncotator.MutationData import MutationData
from oncotator.utils.PhasingUtils import PhasingUtils
import unittest


class PhasingUtilsTest(unittest.TestCase):
    def test_phasing_info_missing(self):
        """Test whether we accurately say whether the phasing info present test works"""
        m1 = MutationData()
        m2 = MutationData()
        m3 = MutationData()
        m4 = MutationData()

        m1.createAnnotation("phasing_id", "blah")
        m2.createAnnotation("phasing_id", "blah")
        m2.createAnnotation("phasing_genotype", "0|1")
        m4.createAnnotation("phasing_genotype", "0|1")

        # m1 missing gt, m2 complete, m3 missing everything, m4 missing ID
        self.assertFalse(PhasingUtils.has_phasing_information(m1))
        self.assertTrue(PhasingUtils.has_phasing_information(m2))
        self.assertFalse(PhasingUtils.has_phasing_information(m3))
        self.assertFalse(PhasingUtils.has_phasing_information(m4))

    def test_phasing_check(self):
        """
        Test the actual phasing check.
        """
        m1 = MutationData()
        m2 = MutationData()
        m3 = MutationData()
        m4 = MutationData()
        m5 = MutationData()
        m6 = MutationData()
        m7 = MutationData()

        m1.createAnnotation("phasing_id", "blah")
        m2.createAnnotation("phasing_id", "blah")
        m2.createAnnotation("phasing_genotype", "0|1")
        m4.createAnnotation("phasing_genotype", "0|1")
        m5.createAnnotation("phasing_id", "blah")
        m5.createAnnotation("phasing_genotype", "0|1")
        m6.createAnnotation("phasing_id", "blahdifferent")
        m6.createAnnotation("phasing_genotype", "0|1")

        # m1 and m2 should not be in phase, even though they share IDs, since m1 is missing the genotype info
        unknown_val = True
        self.assertFalse(PhasingUtils.is_in_phase(m1, m2, unknown_val))
        self.assertFalse(PhasingUtils.is_in_phase(m2, m1, unknown_val))

        unknown_val = False
        self.assertFalse(PhasingUtils.is_in_phase(m1, m2, unknown_val))
        self.assertFalse(PhasingUtils.is_in_phase(m2, m1, unknown_val))

        # m2 and m4 should not be in phase, since m4 is missing the ID
        unknown_val = True
        self.assertFalse(PhasingUtils.is_in_phase(m4, m2, unknown_val))
        self.assertFalse(PhasingUtils.is_in_phase(m2, m4, unknown_val))

        unknown_val = False
        self.assertFalse(PhasingUtils.is_in_phase(m4, m2, unknown_val))
        self.assertFalse(PhasingUtils.is_in_phase(m2, m4, unknown_val))

        # m3 and m7 should be unknown_val, since phasing info is missing.
        unknown_val = True
        self.assertTrue(PhasingUtils.is_in_phase(m3, m7, unknown_val) == unknown_val)
        self.assertTrue(PhasingUtils.is_in_phase(m7, m3, unknown_val) == unknown_val)

        unknown_val = False
        self.assertTrue(PhasingUtils.is_in_phase(m3, m7, unknown_val) == unknown_val)
        self.assertTrue(PhasingUtils.is_in_phase(m7, m3, unknown_val) == unknown_val)

        # m2 and m5 should be in phase, regardless of the unknown_val parameter
        self.assertTrue(PhasingUtils.is_in_phase(m2, m5, True))
        self.assertTrue(PhasingUtils.is_in_phase(m5, m2, False))
        self.assertTrue(PhasingUtils.is_in_phase(m2, m5, False))
        self.assertTrue(PhasingUtils.is_in_phase(m5, m2, True))

        # m2 and m6 should not be in phase, since the ID is different, regardless of the unknown_val parameter
        self.assertFalse(PhasingUtils.is_in_phase(m2, m6, True))
        self.assertFalse(PhasingUtils.is_in_phase(m6, m2, False))
        self.assertFalse(PhasingUtils.is_in_phase(m2, m6, False))
        self.assertFalse(PhasingUtils.is_in_phase(m6, m2, True))


if __name__ == '__main__':
    unittest.main()
