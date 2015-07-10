import unittest
from oncotator.DuplicateAnnotationException import DuplicateAnnotationException
from oncotator.MutationDataFactory import MutationDataFactory
from test.TestUtils import TestUtils

TestUtils.setupLogging(__file__, __name__)


class MutationFactoryTest(unittest.TestCase):

    _multiprocess_can_split_ = True

    def test_annotation_overwriting_on(self):
        """Test that the factory can produce a mutation that allows overwriting.  Just need to make sure no exception thrown."""
        mdf = MutationDataFactory(allow_overwriting=True)
        mut = mdf.create()

        mut.createAnnotation("blah", "123")
        self.assertTrue(mut['blah'] == "123")

        mut.createAnnotation("blah", "456")
        self.assertTrue(mut['blah'] == "456")

    def test_annotation_overwriting_off(self):
        """Test that the factory can produce a mutation that allows overwriting.  Just need to make sure no exception thrown."""
        mdf = MutationDataFactory(allow_overwriting=False)
        mut = mdf.create()

        mut.createAnnotation("blah", "123")
        self.assertTrue(mut['blah'] == "123")

        is_exception_raised = False
        try:
            mut.createAnnotation("blah", "456")
        except DuplicateAnnotationException as dae:
            is_exception_raised = True

        self.assertTrue(is_exception_raised, "DuplicateAnnotationException should have been seen, but wasn't")
