from oncotator.MutationDataFactory import MutationDataFactory
from oncotator.utils.FieldMapCreator import FieldMapCreator

import unittest


class FieldMapCreatorTest(unittest.TestCase):
    def test_reannotating_with_prepends(self):
        """Test that we will disregard the prepend when looking for fields to write"""
        m = MutationDataFactory.default_create()
        m.createAnnotation('i_foo', "blah", "INPUT")
        m.createAnnotation('foo', "bloop", "some datasource")

        headers = ['i_foo']
        alt_dict = {'i_foo': ['i_i_foo', 'foo']}
        mapping = FieldMapCreator.create_field_map(headers, m, alt_dict, is_render_internal_fields=True,deprioritize_input_annotations=True)
        self.assertTrue(mapping['i_foo'] == 'foo')

        mapping = FieldMapCreator.create_field_map(headers, m, alt_dict, is_render_internal_fields=True,deprioritize_input_annotations=False)
        self.assertTrue(mapping['i_foo'] == 'i_foo')


if __name__ == '__main__':
    unittest.main()
