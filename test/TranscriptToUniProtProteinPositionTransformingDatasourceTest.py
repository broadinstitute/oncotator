# LICENSE_GOES_HERE


from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MutationData import MutationData
from oncotator.datasources.TranscriptToUniProtProteinPositionTransformingDatasource import TranscriptToUniProtProteinPositionTransformingDatasource
from TestUtils import TestUtils

__author__ = 'lichtens'

import unittest

TestUtils.setupLogging(__file__, __name__)
class TranscriptToUniProtProteinPositionTransformingDatasourceTest(unittest.TestCase):

    _multiprocess_can_split_ = True

    def tearDown(self):
        pass

    def testBasicAnnotationNoChange(self):
        """ Test whether we can translate from one coordinate system to another.  This tests no change.
        """
        tDS = TranscriptToUniProtProteinPositionTransformingDatasource(title="UniProt", version="test", src_file="file://testdata/small_uniprot_prot_seq_ds/db")

        # Must correspond to what the datasource is going to generate
        outputAnnotation = "UniProt_aapos"
        m = MutationData()
        m.createAnnotation('transcript_id', 'uc003tqk.3')
        m.createAnnotation('protein_change', 'p.S50T')
        m = tDS.annotate_mutation(m)
        self.assertTrue(m[outputAnnotation] == "50", "Did not get proper value (50): " + m[outputAnnotation])

    def testBasicAnnotationWithChange(self):
        """ Test whether we can translate from one coordinate system to another.  This tests a known change.
        """
        tDS = TranscriptToUniProtProteinPositionTransformingDatasource(title="UniProt", version="test", src_file="file://testdata/small_uniprot_prot_seq_ds/db")

        # Must correspond to what the datasource is going to generate.
        outputAnnotation = "UniProt_aapos"
        m = MutationData()
        m.createAnnotation('transcript_id', 'uc009vvt.1')
        m.createAnnotation('protein_change', 'p.T1105A')
        m = tDS.annotate_mutation(m)
        self.assertTrue(m[outputAnnotation] == "969", "Did not get proper value (969): " + m[outputAnnotation])

    def testDatasourceCreator(self):
        """ Test that the datasource creator process will work for  TranscriptToUniProtProteinPositionTransformingDatasource.  NOTE: This test needs to be updated to use sqlite instead of filesystem file.
        """

        tDS = DatasourceFactory.createDatasource("testdata/small_uniprot_prot_seq_ds/small_uniprot_prot_seq_ds.config", "testdata/small_uniprot_prot_seq_ds/")
        outputAnnotation = "UniProt_aapos"
        m = MutationData()
        m.createAnnotation('transcript_id', 'uc009vvt.1')
        m.createAnnotation('protein_change', 'p.T1105A')
        m = tDS.annotate_mutation(m)
        self.assertTrue(m[outputAnnotation] == "969", "Did not get proper value (969): " + m[outputAnnotation])


if __name__ == '__main__':
    unittest.main()
