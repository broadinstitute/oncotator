# LICENSE_GOES_HERE

import unittest
import logging

import pysam

from oncotator.MutationData import MutationData
from oncotator.datasources.Cosmic import Cosmic
from TestUtils import TestUtils


TestUtils.setupLogging(__file__, __name__)
class CosmicDatasourceTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.config = TestUtils.createUnitTestConfig()

    def testBasicFetching(self):
        '''Test fetching records with only a gene protein position from an existing tabix file.'''
        genePPFilename = "testdata/small_cosmic_with_gp_and_gpp/small_cosmic_trimmed_for_sorting.txt.tbi.byAA.sorted.tsv.gz"
        tabixFile = pysam.Tabixfile(genePPFilename)
        headers=["Gene_name","HGNC_ID","Sample_name","Primary_site","Site_subtype","Primary_histology","Histology_subtype","Mutation_ID","Mutation_CDS","Mutation_AA","Mutation_Description","Mutation_zygosity","Mutation_NCBI36_genome_position","Mutation_GRCh37_genome_position","Pubmed_PMID","startAA","endAA"]

        results = tabixFile.fetch(reference="EGFR", start=748, end=749)
        resultDicts = []
        for result in results:
            resultDicts.append(dict(zip(headers, result.strip().split('\t'))))
        self.assertTrue(len(resultDicts) == 2, "Should have only had two entries, but found: " + str(resultDicts))


        #A1BG	5	ME024T	NS	NS	malignant_melanoma	NS	226401	c.1132G>A	p.D378N	Substitution - Missense	unk		19:58861796-58861796	22622578
        results = tabixFile.fetch(reference="A1BG", start=377, end=378)
        resultDicts = []
        for result in results:
            resultDicts.append(dict(zip(headers, result.strip().split('\t'))))
        self.assertTrue(len(resultDicts) == 1, "Should have only had one entry, but found: " + str(resultDicts))

    def testBasicAnnotate(self):
        '''Test that the COSMIC datasource can be initialized with two index files (gp and gpp) and a simple annotation performed'''
        tabixDir = "testdata/small_cosmic_with_gp_and_gpp/"
        cosmicDS = Cosmic(src_file=tabixDir + "small_cosmic_trimmed_for_sorting.txt.tbi.gz", title="Cosmic", version="test", gpp_tabix_file= tabixDir + "small_cosmic_trimmed_for_sorting.txt.tbi.byAA.sorted.tsv.gz")

        # These values are not taken from a real world scenario, but are cooked for this test.
        m = MutationData()
        m.createAnnotation("gene", "EGFR")
        m.createAnnotation("transcript_protein_position_start", "747")
        m.createAnnotation("transcript_protein_position_end", "747")
        m.chr = '7'
        m.start = '55259560'
        m.end = '55259560'
        m = cosmicDS.annotate_mutation(m)

        self.assertTrue(m['COSMIC_n_overlapping_mutations'] == '2')

    def testMixedAnnotation(self):
        """Test that the COSMIC datasource can retrieve entries by both gp and gpp."""
        tabixDir = "testdata/small_cosmic_with_gp_and_gpp/"
        cosmicDS = Cosmic(src_file=tabixDir + "small_cosmic_trimmed_for_sorting.txt.tbi.gz", title="Cosmic", version="test", gpp_tabix_file= tabixDir + "small_cosmic_trimmed_for_sorting.txt.tbi.byAA.sorted.tsv.gz")

        # These values are not taken from a real world scenario, but are cooked for this test.
        # Line 9 should get picked up genomic coords
        # Lines 7,8 should get picked up by the protein position
        m = MutationData()
        m.createAnnotation("gene", "A2M")
        m.createAnnotation("transcript_protein_position_start", "1300")
        m.createAnnotation("transcript_protein_position_end", "1400")
        m.chr = '12'
        m.start = '9227220'
        m.end = '9227230'
        m = cosmicDS.annotate_mutation(m)

        self.assertTrue(m['COSMIC_n_overlapping_mutations'] == '3')
        self.assertTrue(m['COSMIC_overlapping_mutation_AAs'].find('1229') != -1, "Could not find the entry specified by genomic coords.")
        self.assertTrue(m['COSMIC_overlapping_primary_sites'] == "lung(3)", "Did not have the correct primary sites annotation (lung(3)): " + m['COSMIC_overlapping_primary_sites'])

    def testRealWorld(self):
        """Test that the full COSMIC datasource can retrieve entries by both gp and gpp."""
        gafDS = TestUtils.createTranscriptProviderDatasource(self.config)
        cosmicDS = TestUtils.createCosmicDatasource(self.config)

        # These values are not taken from a real world scenario, but are cooked for this test.

        m = MutationData()
        m.chr = '1'
        m.start = '12941796'
        m.end = '12941796'
        m.ref_allele = "G"
        m.alt_allele = "T"
        m = gafDS.annotate_mutation(m)
        m = cosmicDS.annotate_mutation(m)

        self.assertTrue(m['COSMIC_n_overlapping_mutations'] == '0')

        #1	150483621	150483621
        m = MutationData()
        m.chr = '1'
        m.start = '150483621'
        m.end = '150483621'
        m.ref_allele = "G"
        m.alt_allele = "T"
        m = gafDS.annotate_mutation(m)
        m = cosmicDS.annotate_mutation(m)

if __name__ == '__main__':
    unittest.main()
