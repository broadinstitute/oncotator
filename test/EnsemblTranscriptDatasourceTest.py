# LICENSE_GOES_HERE
import ConfigParser
import logging
import shutil

import unittest
from oncotator.datasources.EnsemblTranscriptDatasource import EnsemblTranscriptDatasource
from oncotator.DatasourceFactory import DatasourceFactory
from oncotator.MutationData import MutationData
from oncotator.datasources.TranscriptProvider import TranscriptProvider
from oncotator.utils.ConfigUtils import ConfigUtils
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from test.TestUtils import TestUtils

TestUtils.setupLogging(__file__, __name__)
class EnsemblTranscriptDatasourceTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_intitialize(self):
        """Test a simple initialization of an ensembl datasource """
        base_config_location = "testdata/ensembl/saccer/"
        config_parser = ConfigUtils.createConfigParser(base_config_location + "ensembl.config")
        title = config_parser.get("general", "title")
        version = config_parser.get("general", "version")
        src_file = config_parser.get("general", "src_file")

        ensembl_ds = EnsemblTranscriptDatasource(title=title, version=version, src_file=src_file)
        self.assertIsNotNone(ensembl_ds)
        ensembl_ds.set_tx_mode(TranscriptProvider.TX_MODE_BEST_EFFECT)
        self.assertTrue(TranscriptProvider.TX_MODE_BEST_EFFECT == ensembl_ds.get_tx_mode())

    def test_overlapping_single_transcripts(self):
        base_config_location = "testdata/ensembl/saccer/"

        ensembl_ds = DatasourceFactory.createDatasource(base_config_location + "ensembl.config", base_config_location)
        recs = ensembl_ds.get_overlapping_transcripts("I", "500", "500")
        self.assertTrue(len(recs) == 1)
        self.assertTrue(recs[0].get_gene() == 'YAL069W')

    def test_overlapping_multiple_transcripts_snp(self):
        base_config_location = "testdata/ensembl/saccer/"

        ensembl_ds = DatasourceFactory.createDatasource(base_config_location + "ensembl.config", base_config_location)
        recs = ensembl_ds.get_overlapping_transcripts("I", "550", "550")
        self.assertTrue(len(recs) == 2)
        ids = set()
        for r in recs:
            ids.add(r.get_transcript_id())

        self.assertTrue(len(ids - set(['YAL069W', 'YAL068W-A'])) == 0)

    def test_overlapping_multiple_transcripts_indel(self):
        base_config_location = "testdata/ensembl/saccer/"

        ensembl_ds = DatasourceFactory.createDatasource(base_config_location + "ensembl.config", base_config_location)
        recs = ensembl_ds.get_overlapping_transcripts("I", "2500", "8000")
        self.assertTrue(len(recs) == 2)
        ids = set()
        for r in recs:
            ids.add(r.get_transcript_id())

        self.assertTrue(len(ids - set(['YAL067W-A', 'YAL067C'])) == 0)

    def _create_ensembl_ds_from_saccer(self):
        gencode_input_gtf = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf"
        gencode_input_fasta = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa"
        base_output_filename = "out/test_saccer_ds"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([gencode_input_gtf], [gencode_input_fasta], base_output_filename)
        ensembl_ds = EnsemblTranscriptDatasource(base_output_filename, title="ensembl", version="71")
        return ensembl_ds

    def test_simple_annotate_with_nonhuman(self):
        """Test a very simple annotation with a nonhuman genome (saccer)"""
        ensembl_ds = self._create_ensembl_ds_from_saccer()

        m = MutationData()
        m.chr = "I"
        m.start = "500"
        m.end = "500"
        m.ref_allele = "C"
        m.alt_allele = "A"

        m2 = ensembl_ds.annotate_mutation(m)

        self.assertTrue(m2['annotation_transcript'] == "YAL069W")
        self.assertTrue(m2['gene'] == "YAL069W")

    def test_simple_annotate(self):
        """ Annotate a simple example.
        """
        base_config_location = "testdata/ensembl/saccer/"
        config_parser = ConfigUtils.createConfigParser(base_config_location + "ensembl.config")
        title = config_parser.get("general", "title")
        version = config_parser.get("general", "version")
        src_file = config_parser.get("general", "src_file")

        ensembl_ds = EnsemblTranscriptDatasource(title=title, version=version, src_file=src_file)

        m = MutationData()
        m.chr = "22"
        m.start = "22161963"
        m.end = "22161963"
        m.ref_allele = "C"
        m.alt_allele = "A"

        m2 = ensembl_ds.annotate_mutation(m)

    def test_overlapping_multiple_genes(self):
        """Test that we can collect multiple overlapping genes """
        ds = TestUtils._create_test_gencode_ds("out/overlapping_genes_multiple_")
        genes = ds.get_overlapping_genes("22", 22080000, 22120000)
        self.assertTrue(len(set(["MAPK1", "YPEL1"]) - genes) ==0 )

    def test_overlapping_gene(self):
        """Test that we can collect an overlapping gene """
        ds = TestUtils._create_test_gencode_ds("out/overlapping_genes_")
        genes = ds.get_overlapping_genes("22", 22115000, 22120000)
        self.assertTrue(len(set(["MAPK1"]) - genes) == 0)

    def test_overlapping_gene_5flank(self):
        """Test that we can collect an overlapping gene on its 5' Flank """
        ds = TestUtils._create_test_gencode_ds("out/overlapping_genes_flank")
        txs = ds.get_overlapping_transcripts("22", 22222050, 22222050, padding=100)
        self.assertTrue( len(txs) == 1)
        self.assertTrue(txs[0].get_transcript_id() == "ENST00000398822.3")

        txs = ds.get_overlapping_transcripts("22", 22224920, 22224920)
        self.assertTrue(len(txs) == 0)


    def test_small_positive_strand_transcript_change(self):
        """Test one location on a transcript and make sure that the transcript change rendered properly """
        ds = TestUtils._create_test_gencode_ds("out/small_positive_strand_")

        # Now for a negative strand
        m = MutationData()
        m.chr = "22"
        m.start = "22221730"
        m.end = "22221730"
        m.ref_allele = "T"
        m.alt_allele = "G"
        m2 = ds.annotate_mutation(m)
        self.assertTrue(m2['transcript_change'] == "c.1A>C", "Incorrect transcript change: " + m2['transcript_change'])

        # positive strand
        m = MutationData()
        m.chr = "3"
        m.start = "178916614"
        m.end = "178916614"
        m.ref_allele = "G"
        m.alt_allele = "T"
        m2 = ds.annotate_mutation(m)
        self.assertTrue(m2['transcript_change'] == "c.1G>T", "Incorrect transcript change: " + m2['transcript_change'])

    def test_hgvs_annotations_simple_SNP(self):
        """Test that HGVS annotations appear (incl. protein change) in a mutation, so we believe that the Transcript objects are populated properly."""
        ds = TestUtils._create_test_gencode_ds("out/test_hgvs_annotations_")

        # Now for a negative strand
        m = MutationData()
        m.chr = "22"
        m.start = "22221730"
        m.end = "22221730"
        m.ref_allele = "T"
        m.alt_allele = "G"
        m.build = "hg19"
        m2 = ds.annotate_mutation(m)
        self.assertEqual(m2.get('HGVS_genomic_change', None), 'chr22.hg19:g.22221730T>G')
        self.assertEqual(m2.get('HGVS_coding_DNA_change', None), 'ENST00000215832.6:c.1A>C')
        self.assertEqual(m2.get('HGVS_protein_change', None), 'ENSP00000215832:p.Met1Leu')

    def test_hgvs_annotations_IGR(self):
        """Test that the HGVS annotations appear for IGR"""
        ds = TestUtils._create_test_gencode_ds("out/test_hgvs_annotations_IGR_")
        m = MutationData()
        m.createAnnotation('variant_type', 'SNP')
        m.createAnnotation('build', 'hg19')
        m.createAnnotation('variant_classification', 'IGR')
        m.createAnnotation('chr', '15')
        m.createAnnotation('start', 30938316)
        m.createAnnotation('end', 30938316)
        m.createAnnotation('ref_allele', 'G')
        m.createAnnotation('alt_allele', 'A')
        m2 = ds.annotate_mutation(m)
        self.assertEqual(m2.get('HGVS_genomic_change', None), 'chr15.hg19:g.30938316G>A')
        self.assertEqual(m2.get('HGVS_coding_DNA_change', None), '')
        self.assertEqual(m2.get('HGVS_protein_change', None), '')

    def test_no_mapping_file(self):
        """Test that we can still create (from scratch) and instantiate a EnsemblDatasource when no protein mapping is specified (i.e. limited HGVS support)"""
        """Test that HGVS annotations appear (incl. protein change) in a mutation, so we believe that the Transcript objects are populated properly."""
        ds = TestUtils._create_test_gencode_ds("out/test_hgvs_annotations_no_mapping_", protein_id_mapping_file=None)

        # Now for a negative strand
        m = MutationData()
        m.chr = "22"
        m.start = "22221730"
        m.end = "22221730"
        m.ref_allele = "T"
        m.alt_allele = "G"
        m.build = "hg19"
        m2 = ds.annotate_mutation(m)
        self.assertEqual(m2.get('HGVS_genomic_change', None), 'chr22.hg19:g.22221730T>G')
        self.assertEqual(m2.get('HGVS_coding_DNA_change', None), 'ENST00000215832.6:c.1A>C')
        self.assertEqual(m2.get('HGVS_protein_change', None), 'unknown_prot_seq_id:p.Met1Leu')

if __name__ == '__main__':
    unittest.main()
