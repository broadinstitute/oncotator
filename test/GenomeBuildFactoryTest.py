# LICENSE_GOES_HERE
import os
import shutil
import unittest
from shove.core import Shove
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from oncotator.index.gaf import region2bin,region2bins
from oncotator.utils.MutUtils import MutUtils

class GenomeBuildFactoryTest(unittest.TestCase):

    def test_build_ensembl_transcript_index(self):
        """Build the gtf portion of the ensembl transcript db
        """
        # cat ~/oncotator_pycharm/oncotator/test/testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf | cut -f 9 | cut -f 5 --delimiter=" " | sort | uniq | sed -r "s/;//g" | sed -r "s/\"//g"
        #  snR84, tK(UUU)K, YAL067C, YAL067W-A, YAL068C, YAL068W-A, YAL069W, YBR278W, YBR279W, YBR280C, YBR281C, YDR528W, YDR529C, YKR074W,
        #
        # grep -Pzo  ">(snR84|tK\(UUU\)K|YAL067C|YAL067W-A|YAL068C|YAL068W-A|YAL069W|YBR278W|YBR279W|YBR280C|YBR281C|YDR528W|YDR529C|YKR074W)([A-Za-z_0-9 \:\-\n]+)" Saccharomyces_cerevisiae.EF4.71.cdna.all.fa >Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa
        #
        ensembl_input_gtf = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf"
        ensembl_input_fasta = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa"

        output_filename = "out/test_ensembl_gtf.db"
        protocol = "file"
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.build_ensembl_transcript_index([ensembl_input_gtf], [ensembl_input_fasta], output_filename, protocol=protocol)
        self.assertTrue(os.path.exists(output_filename))

        shove = Shove(protocol + "://" + output_filename, "memory://")
        self.assertTrue(len(shove.keys()) > 0)
        self.assertTrue("YDR529C" in shove.keys())
        t = shove["YDR529C"]
        self.assertTrue(t.get_seq() is not None)
        self.assertTrue(t.get_seq() is not "")
        self.assertTrue(len(t.get_cds()) > 0)
        self.assertTrue(len(t.get_exons()) > 0)
        MutUtils.removeDir(output_filename)

    def test_region2bin(self):
        """Simple test that the region2bin works for genomic position indexing """

        # Footprint for PIK3CA transcript chr3:178,866,311-178,952,497  uc003fjk.3
        guess = region2bin(178866311, 178952497)

        self.assertTrue(guess == 243)

    def test_build_ensembl_transcripts_by_gene_index(self):
        """Test building an index for getting a transcript given a gene."""
        protocol = "file"
        transcript_index_filename = "out/test_ensembl_gtf_for_gene.db"
        output_filename = "out/test_ensembl_gtf_for_gene.db.gene.idx"
        shutil.rmtree(output_filename,ignore_errors=True)

        ensembl_input_gtf = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf"
        ensembl_input_fasta = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa"

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.build_ensembl_transcript_index([ensembl_input_gtf], [ensembl_input_fasta], transcript_index_filename, protocol=protocol)
        genome_build_factory.build_ensembl_transcripts_by_gene_index(transcript_index_filename, output_filename)

        # Now load the index and look something up.
        gene_index = Shove(protocol + "://" + output_filename, optimize=False)
        self.assertTrue(len(gene_index['SEO1']) == 1)
        tx = gene_index['SEO1'][0]

        self.assertTrue(tx.get_transcript_id()=="YAL067C")

    def test_build_ensembl_transcripts_by_genomic_location_index(self):
        """Test that we can get an ensembl transcript from a genomic position"""
        protocol = "file"
        transcript_index_filename = "out/test_ensemble_gtf_for_gp.db"
        output_filename = "out/test_ensemble_gtf_for_gp.db.idx"
        shutil.rmtree(output_filename, ignore_errors=True)

        ensembl_input_gtf = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf"
        ensembl_input_fasta = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa"

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.build_ensembl_transcript_index([ensembl_input_gtf], [ensembl_input_fasta], transcript_index_filename, protocol=protocol)
        genome_build_factory.build_ensembl_transcripts_by_genomic_location_index(transcript_index_filename, output_filename, protocol=protocol)

        # Now load the index and look something up.
        gp_index = Shove(protocol + "://" + output_filename)
        gt_transcript_id = "YAL067C"
        bins = region2bins(1496172, 1496400)

        for bin in bins:
            key = 'I_' + str(bin)
            if key in gp_index.keys():
                self.assertTrue(gp_index[key] == gt_transcript_id)

    def test_construct_full_indices(self):
        """Attempt to construct all three ensembl indices with one command. """
        ensembl_input_gtf = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf"
        ensembl_input_fasta = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa"
        base_output_filename = "out/test_full_indices_ensembl"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([ensembl_input_gtf], [ensembl_input_fasta], base_output_filename)

        self.assertTrue(os.path.exists(base_output_filename + ".transcript.idx"))
        self.assertTrue(os.path.exists(base_output_filename + ".transcript_by_gene.idx"))
        self.assertTrue(os.path.exists(base_output_filename + ".transcript_by_gp_bin.idx"))

    def test_retrieving_sequence(self):
        """Ensure we can retrieve a sequence from an ensembl transcript given a gene.  """

        ensembl_input_gtf = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf"
        ensembl_input_fasta = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa"
        base_output_filename = "out/test_retrieving_full_indices_ensembl"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([ensembl_input_gtf], [ensembl_input_fasta], base_output_filename)

        seq_index = Shove("file://" + base_output_filename + ".transcript_by_gene.idx", optimize=False)
        transcripts = seq_index['SEO1']
        transcript = transcripts[0]
        for i in xrange(len(transcripts)):
            transcript = transcripts[i]
            if transcript._transcript_id == "YAL067C":
                break
        self.assertTrue(transcript.get_seq().startswith('ATGTATTCAATTGTTAAAGAGATTATTGTAGATCCTTACAAAAGACTAAAATGGGGTTTT'))

        transcripts = seq_index['PAU8']
        transcript = transcripts[0]
        for i in xrange(len(transcripts)):
            transcript = transcripts[i]
            if transcript._transcript_id == "YAL068C":
                break
        self.assertTrue(transcript.get_strand() == "-")

        seq_index_gp = Shove("file://" + base_output_filename + ".transcript_by_gp_bin.idx", "memory://")
        transcripts = seq_index_gp["I_585"]
        self.assertTrue(len(transcripts) == 5, "There should be 5 transcripts.")
        transcript = transcripts[0]
        for i in xrange(len(transcripts)):
            transcript = transcripts[i]
            if transcript._transcript_id == "YAL069W":
                break
        self.assertTrue(transcript.get_strand() == "+")

    def test_gencode_small(self):
        """Test that we can create Transcript instances from a small gencode gtf and fasta."""
        gencode_input_gtf = "testdata/gencode/MAPK1.gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/MAPK1.gencode.v18.pc_transcripts.fa"
        base_output_filename = "out/test_small_gencode"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([gencode_input_gtf], [gencode_input_fasta], base_output_filename)

        seq_index = Shove("file://" + base_output_filename + ".transcript_by_gene.idx", "memory://", optimize=False)
        transcripts = seq_index["MAPK1"]
        self.assertTrue(len(transcripts) == 4)

        seq_index_gp = Shove("file://" + base_output_filename + ".transcript_by_gp_bin.idx", "memory://", optimize=False)
        transcripts = seq_index_gp["22_753"]
        self.assertTrue(transcripts[0].get_strand() == "-")
        self.assertTrue(len(transcripts) == 1)

        for tx in transcripts:
            if tx.get_transcript_id() != "ENST00000215832.6":
                continue
            self.assertTrue(tx.get_seq().startswith("AGGCAATCGGTCCGAG"))

    def test_gencode_cp(self):
        """Test the indexing of a gene that was causing problems and make sure that it can be indexed."""
        gencode_input_gtf = "testdata/gencode/CP.gencode.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/CP.gencode.pc_transcripts.fa"
        base_output_filename = "out/test_cp_gencode"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices([gencode_input_gtf], [gencode_input_fasta], base_output_filename)
        seq_index = Shove("file://" + base_output_filename + ".transcript_by_gene.idx", "memory://", optimize=False)
        transcripts = seq_index["CP"]

        self.assertTrue(len(transcripts) == 15)
        troubled_transcript = "ENST00000474204.1"
        is_troubled_transcript_seen = False
        for tx in transcripts:
            if tx.get_transcript_id() == troubled_transcript:
                is_troubled_transcript_seen = True
                break
        self.assertTrue(is_troubled_transcript_seen)

    def test_multiple_gtf_initialization(self):
        """Test that we can create a datasource from multiple gtf & fastas"""
        gencode_input_gtfs = ["testdata/gencode/CP.gencode.annotation.gtf", "testdata/gencode/MAPK1.gencode.v18.annotation.gtf"]
        gencode_input_fastas = ["testdata/gencode/CP.gencode.pc_transcripts.fa", "testdata/gencode/MAPK1.gencode.v18.pc_transcripts.fa"]
        base_output_filename = "out/test_multi_gencode"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(gencode_input_gtfs, gencode_input_fastas, base_output_filename)
        seq_index = Shove("file://" + base_output_filename + ".transcript_by_gene.idx", "memory://", optimize=False)
        transcripts = seq_index["CP"]
        self.assertTrue(len(transcripts) == 15)
        transcripts = seq_index["MAPK1"]
        self.assertTrue(len(transcripts) == 4)
        for tx in transcripts:
            self.assertTrue(tx.get_transcript_id() == "ENST00000491588.1" or len(tx.get_seq()) > 100, "No seq data for " + tx.get_transcript_id() )
