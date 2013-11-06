"""
# By downloading the PROGRAM you agree to the following terms of use:
#
# BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
# FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
#
# This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
# WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
# WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
# NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
#
# 1. DEFINITIONS
# 1.1	"PROGRAM" shall mean copyright in the object code and source code known as Oncotator and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/cancer/cga/oncotator on the EFFECTIVE DATE.
#
# 2. LICENSE
# 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
# LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
# The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
# 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
# 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
#
# 3. OWNERSHIP OF INTELLECTUAL PROPERTY
# LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
#
# Copyright 2012 Broad Institute, Inc.
# Notice of attribution:  The Oncotator program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc.
#
# LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
#
# 4. INDEMNIFICATION
# LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
#
# 5. NO REPRESENTATIONS OR WARRANTIES
# THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
# IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
#
# 6. ASSIGNMENT
# This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
#
# 7. MISCELLANEOUS
# 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
# 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
# 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
# 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
# 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
# 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
# 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
#"""
import os
import shutil
import unittest
from shove.core import Shove
from oncotator.utils.install.GenomeBuildFactory import GenomeBuildFactory
from oncotator.index.gaf import region2bin,region2bins

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
        genome_build_factory.build_ensembl_transcript_index(ensembl_input_gtf, ensembl_input_fasta, output_filename, protocol=protocol)
        self.assertTrue(os.path.exists(output_filename))

        shove = Shove(protocol + "://" + output_filename, "memory://")
        self.assertTrue(len(shove.keys()) > 0)
        self.assertTrue("YDR529C" in shove.keys())
        t = shove["YDR529C"]
        self.assertTrue(t.get_seq() is not None)
        self.assertTrue(t.get_seq() is not "")
        self.assertTrue(len(t.get_cds()) > 0)
        self.assertTrue(len(t.get_exons()) > 0)
        shutil.rmtree(output_filename)

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
        genome_build_factory.build_ensembl_transcript_index(ensembl_input_gtf, ensembl_input_fasta, transcript_index_filename, protocol=protocol)
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
        genome_build_factory.build_ensembl_transcript_index(ensembl_input_gtf, ensembl_input_fasta, transcript_index_filename, protocol=protocol)
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
        genome_build_factory.construct_ensembl_indices(ensembl_input_gtf, ensembl_input_fasta, base_output_filename)

        self.assertTrue(os.path.exists(base_output_filename + ".transcript.idx"))
        self.assertTrue(os.path.exists(base_output_filename + ".transcript_by_gene.idx"))
        self.assertTrue(os.path.exists(base_output_filename + ".transcript_by_gp_bin.idx"))

    def test_retrieving_sequence(self):
        """Ensure we can retrieve a sequence from an ensembl transcript given a gene.  """

        ensembl_input_gtf = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.gtf"
        ensembl_input_fasta = "testdata/Saccharomyces_cerevisiae.EF4.71_trim.cdna.all.fa"
        base_output_filename = "out/test_full_indices_ensembl"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp.idx", ignore_errors=True)
        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(ensembl_input_gtf, ensembl_input_fasta, base_output_filename)

        seq_index = Shove("file://" + base_output_filename + ".transcript_by_gene.idx", "memory://", optimize=False)
        transcripts = seq_index['SEO1']

        self.assertTrue(transcripts[0].get_seq().startswith('ATGTATTCAATTGTTAAAGAGATTATTGTAGATCCTTACAAAAGACTAAAATGGGGTTTT'))

        transcripts = seq_index['PAU8']
        self.assertTrue(transcripts[0].get_strand() == "-")

        seq_index_gp = Shove("file://" + base_output_filename + ".transcript_by_gp_bin.idx", "memory://")
        transcripts = seq_index_gp["I_585"]
        self.assertTrue(transcripts[0].get_strand() == "+")

    def test_gencode_small(self):
        """Test that we can create Transcript instances from a small gencode gtf and fasta."""
        gencode_input_gtf = "testdata/gencode/MAPK1.gencode.v18.annotation.gtf"
        gencode_input_fasta = "testdata/gencode/MAPK1.gencode.v18.pc_transcripts.fa"
        base_output_filename = "out/test_small_gencode"
        shutil.rmtree(base_output_filename + ".transcript.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gene.idx", ignore_errors=True)
        shutil.rmtree(base_output_filename + ".transcript_by_gp_bin.idx", ignore_errors=True)

        genome_build_factory = GenomeBuildFactory()
        genome_build_factory.construct_ensembl_indices(gencode_input_gtf, gencode_input_fasta, base_output_filename)

        seq_index = Shove("file://" + base_output_filename + ".transcript_by_gene.idx", "memory://", optimize=False)
        transcripts = seq_index["MAPK1"]
        self.assertTrue(len(transcripts) == 4)

        seq_index_gp = Shove("file://" + base_output_filename + ".transcript_by_gp_bin.idx", "memory://", optimize=False)
        transcripts = seq_index_gp["chr22_753"]
        self.assertTrue(transcripts[0].get_strand() == "-")
        self.assertTrue(len(transcripts) == 4)
