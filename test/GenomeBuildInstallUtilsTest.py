# LICENSE_GOES_HERE
import os
import unittest
from oncotator.utils.install.GenomeBuildInstallUtils import GenomeBuildInstallUtils
from TestUtils import TestUtils
from oncotator.utils.MutUtils import MutUtils

TestUtils.setupLogging(__file__, __name__)


class GenomeBuildInstallUtilsTest(unittest.TestCase):
    _multiprocess_can_split_ = True

    def test_current_download(self):
        """Download a current ensembl transcript package.  This test needs an internet connection and can be slow."""

        #ftp://ftp.ensembl.org/pub/release-71/fasta/saccharomyces_cerevisiae/cdna/
        download_dir = "out/test_ensembl_download/"
        MutUtils.removeDir(download_dir)
        os.mkdir(download_dir)
        GenomeBuildInstallUtils.download_reference_data_from_ensembl(download_dir, "saccharomyces_cerevisiae")

        downloaded_files = os.listdir(download_dir)

        transcript_file = None
        for f in downloaded_files:
            if f.find(".cdna.") != -1:
                transcript_file = f
                break

        self.assertIsNotNone(transcript_file)

        statinfo = os.stat(download_dir + transcript_file)
        self.assertTrue(statinfo.st_size > 0, "downloaded transcript file (" + transcript_file + ") is empty.")

    def test_previous_release_download(self):
        """Download an older ensembl transcript package.  This test needs an internet connection and will fail w/o one.
        """
        download_dir = "out/test_ensembl_download_previous/"
        release_num = "68"
        MutUtils.removeDir(download_dir)
        os.mkdir(download_dir)
        GenomeBuildInstallUtils.download_reference_data_from_ensembl(download_dir, "saccharomyces_cerevisiae", release=release_num)

        downloaded_files = os.listdir(download_dir)

        transcript_file = None
        for f in downloaded_files:
            if f.find("." + release_num + ".cdna.") != -1:
                transcript_file = f
                break

        self.assertIsNotNone(transcript_file)

        statinfo = os.stat(download_dir + transcript_file)
        self.assertTrue(statinfo.st_size > 0, "downloaded transcript file (" + transcript_file + ") is empty.")

    @unittest.skip("Backing code not implemented yet.")
    def test_index_ensembl_files(self):
        """Test the instantiation of some files for creating a datasource file."""
        output_dir = "out/"
        GenomeBuildInstallUtils.create_ensembl_transcript_datasource(ensembl_species="saccharomyces_cerevisiae", genome_build="sacCer3")

        self.assertTrue(os.path.exists(output_dir + "ensembl/sacCer3"))
        self.assertTrue(os.path.exists(output_dir + "ensembl/ensembl.config"))
        statinfo = os.stat(output_dir + "ensembl/ensembl.config")
        self.assertTrue(statinfo.st_size > 0, "generated config file (" + output_dir + "ensembl/ensembl.config) is empty.")



if __name__ == '__main__':
    unittest.main()
