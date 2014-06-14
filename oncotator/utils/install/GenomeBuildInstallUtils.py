# LICENSE_GOES_HERE
import logging

import subprocess
import os
import ftplib
import re
import tarfile
import subprocess
from shove.core import Shove
from oncotator.Transcript import Transcript
from oncotator.utils.GenericTsvReader import GenericTsvReader

VALID_ENSEMBL_SPECIES = [
'ailuropoda_melanoleuca',
'anolis_carolinensis',
'bos_taurus',
'caenorhabditis_elegans',
'callithrix_jacchus',
'canis_familiaris',
'cavia_porcellus',
'choloepus_hoffmanni',
'ciona_intestinalis',
'ciona_savignyi',
'danio_rerio',
'dasypus_novemcinctus',
'dipodomys_ordii',
'drosophila_melanogaster',
'echinops_telfairi',
'equus_caballus',
'erinaceus_europaeus',
'felis_catus',
'gadus_morhua',
'gallus_gallus',
'gasterosteus_aculeatus',
'gorilla_gorilla',
'homo_sapiens',
'loxodonta_africana',
'macaca_mulatta',
'macropus_eugenii',
'meleagris_gallopavo',
'microcebus_murinus',
'monodelphis_domestica',
'mus_musculus',
'myotis_lucifugus',
'nomascus_leucogenys',
'ochotona_princeps',
'ornithorhynchus_anatinus',
'oryctolagus_cuniculus',
'oryzias_latipes',
'otolemur_garnettii',
'pan_troglodytes',
'petromyzon_marinus',
'pongo_abelii',
'procavia_capensis',
'pteropus_vampyrus',
'rattus_norvegicus',
'saccharomyces_cerevisiae',
'sarcophilus_harrisii',
'sorex_araneus',
'spermophilus_tridecemlineatus',
'sus_scrofa',
'taeniopygia_guttata',
'takifugu_rubripes',
'tarsius_syrichta',
'tetraodon_nigroviridis',
'tupaia_belangeri',
'tursiops_truncatus',
'vicugna_pacos',
'xenopus_tropicalis'
]


class GenomeBuildInstallUtils(object):

    @staticmethod
    def download_reference_data_from_ensembl(dl_dir, species, release=""):
        """Download reference GTF file and transcript FASTA for given species from Ensembl.
        If release is empty string, get the latest (current).  Release must be a number as a string e.g. "71" """
        if species not in VALID_ENSEMBL_SPECIES:
            raise Exception('Not a valid Ensembl species ID!')


        release_gtf_dir = "/pub/current_gtf/"
        release_fasta_dir = '/pub/current_fasta/'
        if (release is not None) and (release.strip() != ""):
            release_gtf_dir = "/pub/release-" + release + "/gtf/"
            release_fasta_dir = "/pub/release-" + release + "/fasta/"

        ensembl_ftp_url = 'ftp.ensembl.org'

        ftp = ftplib.FTP(ensembl_ftp_url)
        ftp.connect()
        ftp.login()

        ftp_data_dir = os.path.join(release_gtf_dir, species)

        ftp.cwd(ftp_data_dir)
        fnames = ftp.nlst()

        gtf_file = [f for f in fnames if f.endswith('.gtf.gz')][0]
        gtf_version = gtf_file.partition('.gtf.gz')[0]

        print "Now Downloading {0}...".format(gtf_file)
        gtf_fname = os.path.join(dl_dir, gtf_file)
        ftp.retrbinary('RETR ' + gtf_file, open(gtf_fname, 'wb').write)
        print "Done!"

        fnames_list = list()
        fnames_list.append(gtf_fname)

        for fasta_type in ('cdna', 'ncrna', 'pep'):
            ftp_data_dir = os.path.join(release_fasta_dir, species, fasta_type )
            ftp.cwd(ftp_data_dir)
            fnames = ftp.nlst()

            if fasta_type == 'ncrna':
                ext = fasta_type + '.fa.gz'
            else:
                ext = fasta_type + '.all.fa.gz'

            fasta_file = [f for f in fnames if f.endswith(ext)][0]
            fasta_version = fasta_file.partition(ext)[0][:-1]

            if fasta_version != gtf_version:
                raise Exception('Fasta file version does not match gtf file version!  ' + fasta_version + ' and ' + gtf_version)

            print "Now Downloading {0}...".format(fasta_file)
            fasta_fname = os.path.join(dl_dir, fasta_file)
            ftp.retrbinary('RETR ' + fasta_file, open(fasta_fname, 'wb').write)
            print "Done!"

            fnames_list.append(fasta_fname)

        ftp.close()

        return gtf_version

