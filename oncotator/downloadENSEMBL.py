# LICENSE_GOES_HERE


from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os
import pprint
__author__ = 'lichtens'
from oncotator.utils.install.GenomeBuildInstallUtils import GenomeBuildInstallUtils
from oncotator.utils.install.GenomeBuildInstallUtils import VALID_ENSEMBL_SPECIES

def parseOptions():
    # Setup argument parser
    epilog= '''
    Supported species:

    ''' + "\n\t".join(VALID_ENSEMBL_SPECIES)

    desc = ''' Download ENSEMBL files from FTP server.  Please see -h for supported species list.
    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument('speciesName', type=str, help="Species name")
    parser.add_argument('downloadDir', type=str, help="Where to download files.")

    # Process arguments
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parseOptions()
    species = args.speciesName
    destDir = args.downloadDir

    if not os.path.exists(destDir):
        os.makedirs(destDir)
    GenomeBuildInstallUtils.download_reference_data_from_ensembl(destDir, species)