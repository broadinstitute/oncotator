from argparse import ArgumentParser, RawTextHelpFormatter
from oncotator.utils.install.DatasourceInstallUtils import DatasourceInstallUtils


def parseOptions():
    desc = """
    Create a md5 hash within a datasource directory.  The genome direcotory should be specified.
    For example:  oncotator_add_DatasourceMd5 /home/user/db_dir/gaf/hg19
    This will create a gaf/hg19.md5 file that is based on every file (incl. subdirectories) in the gaf/hg19/.

    """

    epilog=""" """
    parser = ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter, epilog=epilog)
    parser.add_argument("ds_dir", type=str, help="datasource directory.  This should be the genome_build dir inside a datasource dir.")

    # Process arguments
    args = parser.parse_args()
    return args

def main():
    args = parseOptions()
    ds_dir = args.ds_dir
    DatasourceInstallUtils.create_datasource_md5_file(ds_dir)

if __name__ == '__main__':
    main()