import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Given a directory, this script creates a concatenated file that "
                                                 "consists .")

    parser.add_argument("--dir", action="store", type=str, dest="dir", required=True,
                        help="Directory where the tsv files are located")
    parser.add_argument("--prefix", action="store", type=str, dest="prefix", required=True,
                        help="Prefix for the tsv files")
    parser.add_argument("--out", action="store", type=str, dest="out", required=True,
                        help="Name of the concatenated output file")

    args, _ = parser.parse_known_args()

    # Collect the names of all tsv files with the prefix


    # Concatenate files



if __name__ == "__main__":
    main()