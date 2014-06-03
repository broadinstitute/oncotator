import argparse
import os
import csv
import string
import natsort


def main():
    parser = argparse.ArgumentParser(description="Concatenate TSVs in a given directory. Consistent column ordering, "
                                                 "commenting and header location are assumed.")

    parser.add_argument("--dir", action="store", type=str, dest="dir", required=True,
                        help="Directory where the tsv files are located")
    parser.add_argument("--prefix", action="store", type=str, dest="prefix", required=True,
                        help="Prefix for the tsv files")
    parser.add_argument("--out", action="store", type=str, dest="out",
                        help="Name of the concatenated output file")

    args, _ = parser.parse_known_args()

    # Collect the names of all tsv files with the prefix
    filenames = [os.path.join(args.dir, entry) for entry in os.listdir(args.dir) if entry.startswith(args.prefix)]
    filenames = natsort.natsorted([filename for filename in filenames if os.path.isfile(filename)])

    chunk_size = 100000
    filename = args.out if args.out is not None else string.join([args.prefix, ".tsv"], ".")
    with open(filename, "w") as outfile:  # Iterate through files and concatenate them
        writer = None
        for filename in filenames:
            with open(filename, "r") as infile:
                reader = csv.DictReader(infile, delimiter="\t")
                if writer is None:
                    writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=reader.fieldnames)
                    writer.writeheader()
                rowdicts = []
                num_rows = 1
                for rowdict in reader:
                    rowdicts += [rowdict]
                    if num_rows % chunk_size == 0:
                        writer.writerows(rowdicts)
                        rowdicts = []
                    num_rows += 1

                if num_rows != 0:
                    writer.writerows(rowdicts)

if __name__ == "__main__":
    main()