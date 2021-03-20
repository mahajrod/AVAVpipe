#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import bz2
import gzip
import sys

import argparse

from pandas import read_csv


if sys.version_info[0] == 3:
    from io import TextIOWrapper as file


def metaopen(filename, flags, buffering=None):
    if not isinstance(filename, str):  # or isinstance(filename, gzip.GzipFile) or isinstance(filename, bz2.BZ2File):
        if isinstance(filename, file):
            return filename
        else:
            raise ValueError("ERROR!!! Not file object or str: {}".format(str(filename)))
    elif filename[-3:] == ".gz":
        return gzip.open(filename, flags + ("t" if "b" not in flags else ""))
    elif filename[-4:] == ".bz2":
        return bz2.open(filename, flags + ("t" if "b" not in flags else ""))
    else:
        if buffering is not None:
            return open(filename, flags, buffering=buffering)
        else:
            return open(filename, flags)


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--path_file", action="store", dest="path_file", required=True,
                    type=lambda path_file: read_csv(path_file, sep="\t", squeeze=True, header=None),
                    help="File with list of paths to vcf files (one path per line)")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file with combined data")

args = parser.parse_args()

print("Handling file %s ..." % args.path_file.iloc[0])
with metaopen(args.output, "w") as out_fd:

    with metaopen(args.path_file.iloc[0], "r") as in_fd:
        for line in in_fd:
            out_fd.write(line)

    if len(args.path_file) > 1:
        for entry in args.path_file.iloc[1:, ]:
            print("Handling file %s ..." % entry)
            with metaopen(entry, "r") as in_fd:
                for line in in_fd:
                    if line[0] == "#":
                        continue
                    else:
                        out_fd.write(line)
                        break
                for line in in_fd:
                    out_fd.write(line)
