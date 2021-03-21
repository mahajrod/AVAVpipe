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

parser.add_argument("-f", "--fai", action="store", dest="fai", required=True,
                    help=".fai file ")
parser.add_argument("-w", "--scaffold_whitelist", action="store", dest="scaffold_whitelist", required=True,
                    type=lambda path: read_csv(path, sep="\t", squeeze=True, header=None),
                    help="File or comma-separated list with scaffolds from white list")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file with intervals")

args = parser.parse_args()

with metaopen(args.fai, "r") as in_fd, metaopen(args.output, "w") as out_fd:
    for line in in_fd:
        scaf, end = line.split("\t")[:2]
        if scaf in args.scaffold_whitelist:
            out_fd.write("{0}:1-{1}\n".format(scaf, end))
