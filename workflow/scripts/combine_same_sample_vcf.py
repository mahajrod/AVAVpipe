#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from pandas import read_csv


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--path_file", action="store", dest="path_file", required=True,
                    type=lambda path_file: read_csv(path_file, sep="\t", squeeze=True, header=None),
                    help="File with list of paths to vcf files (one path per line)")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file")

args = parser.parse_args()

print("Handling file %s ..." % args.path_file.iloc[0])
with open(args.output, "w") as out_fd:

    with open(args.path_file.iloc[0], "r") as in_fd:
        for line in in_fd:
            out_fd.write(line)

    if len(args.path_file) > 1:
        for filename in args.path_file.iloc[1:, ]:
            print("Handling file %s ..." % filename)
            with open(filename, "r") as in_fd:
                for line in in_fd:
                    if line[0] == "#":
                        continue
                    else:
                        out_fd.write(line)
                        break
                for line in in_fd:
                    out_fd.write(line)
