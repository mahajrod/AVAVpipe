#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import bz2
import gzip
import sys

import argparse

from copy import deepcopy
from pandas import read_csv, Series
from collections import OrderedDict


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


def check_gvcf_integrity(gvcf_file, reference_fai, output_prefix):
    len_df = read_csv(reference_fai, sep="\t", header=None, usecols=[0, 1], names=["scaffold", "length"], index_col=0)

    scaffold_dict = OrderedDict()

    with metaopen(gvcf_file, "r") as gvcf_fd:
        prev_scaffold = ""

        for line in gvcf_fd:
            # print line
            if line[0] == "#":
                continue

            line_list = line.split("\t")
            scaffold = line_list[0]
            start = int(line_list[1])
            format = line_list[7].split(";")

            if (len(format) == 1) and (format[0][0:3] == "END"):
                end = int(format[0].split("=")[1])
            else:
                end = start + len(line_list[3]) - 1
            # print line_list
            # print scaffold, start, end, format

            if scaffold not in scaffold_dict:
                scaffold_dict[scaffold] = []

            if scaffold != prev_scaffold:
                scaffold_dict[scaffold].append([deepcopy(start), deepcopy(end)])
            else:
                # print scaffold_dict[scaffold][-1][1]
                if scaffold_dict[scaffold][-1][1] + 1 >= start:
                    scaffold_dict[scaffold][-1][1] = deepcopy(max(end, scaffold_dict[scaffold][-1][1]))
                else:
                    print(scaffold_dict[scaffold])
                    print(line)
                    scaffold_dict[scaffold].append([deepcopy(start), deepcopy(end)])
            prev_scaffold = scaffold

    complete_scaffolds = []
    fragmented_scaffolds = []
    scaffolds_with_absent_fragments = []

    with open("%s.scaffold_regions" % output_prefix, "w") as scaf_reg_fd:

        for scaffold in scaffold_dict:
            if len(scaffold_dict[scaffold]) > 1:
                fragmented_scaffolds.append(scaffold)

            scaffold_length = sum(map(lambda s: s[1] - s[0] + 1, scaffold_dict[scaffold]))
            if scaffold_length != len_df.loc[scaffold, "length"]:
                scaffolds_with_absent_fragments.append(scaffold)
            else:
                complete_scaffolds.append(scaffold)
            scaf_reg_fd.write(
                "%s\t%s\n" % (scaffold, ",".join(map(lambda s: "-".join(map(str, s)), scaffold_dict[scaffold]))))

    Series(complete_scaffolds).to_csv("%s.complete_scaffolds" % output_prefix, sep="\t", header=False, index=False)
    Series(fragmented_scaffolds).to_csv("%s.fragmented_scaffolds" % output_prefix, sep="\t", header=False, index=False)
    Series(scaffolds_with_absent_fragments).to_csv("%s.scaffolds_with_absent_fragments" % output_prefix, sep="\t", header=False, index=False)


parser = argparse.ArgumentParser()

parser.add_argument("-g", "--gvcf", action="store", dest="gvcf", required=True,
                    help="Input gvcf file")
parser.add_argument("-f", "--fai", action="store", dest="fai", required=True,
                    help="Reference .fai file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")

args = parser.parse_args()

check_gvcf_integrity(args.gvcf, args.fai, args.output_prefix)
