#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import sys
import argparse

from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd


def parse_file_to_series(path):
    if Path(path).is_file():
        try:
            return pd.read_csv(path, sep="\t", squeeze=True, header=None)
        except:
            return pd.Series(dtype=str)
    else:
        return pd.Series(path.split(","))


def safe_mkdir(path):
    Path(path).mkdir(exist_ok=True)


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fai", action="store", dest="fai", required=True,
                    help=".fai file ")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    help="Output directory. If not set output will be written to stdout")
parser.add_argument("-s", "--split_scaffolds", action="store_true", dest="split_scaffolds", default=False,
                    help="Split scaffolds. Default: False")
parser.add_argument("-m", "--max_length", action="store", dest="max_length", type=int,
                    help="Soft maximum length of region(1.5x longer regions are allowed). Default: not set")
parser.add_argument("-n", "--max_seq_number", action="store", dest="max_seq_number", type=int, default=1,
                    help="Maximum number of sequences per region. Default: 1")
parser.add_argument("-b", "--scaffold_blacklist", action="store", dest="scaffold_blacklist",
                    default=pd.Series(dtype=str),
                    type=parse_file_to_series,
                    help="File or comma-separated list with scaffolds from black list")
parser.add_argument("-w", "--scaffold_whitelist", action="store", dest="scaffold_whitelist",
                    default=pd.Series(dtype=str),
                    type=parse_file_to_series,
                    help="File or comma-separated list with scaffolds from white list")
parser.add_argument("-x", "--min_scaffold_len", action="store", dest="min_scaffold_len", type=int, default=None,
                    help="Minimum length of scaffold to be included in regions. Default: not set")
parser.add_argument("-g", "--region_file_format", action="store", dest="region_file_format", default='simple',
                    help="Output region file format. "
                         "Allowed: 'simple' (default, not appliable for stdout), 'GATK', 'samtools'")
args = parser.parse_args()


def prepare_region_list_by_length(fai_file, max_length=500000, max_seq_number=10,
                                  output_dir=None,
                                  split_scaffolds=True, min_scaffold_length=None, blacklist_scaffolds=None,
                                  whitelist_scaffolds=None,
                                  region_file_format='simple'):

    len_df = pd.read_csv(fai_file, sep="\t", usecols=[0, 1], index_col=0, names=["scaffold", "length"])

    if min_scaffold_length:
        len_df = len_df[len_df["length"] >= min_scaffold_length]

    if not blacklist_scaffolds.empty:
        len_df = len_df[~len_df.index.isin(blacklist_scaffolds)]
    if not whitelist_scaffolds.empty:
        len_df = len_df[len_df.index.isin(whitelist_scaffolds)]
    number_of_scaffolds = len(len_df)

    max_length_soft_threshold = None if max_length is None else int(1.5 * max_length)
    region_list = []

    remnant_seq_list = []
    remnant_seq_length = 0

    scaffold_to_region_correspondence_dict = OrderedDict()
    region_index = 0

    if max_length_soft_threshold is None:
        key_list = list(len_df.index)
        bins = np.arange(0, number_of_scaffolds, max_seq_number)
        bins = bins if bins[-1] == number_of_scaffolds else np.append(bins, number_of_scaffolds)

        for i in range(0, len(bins) - 1):
            bunch_list = []
            for region in key_list[bins[i]:bins[i + 1]]:
                #bunch_list.append([region, 1, len_dict[region]])
                bunch_list.append([region, 1, len_df.loc[region, "length"]])
            region_list.append(bunch_list)
            for scaffold in key_list[bins[i]:bins[i + 1]]:
                scaffold_to_region_correspondence_dict[scaffold] = [i]

    elif not split_scaffolds:
        bunch_length = 0
        bunch_list = []

        for region in len_df.index:
            if len_df.loc[region, "length"] >= max_length:
                region_list.append([[region, 1, len_df.loc[region, "length"]]])
                scaffold_to_region_correspondence_dict[region] = [region_index]
                region_index += 1
            else:
                bunch_list.append([region, 1, len_df.loc[region, "length"]])
                bunch_length += len_df.loc[region, "length"]
                if (bunch_length >= max_length) or (len(bunch_list) == max_seq_number):
                    region_list.append(bunch_list)
                    for scaffold in bunch_list:
                        scaffold_to_region_correspondence_dict[scaffold[0]] = [region_index]
                    region_index += 1
                    bunch_length = 0
                    bunch_list = []
        if bunch_list:
            region_list.append(bunch_list)
            for scaffold in bunch_list:
                scaffold_to_region_correspondence_dict[scaffold[0]] = [region_index]
            region_index += 1
            bunch_length = 0
            bunch_list = []
    else:
        remnant_seq_length = 0
        for region in len_df.index:
            if len(remnant_seq_list) == max_seq_number:
                region_list.append(remnant_seq_list)
                for remnant_entry in remnant_seq_list:
                    if remnant_entry[0] in scaffold_to_region_correspondence_dict:
                        scaffold_to_region_correspondence_dict[remnant_entry[0]].append(region_index)
                    else:
                        scaffold_to_region_correspondence_dict[remnant_entry[0]] = [region_index]
                region_index += 1
                remnant_seq_list = []
                remnant_seq_length = 0
            #print(len_df)
            #print(len_df.loc[region, "length"])
            if len_df.loc[region, "length"] > max_length:
                points = np.arange(0, len_df.loc[region, "length"], max_length)
                if len(points) > 1:
                    scaffold_to_region_correspondence_dict[region] = list(
                        np.arange(region_index, region_index + len(points) - 1))
                    region_index += len(points) - 1
                    for i in range(0, len(points) - 1):
                        region_list.append([[region, points[i] + 1, points[i + 1]]])

                    remnant_length = len_df.loc[region, "length"] - points[-1]
                    if remnant_length + max_length <= max_length_soft_threshold:
                        region_list[-1][0][2] = len_df.loc[region, "length"]
                        remnant_length = 0
                    else:
                        region_list.append([[region, points[-1] + 1, len_df.loc[region, "length"]]])
                        scaffold_to_region_correspondence_dict[region].append(region_index)
                        region_index += 1
                else:

                    region_list.append([[region, 1, len_df.loc[region, "length"]]])
                    scaffold_to_region_correspondence_dict[region] = [region_index]
                    region_index += 1
                continue
            else:
                remnant = [region, 1, len_df.loc[region, "length"]]
                remnant_length = len_df.loc[region, "length"]

            if remnant_seq_length + remnant_length <= max_length_soft_threshold:
                remnant_seq_list.append(remnant)
                remnant_seq_length += remnant_length
            else:
                region_list.append(remnant_seq_list)
                for remnant_entry in remnant_seq_list:
                    if remnant_entry[0] in scaffold_to_region_correspondence_dict:
                        scaffold_to_region_correspondence_dict[remnant_entry[0]].append(region_index)
                    else:
                        scaffold_to_region_correspondence_dict[remnant_entry[0]] = [region_index]
                region_index += 1
                remnant_seq_list = [remnant]
                remnant_seq_length = remnant_length
        else:
            if remnant_seq_list:
                region_list.append(remnant_seq_list)
                for remnant_entry in remnant_seq_list:
                    if remnant_entry[0] in scaffold_to_region_correspondence_dict:
                        scaffold_to_region_correspondence_dict[remnant_entry[0]].append(region_index)
                    else:
                        scaffold_to_region_correspondence_dict[remnant_entry[0]] = [region_index]
                region_index += 1

    if output_dir:
        for directory in output_dir, "%s/intervals/" % output_dir:
            safe_mkdir(directory)
        scaffold_ids = pd.Series(len_df.index)
        #print(scaffold_ids)
        scaffold_ids.to_csv("%s/scaffold.ids" % output_dir, sep="\t", header=False, index=False)
        len_df.to_csv("%s/scaffold.len" % output_dir, sep="\t", index=True, header=False)
        index = 1
        for regions in region_list:
            with open("%s/intervals/region_%i.list" % (output_dir, index), "w") as out_fd:
                for region in regions:
                    if isinstance(region, str):
                        out_fd.write(region)
                        out_fd.write("\n")
                    else:
                        if region_file_format == 'simple':
                            if len(region) == 3:
                                out_fd.write("{0}\t{1}\t{2}\n".format(*region))
                            elif len(region) == 1:
                                out_fd.write(region[0])
                                out_fd.write("\n")
                        elif region_file_format == 'GATK':
                            if len(region) == 3:
                                out_fd.write("{0}:{1}-{2}".format(*region))
                            elif len(region) == 1:
                                out_fd.write(region[0])
                                out_fd.write("\n")
                        elif region_file_format == 'samtools':
                            if len(region) == 3:
                                out_fd.write("\t".join(list(map(str, region))))
                            elif len(region) == 1:
                                out_fd.write(region[0])
                                out_fd.write("\n")
            index += 1
        with open("%s/scaffold_to_region.correspondence" % output_dir, "w") as cor_fd:
            for region in scaffold_to_region_correspondence_dict:
                cor_fd.write("{0}\t{1}\n".format(region, ",".join(map(str, scaffold_to_region_correspondence_dict[region]))))
        #scaffold_to_region_correspondence_dict.write("%s/SCAFFOLD_TO_REGION.correspondence" % output_dir,
        #                                             splited_values=True)
    else:
        if region_file_format == "samtools":
            for regions in region_list:
                sys.stdout.write(",".join(["{0}:{1}-{2}".format(*region) if len(region) == 3 else region[0] for region in regions]) + "\n")
        elif region_file_format == "GATK":
            for regions in region_list:
                sys.stdout.write(" -L " + " -L ".join(["{0}:{1}-{2}".format(*region) if len(region) == 3 else region[0]  for region in regions]) + "\n")
    return region_list, scaffold_to_region_correspondence_dict


prepare_region_list_by_length(args.fai, max_length=args.max_length, max_seq_number=args.max_seq_number,
                              output_dir=args.output_dir,
                              split_scaffolds=args.split_scaffolds, min_scaffold_length=args.min_scaffold_len,
                              blacklist_scaffolds=args.scaffold_blacklist,
                              whitelist_scaffolds=args.scaffold_whitelist,
                              region_file_format=args.region_file_format)
