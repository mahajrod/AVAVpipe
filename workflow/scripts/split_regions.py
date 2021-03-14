def prepare_region_list_by_length(self, max_length=500000, max_seq_number=10,
                                  length_dict=None, reference=None, parsing_mode="parse", output_dir=None,
                                  split_scaffolds=True, min_scaffold_length=None, black_list_scaffolds=None,
                                  white_list_scaffolds=None,
                                  region_file_format='simple'):
    if reference[-4:] == ".fai":
        # len_df = pd.read_csv(sep="\t", filename=reference, use_cols=(0, 1), index_col=0,
        #                     col_names=("scaffold", "length"), header=None)
        raw_len_dict = SynDict(filename=reference, key_index=0, value_index=1, expression=int)
    else:
        raw_len_dict = length_dict if length_dict else self.get_lengths(record_dict=self.parse_seq_file(reference,
                                                                                                        mode=parsing_mode),
                                                                        out_file=None,
                                                                        close_after_if_file_object=False)

        # len_df = pd.DataFrame.from_dict(raw_len_dict, orient="index", columns=("scaffold", "length"))

    len_dict = OrderedDict()
    """
    if black_list_scaffolds:
        len_df = len_df[~len_df.index.isin(black_list_scaffolds)]

    if white_list_scaffolds:
        len_df = len_df[len_df.index.isin(white_list_scaffolds)]

    if min_scaffold_length:
        len_df = len_df[len_df["length"] >= min_scaffold_length]

    number_of_scaffolds = len(len_df)
    """
    if black_list_scaffolds and min_scaffold_length:
        for scaffold_id in raw_len_dict:
            if (scaffold_id in black_list_scaffolds) or (raw_len_dict[scaffold_id] < min_scaffold_length):
                continue
            len_dict[scaffold_id] = raw_len_dict[scaffold_id]

    elif black_list_scaffolds:
        for scaffold_id in raw_len_dict:
            if scaffold_id in black_list_scaffolds:
                continue
            len_dict[scaffold_id] = raw_len_dict[scaffold_id]

    elif min_scaffold_length:
        for scaffold_id in raw_len_dict:
            if raw_len_dict[scaffold_id] < min_scaffold_length:
                continue
            len_dict[scaffold_id] = raw_len_dict[scaffold_id]
    else:
        len_dict = raw_len_dict

    number_of_scaffolds = len(len_dict)

    if white_list_scaffolds:
        tmp_dict = OrderedDict()
        for scaffold in len_dict:
            if scaffold in white_list_scaffolds:
                tmp_dict[scaffold] = len_dict[scaffold]
        len_dict = tmp_dict

    max_length_soft_threshold = None if max_length is None else int(1.5 * max_length)
    region_list = []

    remnant_seq_list = []
    remnant_seq_length = 0

    scaffold_to_region_correspondence_dict = SynDict()
    region_index = 0

    if max_length_soft_threshold is None:
        key_list = list(len_dict.keys())
        bins = np.arange(0, number_of_scaffolds, max_seq_number)
        bins = bins if bins[-1] == number_of_scaffolds else np.append(bins, number_of_scaffolds)

        for i in range(0, len(bins) - 1):
            bunch_list = []
            for region in key_list[bins[i]:bins[i + 1]]:
                bunch_list.append([region, 1, len_dict[region]])
            region_list.append(bunch_list)
            for scaffold in key_list[bins[i]:bins[i + 1]]:
                scaffold_to_region_correspondence_dict[scaffold] = [i]

    elif not split_scaffolds:
        bunch_length = 0
        bunch_list = []

        for region in len_dict:
            if len_dict[region] >= max_length:
                region_list.append([[region, 1, len_dict[region]]])
                scaffold_to_region_correspondence_dict[region] = [region_index]
                region_index += 1
            else:
                bunch_list.append([region, 1, len_dict[region]])
                bunch_length += len_dict[region]
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
        for region in len_dict:
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

            if len_dict[region] > max_length:
                points = np.arange(0, len_dict[region], max_length)
                if len(points) > 1:
                    scaffold_to_region_correspondence_dict[region] = list(
                        np.arange(region_index, region_index + len(points) - 1))
                    region_index += len(points) - 1
                    for i in range(0, len(points) - 1):
                        region_list.append([[region, points[i] + 1, points[i + 1]]])

                    remnant_length = len_dict[region] - points[-1]
                    if remnant_length + max_length <= max_length_soft_threshold:
                        region_list[-1][0][2] = len_dict[region]
                        remnant_length = 0
                    else:
                        region_list.append([[region, points[-1] + 1, len_dict[region]]])
                        scaffold_to_region_correspondence_dict[region].append(region_index)
                        region_index += 1
                else:

                    region_list.append([[region, 1, len_dict[region]]])
                    scaffold_to_region_correspondence_dict[region] = [region_index]
                    region_index += 1
                continue
            else:
                remnant = [region, 1, len_dict[region]]
                remnant_length = len_dict[region]

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
        for directory in output_dir, "%s/splited/" % output_dir:
            self.safe_mkdir(directory)
        scaffold_ids = IdList(len_dict)
        scaffold_ids.write("%s/scaffold.ids" % output_dir)
        len_dict.write("%s/scaffold.len" % output_dir)
        index = 1
        for regions in region_list:
            with open("%s/splited/region_%i.list" % (output_dir, index), "w") as out_fd:
                for region in regions:
                    if isinstance(region, str):
                        out_fd.write(region)
                        out_fd.write("\n")
                    else:
                        if region_file_format == 'simple':
                            if len(region) == 3:
                                out_fd.write("%s\t%s\t%s\n" % (region[0], region[1], region[2]))
                            elif len(region) == 1:
                                out_fd.write(region[0])
                                out_fd.write("\n")
                        elif region_file_format == 'GATK':
                            if len(region) == 3:
                                out_fd.write("%s:%s-%s\n" % (region[0], region[1], region[2]))
                            elif len(region) == 1:
                                out_fd.write(region[0])
                                out_fd.write("\n")
            index += 1
        scaffold_to_region_correspondence_dict.write("%s/SCAFFOLD_TO_REGION.correspondence" % output_dir,
                                                     splited_values=True)
    else:
        if region_file_format == "samtools":
            for regions in region_list:
                sys.stdout.write(",".join(["{0}:{1}-{2}".format(*region) for region in regions]) + "\n")
        elif region_file_format == "GATK":
            for regions in region_list:
                sys.stdout.write(" -L " + " -L ".join(["{0}:{1}-{2}".format(*region) for region in regions]) + "\n")
    return region_list, scaffold_to_region_correspondence_dict