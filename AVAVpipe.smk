import os

#conda: "%s/config/conda.yaml" % __file__
#configfile: "%s/config/default.yaml" % __file__

fastqc_dir = "{0}/{1}".format(config["out_dir"], config["fastqc_dir"])
kmer_dir = "{0}/{1}".format(config["out_dir"], config["kmer_dir"])
filtered_read_dir = "{0}/{1}".format(config["out_dir"], config["filtered_read_dir"])
alignment_dir = "{0}/{1}".format(config["out_dir"], config["alignment_dir"])
snpcall_dir = "{0}/{1}".format(config["out_dir"], config["snpcall_dir"])
log_dir = "{0}/{1}".format(config["out_dir"], config["log_dir"])
error_dir  = "{0}/{1}".format(config["out_dir"], config["error_dir"])
benchmark_dir =  "{0}/{1}".format(config["out_dir"], config["benchmark_dir"])

# if "sample_list" key is absent in config variable, use folder names from config["sample_dir"] as sample ids
if "sample_list" not in config:
    samples = []
    for entry in os.scandir(config["sample_dir"]):
        if entry.is_dir():
            samples.append(entry.name)

    config["sample_list"] = samples

localrules: all

rule all:
    input:
        expand("%s/{sample_id}/fastqc_raw.log" % log_dir, sample_id=config["sample_list"]),
        expand("%s/{sample_id}/fastqc_filtered.log" % log_dir, sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.trimmed_1.fastq.gz" % filtered_read_dir, sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.trimmed_2.fastq.gz" % filtered_read_dir, sample_id=config["sample_list"]),
        expand("%s/{sample_id}/{sample_id}.sorted.mkdup.bam" % alignment_dir, sample_id=config["sample_list"]),
        expand("%s/{sample_id}/{sample_id}.sorted.mkdup.bam.bai" % alignment_dir, sample_id=config["sample_list"]),
        expand("%s/{sample_id}/{sample_id}.coverage.per-base.bed.gz" % alignment_dir, sample_id=config["sample_list"]),
        rules.jellyfish_histo.output


include: "workflow/QCFiltering/FastQC_raw.smk"
include: "workflow/QCFiltering/Trimmomatic.smk"
include: "workflow/QCFiltering/FastQC_filtered.smk"
include: "workflow/QCFiltering/Kmer.smk"
#include: "workflow/Alignment/Alignment.smk"
#include: "workflow/Alignment/Coverage.smk"
#include: "workflow/VariantCall/BSQR.smk"
#include: "workflow/VariantCall/Genotyping.smk"

