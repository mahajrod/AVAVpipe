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
        "{0}.fai".format(config["reference"]),
        "{0}.dict".format(os.path.splitext(config["reference"])[0]),
        expand("%s/{sample_id}/fastqc_raw.log" % log_dir, sample_id=config["sample_list"]),
        expand("%s/{sample_id}/fastqc_filtered.log" % log_dir, sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.sorted.mkdup.bam" % alignment_dir, sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.sorted.mkdup.bam.bai" % alignment_dir, sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.coverage.per-base.bed.gz" % alignment_dir, sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.%i.histo" % (kmer_dir, config["jellyfish_kmer_length"]), sample_id=config["sample_list"])

include: "workflow/rules/Preprocessing/Reference.smk"
include: "workflow/rules/QCFiltering/FastQC_raw.smk"
include: "workflow/rules/QCFiltering/Trimmomatic.smk"
include: "workflow/rules/QCFiltering/FastQC_filtered.smk"
#include: "workflow/rules/QCFiltering/Kmer.smk"  # not tested
include: "workflow/rules/Alignment/Alignment.smk"
include: "workflow/rules/Alignment/Coverage.smk"
#include: "workflow/rules/VariantCall/BSQR.smk"
#include: "workflow/rules/VariantCall/Genotyping.smk"

