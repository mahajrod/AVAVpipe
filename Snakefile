from pathlib import Path

out_dir_path = Path(config["out_dir"])
"""
fastqc_dir = out_dir_path / config["fastqc_dir"]
kmer_dir = out_dir_path / config["kmer_dir"]
filtered_read_dir = out_dir_path / config["filtered_read_dir"]
alignment_dir = out_dir_path / config["alignment_dir"]
snpcall_dir = out_dir_path / config["snpcall_dir"]
joint_snpcall_dir = out_dir_path / config["joint_snpcall_dir"]
log_dir = out_dir_path / config["log_dir"]
error_dir  = out_dir_path / config["error_dir"]
benchmark_dir =  out_dir_path / config["benchmark_dir"]
"""

fastqc_dir_path = out_dir_path / config["fastqc_dir"]
kmer_dir_path = out_dir_path / config["kmer_dir"]
filtered_read_dir_path = out_dir_path / config["filtered_read_dir"]
alignment_dir_path = out_dir_path / config["alignment_dir"]
snpcall_dir_path = out_dir_path / config["snpcall_dir"]
joint_snpcall_dir_path = out_dir_path / config["joint_snpcall_dir"]
joint_snpcall_per_sample_dir_path = joint_snpcall_dir_path / "per_sample"
log_dir_path = out_dir_path / config["log_dir"]
error_dir_path  = out_dir_path / config["error_dir"]
benchmark_dir_path =  out_dir_path / config["benchmark_dir"]
cluster_log_dir_path = Path(config["cluster_log_dir"])

known_variants_dir_path = Path(config["known_variants_dir"])
known_variants_vcf_list = list(known_variants_dir_path.glob("*.vcf")) + list(known_variants_dir_path.glob("*.vcf.gz"))
known_variants_mills_path = known_variants_dir_path / config["mills_vcf"] #"Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_variants_axiompoly_path = known_variants_dir_path / config["axiompoly_vcf"] #"Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
known_variants_dbsnp_path = known_variants_dir_path / config["dbsnp_vcf"] #"Homo_sapiens_assembly38.dbsnp138.vcf"
known_variants_hapmap_path = known_variants_dir_path / config["hapmap_vcf"] #"hapmap_3.3.hg38.vcf.gz"
known_variants_omni_path = known_variants_dir_path / config["omni_vcf"] #"1000G_omni2.5.hg38.vcf.gz"
known_variants_1000g_path = known_variants_dir_path / config["g1000_vcf"] #"1000G_phase1.snps.high_confidence.hg38.vcf.gz"

sample_dir_path = Path(config["sample_dir"])
reference_path = Path(config["reference"])
reference_dir_path = reference_path.parent
reference_fai_path = Path(str(reference_path) + ".fai")
reference_dict_path = reference_dir_path / (reference_path.stem + ".dict")
reference_blacklist_path = reference_dir_path / (reference_path.stem + ".blacklist")
reference_whitelist_path = reference_dir_path / (reference_path.stem + ".whitelist")

reference_genotyping_whitelist_path = reference_dir_path / (reference_path.stem + ".genotyping.whitelist")
reference_genotyping_whitelist_intervals_path = reference_dir_path / (reference_path.stem + ".genotyping.whitelist.intervals")

reference_recalibration_region_dir_path = reference_dir_path / "recalibration_regions"
reference_recalibration_region_correspondence_path = reference_recalibration_region_dir_path / "scaffold_to_region.correspondence"

# if "sample_list" key is absent in config variable, use folder names from config["sample_dir"] as sample ids
if "sample_list" not in config:
    config["sample_list"] = [d.name for d in sample_dir_path.iterdir() if d.is_dir()]

localrules: all

rule all:
    input:
        reference_fai_path,
        reference_dict_path,
        reference_recalibration_region_correspondence_path,
        reference_genotyping_whitelist_intervals_path,
        expand(log_dir_path / "{sample_id}/fastqc_raw.log", sample_id=config["sample_list"]),
        expand(log_dir_path / "{sample_id}/fastqc_filtered.log", sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.sorted.mkdup.bam" % alignment_dir, sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.sorted.mkdup.bam.bai" % alignment_dir, sample_id=config["sample_list"]),
        expand(alignment_dir_path / "{sample_id}/{sample_id}.coverage.per-base.bed.gz", sample_id=config["sample_list"]),
        expand(kmer_dir_path / ("{sample_id}/{sample_id}.%i.histo" % config["jellyfish_kmer_length"]), sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.sorted.recal.table" % alignment_dir, sample_id=config["sample_list"]),
        #expand("%s/{sample_id}/{sample_id}.sorted.mkdup.recalibrated.bam" %  alignment_dir, sample_id=config["sample_list"] ),
        expand(snpcall_dir_path / "{sample_id}/{sample_id}.gvcf", sample_id=config["sample_list"] ),
        directory(joint_snpcall_dir_path / "gvcf_database/callset.json"),
        joint_snpcall_dir_path / "all_samples.recalibrated.good.vcf.gz",
        expand(joint_snpcall_per_sample_dir_path/ "{sample_id}.recalibrated.good.vcf.gz", sample_id=config["sample_list"]),
        joint_snpcall_dir_path / "all_samples.recalibrated.good.variant_calling_summary_metrics"
        #joint_snpcall_dir_path / "all_samples.recalibrated.good.varianteval",


include: "workflow/rules/Preprocessing/Reference.smk"
include: "workflow/rules/QCFiltering/FastQC_raw.smk"
include: "workflow/rules/QCFiltering/Trimmomatic.smk"
include: "workflow/rules/QCFiltering/FastQC_filtered.smk"
include: "workflow/rules/QCFiltering/Kmer.smk"
include: "workflow/rules/Alignment/Alignment.smk"
include: "workflow/rules/Alignment/Coverage.smk"
include: "workflow/rules/VariantCall/BQSR.smk"
include: "workflow/rules/VariantCall/Genotyping.smk"
include: "workflow/rules/VariantCall/VQSR.smk"
include: "workflow/rules/Evaluation/Evaluation.smk"
