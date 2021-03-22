localrules: prepare_recalibration_regions, prepare_genotyping_whitelist_intervals

rule ref_faidx:
    input:
        config["reference"]
    output:
        "%s.fai" % config["reference"]
    log:
        std=log_dir_path / "faidx.log",
        cluster_log=cluster_log_dir_path / "faidx.cluster.log",
        cluster_err=cluster_log_dir_path / "faidx.cluster.err"
    benchmark:
        benchmark_dir_path / "faidx.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["faidx_threads"],
        time=config["faidx_time"],
        mem=config["faidx_mem_mb"],
    threads:
        config["faidx_threads"]
    shell:
         "samtools faidx {input} > {log.std} 2>&1"

rule ref_dict:
    input:
        config["reference"]
    output:
        reference_dict_path
    log:
        std=log_dir_path / "/dict.log",
        cluster_log=cluster_log_dir_path / "dict.cluster.log",
        cluster_err=cluster_log_dir_path / "dict.cluster.err"
    benchmark:
        benchmark_dir_path / "dict.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["dict_threads"],
        time=config["dict_time"],
        mem=config["dict_mem_mb"],
    threads:
        config["dict_threads"]
    shell:
         "picard CreateSequenceDictionary -R {input} > {log.std} 2>&1"

rule prepare_genotyping_whitelist_intervals:
    input:
        fai=rules.ref_faidx.output,
        whitelist=reference_genotyping_whitelist_path
    output:
        reference_genotyping_whitelist_intervals_path
    log:
        std=log_dir_path / "prepare_genotyping_whitelist_intervals.log",
        cluster_log=cluster_log_dir_path / "prepare_genotyping_whitelist_intervals.cluster.log",
        cluster_err=cluster_log_dir_path / "prepare_genotyping_whitelist_intervals.cluster.err"
    benchmark:
        benchmark_dir_path / "prepare_genotyping_whitelist_intervals.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["prepare_genotyping_whitelist_intervals_threads"],
        time=config["prepare_genotyping_whitelist_intervals_time"],
        mem=config["prepare_genotyping_whitelist_intervals_mem_mb"],
    threads:
        config["prepare_genotyping_whitelist_intervals_threads"]
    shell:
         "workflow/scripts/fai2intervals.py -f {input.fai} -w {input.whitelist} -o {output} > {log.std} 2>&1"

rule prepare_recalibration_regions:
    input:
         fai=rules.ref_faidx.output,
         blacklist=reference_blacklist_path,
         whitelist=reference_whitelist_path
    output:
         reference_region_correspondence_path
    params:
        max_region_length=config["split_regions_max_region_length"],
        max_seq_number=config["split_regions_max_seq_number"],
        region_file_format=config["split_regions_region_file_format"],
        min_scaffold_length=config["split_regions_min_scaffold_length"],
        output_dir=reference_region_dir_path
    log:
        std=log_dir_path / "prepare_recalibration_regions.log",
        cluster_log=cluster_log_dir_path / "prepare_recalibration_regions.cluster.log",
        cluster_err=cluster_log_dir_path / "prepare_recalibration_regions.cluster.err"
    benchmark:
        benchmark_dir_path / "prepare_recalibration_regions.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["prepare_recalibration_regions_threads"],
        time=config["prepare_recalibration_regions_time"],
        mem=config["prepare_recalibration_regions_mem_mb"],
    threads:
        config["prepare_recalibration_regions_threads"]
    shell:
         "workflow/scripts/split_regions.py -s -f {input.fai} -m {params.max_region_length}"
         " -w {input.whitelist} -b {input.blacklist} "
         " -n {params.max_seq_number} -g {params.region_file_format} -x {params.min_scaffold_length} "
         " -o {params.output_dir} > {log.std} 2>&1"