localrules: prepare_regions

rule ref_faidx:
    input:
        config["reference"]
    output:
        "%s.fai" % config["reference"]
    log:
        std="%s/faidx.log" % log_dir,
        cluster_log="%s/faidx.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/faidx.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/faidx.benchmark.txt" % benchmark_dir
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
        reference_dict_path #"%s.dict" % (os.path.splitext(config["reference"])[0])
    log:
        std="%s/dict.log" % log_dir,
        cluster_log="%s/dict.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/dict.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/dict.benchmark.txt" % benchmark_dir
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

rule prepare_recalibration_regions:
    input:
         fai=rules.ref_faidx.output,
         blacklist=reference_blacklist_path, #"%s.blacklist" % (os.path.splitext(config["reference"])[0]),
         whitelist=reference_whitelist_path, #"%s.whitelist" % (os.path.splitext(config["reference"])[0])
    output:
         reference_region_correspondence_path #"%s/regions/scaffold_to_region.correspondence" % os.path.dirname(config["reference"])
    params:
        max_region_length=config["split_regions_max_region_length"],
        max_seq_number=config["split_regions_max_seq_number"],
        region_file_format=config["split_regions_region_file_format"],
        min_scaffold_length=config["split_regions_min_scaffold_length"],
        #scaffold_whitelist=config["split_regions_scaffold_whitelist"],
        #scaffold_blacklist=config["split_regions_scaffold_blacklist"],
        output_dir=reference_region_dir_path #"%s/regions/" % os.path.dirname(config["reference"])
    log:
        std="%s/prepare_regions.log" % log_dir,
        cluster_log="%s/prepare_regions.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/prepare_regions.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/prepare_regions.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["prepare_regions_threads"],
        time=config["prepare_regions_time"],
        mem=config["prepare_regions_mem_mb"],
    threads:
        config["prepare_regions_threads"]
    shell:
         "workflow/scripts/split_regions.py -s -f {input.fai} -m {params.max_region_length}"
         " -w {input.whitelist} -b {input.blacklist} "
         " -n {params.max_seq_number} -g {params.region_file_format} -x {params.min_scaffold_length} "
         " -o {params.output_dir} > {log.std} 2>&1"