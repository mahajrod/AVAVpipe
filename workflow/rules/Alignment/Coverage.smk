rule mosdepth:
    input:
        bam=rules.bwa_map.output.bam,
        bai=rules.index_bam.output
    output:
        "%s/{sample_id}/{sample_id}.coverage.per-base.bed.gz" % alignment_dir
    params:
        min_mapping_quality=config["mosdepth_min_mapping_quality"],
        output_pefix="%s/{sample_id}/{sample_id}.coverage" % alignment_dir
    log:
        std="%s/{sample_id}/mosdepth.log" % log_dir,
        cluster_log="%s/{sample_id}.mosdepth.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.mosdepth.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/mosdepth.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["mosdepth_threads"],
        time=config["mosdepth_time"],
        mem=config["mosdepth_mem_mb"],
    threads: config["mosdepth_threads"]
    shell:
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {params.output_pefix} {input.bam} > {log.std} 2>&1"


