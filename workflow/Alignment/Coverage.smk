rule mosdepth:
    input:
        rules.bwa_map.output
    output:
        "%s/{sample_id}/{sample_id}.coverage.per-base.bed.gz" % alignment_dir
    params:
        min_mapping_quality=config["mosdpth_min_mapping_quality"],
        output_pefix="%s/{sample_id}/{sample_id}.coverage" % alignment_dir
    log:
        "%s/{sample_id}/mosdepth.log" % log_dir
    benchmark:
        "%s/{sample_id}/mosdepth.benchmark.txt" % benchmark_dir
    threads: config["mosdepth_threads"]
    shell:
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {params.output_pefix} {input}"


