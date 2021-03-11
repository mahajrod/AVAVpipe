rule fastqc_filtered:
    input:
        rules.trimmomatic.output
    output:
        directory("%s/{sample_id}/filtered/" % fastqc_dir)
    params:
        kmer=10
    log:
        "%s/{sample_id}/fastqc_filtered.log" % log_dir
    benchmark:
        "%s/{sample_id}/fastqc_filtered.benchmark.txt" % benchmark_dir
    conda:
        config["conda_config"]
    threads:
        config["fastqc_threads"]
    shell:
        "mkdir -p {output}; fastqc --nogroup -k {params.kmer} -t {threads} -o {output} {input} 1>{log} 2>&1"
