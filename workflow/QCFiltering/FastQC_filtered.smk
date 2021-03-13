rule fastqc_filtered:
    input:
        rules.trimmomatic.output
    output:
        directory("%s/{sample_id}/filtered/" % fastqc_dir)
    params:
        kmer=10
    log:
        std="%s/{sample_id}/fastqc_filtered.log" % log_dir,
        cluster_log="%s/{sample_id}.fastqc_filtered.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.fastqc_filtered.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/fastqc_filtered.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["fastqc_threads"],
        time=config["fastqc_time"],
        mem=config["fastqc_mem_mb"],
    threads:
        config["fastqc_threads"]
    shell:
        "mkdir -p {output}; fastqc --nogroup -k {params.kmer} -t {threads} -o {output} {input} 1>{log.std} 2>&1"
