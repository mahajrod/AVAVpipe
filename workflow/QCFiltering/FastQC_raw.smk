rule fastqc_raw:
    input:
        "%s/{sample_id}/{sample_id}_1.fastq.gz" % (config["sample_dir"]),
        "%s/{sample_id}/{sample_id}_2.fastq.gz" % (config["sample_dir"])
    output:
        directory("%s/{sample_id}/raw/" % fastqc_dir)
    params:
        kmer=10
    log:
        std="%s/{sample_id}/fastqc_raw.log" % log_dir,
        cluster_log="%s/{sample_id}.fastqc_raw.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.fastqc_raw.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/fastqc_raw.benchmark.txt" % benchmark_dir
    #conda:
    #    "../%s" % config["conda_config"]
    resources:
        cpus=config["fastqc_threads"],
        time=config["fastqc_time"],
        mem=config["fastqc_mem_mb"],

    threads:
        config["fastqc_threads"]
    shell:
        "mkdir -p {output}; fastqc --nogroup -k {params.kmer} -t {threads} -o {output} {input} 1>{log.std} 2>&1"
