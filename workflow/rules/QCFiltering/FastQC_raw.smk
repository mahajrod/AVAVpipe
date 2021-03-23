rule fastqc_raw:
    input:
        sample_dir_path / "{sample_id}/{sample_id}_1.fastq.gz",
        sample_dir_path / "{sample_id}/{sample_id}_2.fastq.gz"
    output:
        directory(fastqc_dir_path / "{sample_id}/raw/")
    params:
        kmer=7
    log:
        std=log_dir_path / "{sample_id}/fastqc_raw.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.fastqc_raw.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.fastqc_raw.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/fastqc_raw.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["fastqc_threads"],
        time=config["fastqc_time"],
        mem=config["fastqc_mem_mb"],

    threads:
        config["fastqc_threads"]
    shell:
        "mkdir -p {output}; fastqc --nogroup -k {params.kmer} -t {threads} -o {output} {input} 1>{log.std} 2>&1"
