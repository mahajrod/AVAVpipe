rule fastqc_raw:
    input:
        "%s/{sample_id}/{sample_id}_1.fastq.gz" % (config["sample_dir"]),
        "%s/{sample_id}/{sample_id}_2.fastq.gz" % (config["sample_dir"])
    output:
        directory("%s/{sample_id}/raw/" % fastqc_dir)
    params:
        kmer=10
    log:
        "%s/{sample_id}/fastqc_raw.log" % log_dir
    benchmark:
        "%s/{sample_id}/fastqc_raw.benchmark.txt" % benchmark_dir
    conda:
        config["conda_config"]
    resources:
        cpus=config["fastqc_threads"],
        time=config["fastqc_time"],
        mem=config["fastqc_mem_mb"],
        slurm_log="%s/{sample_id}/fastqc_raw.slurm.log" % log_dir,
        slurm_err="%s/{sample_id}/fastqc_raw.slurm.err" % log_dir
    threads:
        config["fastqc_threads"]
    shell:
        "mkdir -p {output}; fastqc --nogroup -k {params.kmer} -t {threads} -o {output} {input} 1>{log} 2>&1"
