rule trimmomatic:
    input:
        "%s/{sample_id}/{sample_id}_1.fastq.gz" % (config["sample_dir"]),
        "%s/{sample_id}/{sample_id}_2.fastq.gz" % (config["sample_dir"])
    output:
        pe_forward="%s/{sample_id}/{sample_id}.trimmed_1.fastq.gz" % filtered_read_dir,
        se_forward="%s/{sample_id}/{sample_id}.trimmed_1.se.fastq.gz" % filtered_read_dir,
        pe_reverse="%s/{sample_id}/{sample_id}.trimmed_2.fastq.gz" % filtered_read_dir,
        se_reverse="%s/{sample_id}/{sample_id}.trimmed_2.se.fastq.gz" % filtered_read_dir,
    params:
        adapters=config["trimmomatic_adapters"],
        illumina_clip="2:30:10:1",
        window_size=8,
        window_quality=20,
        minlength=50
    log:
        std="%s/{sample_id}/trimmomatic.log" % log_dir,
        slurm_log="{sample_id}_trimmomatic.slurm.log",
        slurm_err="{sample_id}_trimmomatic.slurm.err"
        #slurm_log="%s/{sample_id}/trimmomatic.slurm.log" % log_dir,
        #slurm_err="%s/{sample_id}/trimmomatic.slurm.err" % log_dir
    benchmark:
        "%s/{sample_id}/trimmomatic.benchmark.txt" % benchmark_dir
    conda:
        config["conda_config"]
    resources:
        cpus=config["trimmomatic_threads"],
        time=config["trimmomatic_time"],
        mem=config["trimmomatic_mem_mb"],
    threads:
        config["trimmomatic_threads"]
    shell:
         "trimmomatic PE -threads {threads} -phred33 {input} {output} "
         "ILLUMINACLIP:{params.adapters}:{params.illumina_clip} "
         "SLIDINGWINDOW:{params.window_size}:{params.window_quality} "
         "MINLEN:{params.minlength} > {log.std} 2>&1"