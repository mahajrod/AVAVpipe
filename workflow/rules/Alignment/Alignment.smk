rule bwa_map:
    input:
        forward_reads=rules.trimmomatic.output.pe_forward,
        reverse_reads=rules.trimmomatic.output.pe_reverse,
        reference=config["reference"]
    output:
        bam=temp(alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.bam")
    params:
        fixmate_threads=config["fixmate_threads"],
        sort_threads=config["sort_threads"],
        markdup_threads=config["markdup_threads"],
        bwa_threads=config["bwa_threads"],
        per_thread_sort_mem="%sG" % config["per_thread_sort_mem"],
        tmp_prefix=alignment_dir_path / "{sample_id}/{sample_id}"
    log:
        bwa=log_dir_path / "{sample_id}/bwa.log",
        fixmate=log_dir_path / "{sample_id}/fixmate.log",
        sort=log_dir_path / "{sample_id}/sort.log",
        markdup=log_dir_path / "{sample_id}/markdup.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.map.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.map.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/alignment.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bwa_threads"] + config["sort_threads"] + config["fixmate_threads"] + config["markdup_threads"] + 1,
        time=config["map_time"],
        mem=config["per_thread_sort_mem"] * config["sort_threads"] * 1024 + config["bwa_mem_mb"] + config["fixmate_mem_mb"] + config["markdup_mem_mb"]
    threads: config["bwa_threads"] + config["sort_threads"] + config["fixmate_threads"] + config["markdup_threads"]
    shell:
        "bwa mem  -t {params.bwa_threads} {input.reference} <(gunzip -c {input.forward_reads}) <(gunzip -c {input.reverse_reads}) "
        "-R  \'@RG\\tID:{wildcards.sample_id}\\tPU:x\\tSM:{wildcards.sample_id}\\tPL:Illumina\\tLB:x\' 2>{log.bwa} | "
        "samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate}| "
        "samtools sort -T {params.tmp_prefix} -@ {params.sort_threads} -m {params.per_thread_sort_mem} 2>{log.sort}| "
        "samtools markdup -@ {params.markdup_threads} - {output.bam} 2>{log.markdup}"

rule index_bam:
    input:
        rules.bwa_map.output.bam
    output:
        temp(alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.bam.bai")
    log:
        std=log_dir_path / "{sample_id}/index.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.index.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.index.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/index.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["index_threads"],
        time=config["index_time"],
        mem=config["index_mem_mb"],
    threads: config["index_threads"]
    shell:
        "samtools index -@ {threads} {input} > {log.std} 2>&1"

