rule bwa_map:
    input:
        forward_reads=rules.trimmomatic.output.pe_forward, #"%s/{sample_id}/{sample_id}.trimmed_1.fastq.gz" % filtered_read_dir,
        reverse_reads=rules.trimmomatic.output.pe_reverse,#"%s/{sample_id}/{sample_id}.trimmed_2.fastq.gz" % filtered_read_dir
        reference=config["reference"]
    output:
        bam="%s/{sample_id}/{sample_id}.sorted.mkdup.bam" % alignment_dir,
        #tmp_prefix="%s/{sample_id}/{sample_id}" % alignment_dir
    params:
        fixmate_threads=config["fixmate_threads"],
        sort_threads=config["sort_threads"],
        markdup_threads=config["markdup_threads"],
        bwa_threads=config["bwa_threads"],
        per_thread_sort_mem="%sG" % config["per_thread_sort_mem"]
    log:
        bwa="%s/{sample_id}/bwa.log" % log_dir,
        fixmate="%s/{sample_id}/fixmate.log" % log_dir,
        sort="%s/{sample_id}/sort.log" % log_dir,
        markdup="%s/{sample_id}/markdup.log" % log_dir,
        cluster_log="%s/{sample_id}.map.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.map.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/alignment.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["bwa_threads"] + config["sort_threads"] + config["fixmate_threads"] + config["markdup_threads"] + 1,
        time=config["map_time"],
        mem=config["per_thread_sort_mem"] * config["sort_threads"] * 1024 + config["bwa_mem_mb"] + config["fixmate_mem_mb"] + config["markdup_mem_mb"]
    threads: config["bwa_threads"] + config["sort_threads"] + config["fixmate_threads"] + config["markdup_threads"]
    shell:
        "bwa mem  -t {params.bwa_threads} {input.reference} <(gunzip -c {input.forward_reads}) <(gunzip -c {input.reverse_reads}) "
        "-R  \'@RG\\tID:{wildcards.sample_id}\\tPU:x\\tSM:{wildcards.sample_id}\\tPL:Illumina\\tLB:x\' 2>{log.bwa} | "
        "samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate}| "
        "samtools sort -T {output.tmp_prefix} -@ {params.sort_threads} -m {params.per_thread_sort_mem} 2>{log.sort}| "
        "samtools markdup -@ {params.markdup_threads} - {output.bam} 2>{log.markdup}"

rule index_bam:
    input:
        rules.bwa_map.output.bam
    output:
        "%s/{sample_id}/{sample_id}.sorted.mkdup.bam.bai" % alignment_dir
    log:
        std="%s/{sample_id}/index.log" % log_dir,
        cluster_log="%s/{sample_id}.index.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.index.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/index.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["index_threads"],
        time=config["index_time"],
        mem=config["index_mem_mb"],
    threads: config["index_threads"]
    shell:
        "samtools index -@ {threads} {input} > {log.std} 2>&1"

