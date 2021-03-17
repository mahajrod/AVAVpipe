rule baserecalibrator:
    input:
        bam=rules.bwa_map.output.bam,
        reference=config["reference"],
        known_sites_vcf_list=expand("--known-sites {vcf}", vcf=config["known_sites_vcf_list"])
    output:
        "%s/{sample_id}/{sample_id}.sorted.recal.table" % alignment_dir
    log:
        std="%s/{sample_id}/baserecalibrator.log" % log_dir,
        cluster_log="%s/{sample_id}.baserecalibrator.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.baserecalibrator.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/baserecalibrator.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["baserecalibrator_threads"],
        time=config["baserecalibrator_time"],
        mem=config["baserecalibrator_mem_mb"],
    threads: config["baserecalibrator_threads"]
    shell:
        "gatk BaseRecalibrator -R {input.reference} -I {input.bam} {input.known_sites_vcf_list} -O {output} > {log.std} 2>&1"

rule applybsqr:
    input:
        bam=rules.bwa_map.output.bam,
        reference=config["reference"],
        table=rules.baserecalibrator.output
    output:
        "%s/{sample_id}/{sample_id}.sorted.recalibrated.bam" % alignment_dir
    log:
        std="%s/{sample_id}/applybsqr.log" % log_dir,
        cluster_log="%s/{sample_id}.applybsqr.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.applybsqr.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/applybsqr.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["applybsqr_threads"],
        time=config["applybsqr_time"],
        mem=config["applybsqr_mem_mb"]
    threads: config["applybsqr_threads"]
    shell:
        "gatk ApplyBQSR -R {input.reference} -I {input.bam} --bqsr-recal-file {input.table} -O {output} > {log.std} 2>&1"