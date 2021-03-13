rule baserecalibrator:
    input:
        bam=rules.bwa_map.output,
        reference=config["reference"],
        known_sites_vcf_list=expand("--known-sites {vcf}", vcf=config["known_sites_vcf_list"])
    output:
        "%s/{sample_id}/{sample_id}.sorted.recal.table" % alignment_dir
    log:
        std="%s/{sample_id}/baserecalibrator.log" % log_dir,
        slurm_log="%s/{sample_id}/baserecalibrator.slurm.log" % log_dir,
        slurm_err="%s/{sample_id}/baserecalibrator.slurm.err" % log_dir
    benchmark:
        "%s/{sample_id}/baserecalibrator.benchmark.txt" % benchmark_dir
    resources:
        cpus=config["baserecalibrator_threads"],
        time=config["baserecalibrator_time"],
        mem=config["baserecalibrator_mem_mb"],

    threads: config["baserecalibrator_threads"]
    shell:
        "gatk BaseRecalibrator -R {input.reference} -I {input.bam} {input.known_sites_vcf_list} -O {output} > {log.std} 2>&1"

rule applybsqr:
    input:
        bam=rules.bwa_map.output,
        reference=config["reference"],
        table=rules.baserecalibrator.output
    output:
        "%s/{sample_id}/{sample_id}.sorted.recalibrated.bam" % alignment_dir
    log:
        std="%s/{sample_id}/applybsqr.log" % log_dir,
        slurm_log="%s/{sample_id}/applybsqr.slurm.log" % log_dir,
        slurm_err="%s/{sample_id}/applybsqr.slurm.err" % log_dir
    benchmark:
        "%s/{sample_id}/applybsqr.benchmark.txt" % benchmark_dir
    resources:
        cpus=config["applybsqr_threads"],
        time=config["applybsqr_time"],
        mem=config["applybsqr_mem_mb"]
    threads: config["applybsqr_threads"]
    shell:
        "gatk ApplyBQSR -R {input.reference} -I {input.bam} --bqsr-recal-file {input.table} -O {output} > {log.std} 2>&1"