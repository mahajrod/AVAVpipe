rule baserecalibrator:
    input:
        region="%s/intervals/region_{region_id}.list" % reference_region_dir_path,
        bam="%s/{sample_id}/{sample_id}.sorted.mkdup.bam" % alignment_dir,
        bai="%s/{sample_id}/{sample_id}.sorted.mkdup.bam.bai" % alignment_dir,
        reference=config["reference"],
        known_variants_vcf_list=known_variants_vcf_list
    output:
        table="%s/{sample_id}/baserecalibrator/{sample_id}.region_{region_id}.sorted.recal.table" % alignment_dir,
    log:
        std="%s/{sample_id}/baserecalibrator/baserecalibrator.region_{region_id}.log" % log_dir,
        cluster_log="%s/baserecalibrator/{sample_id}.baserecalibrator.region_{region_id}.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/baserecalibrator/{sample_id}.baserecalibrator.region_{region_id}.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/baserecalibrator/baserecalibrator.region_{region_id}.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["baserecalibrator_threads"],
        time=config["baserecalibrator_time"],
        mem=config["baserecalibrator_mem_mb"],
    threads: config["baserecalibrator_threads"]
    shell:
        "gatk --java-options '-Xmx{resources.mem}m' BaseRecalibrator -R {input.reference} -L {input.region} "
        "-I {input.bam} -O {output.table} --known-sites %s > {log.std} 2>&1" % " --known-sites ".join(list(map(lambda s: str(s), known_variants_vcf_list)))

rule gatherbsqrreports:
    input:
        expand("%s/{sample_id}/baserecalibrator/{sample_id}.region_{region_id}.sorted.recal.table" % alignment_dir,
               region_id=glob_wildcards("%s/intervals/region_{region_id}.list" % reference_region_dir_path)[0],
               sample_id=config["sample_list"])
    output:
        "%s/{sample_id}/{sample_id}.sorted.recal.table" % alignment_dir
    params:
        ids=glob_wildcards("%s/intervals/region_{region_id}.list" % reference_region_dir_path),
        recal_table_list="%s/{sample_id}/{sample_id}.sorted.recal.table.list" % alignment_dir
    log:
        std="%s/{sample_id}.gatherbsqrreports.log" % log_dir,
        cluster_log="%s/{sample_id}.gatherbsqrreports.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.gatherbsqrreports.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/gatherbsqrreports.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["gatherbsqrreports_threads"],
        time=config["gatherbsqrreports_time"],
        mem=config["gatherbsqrreports_mem_mb"],
    threads: config["gatherbsqrreports_threads"]
    shell:
        "echo {input} | sed 's/ /\n/' > {params.recal_table_list}; "
        "gatk --java-options '-Xmx{resources.mem}m' GatherBQSRReports -I {params.recal_table_list} -O {output}"

"""
rule applybsqr:
    input:
        bam=rules.bwa_map.output.bam,
        reference=config["reference"],
        table=rules.baserecalibrator.output
    output:
        "%s/{sample_id}/{sample_id}.sorted.recalibrated.bam" %  alignment_dir
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
"""