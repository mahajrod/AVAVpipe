rule baserecalibrator:
    input:
        region=reference_recalibration_region_dir_path / "intervals/region_{region_id}.list",
        bam=rules.bwa_map.output.bam,
        bai=rules.index_bam.output,
        reference=config["reference"],
        known_variants_vcf_list=known_variants_vcf_list
    output:
        table=alignment_dir_path / "{sample_id}/baserecalibrator/{sample_id}.region_{region_id}.sorted.mkdup.recal.table"
    log:
        std=log_dir_path / "{sample_id}/baserecalibrator/baserecalibrator.region_{region_id}.log",
        cluster_log=cluster_log_dir_path / "baserecalibrator/{sample_id}.baserecalibrator.region_{region_id}.cluster.log",
        cluster_err=cluster_log_dir_path / "baserecalibrator/{sample_id}.baserecalibrator.region_{region_id}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/baserecalibrator/baserecalibrator.region_{region_id}.benchmark.txt"
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
        lambda wildcards:  expand(alignment_dir_path / "{sample_id}/baserecalibrator/{sample_id}.region_{region_id}.sorted.mkdup.recal.table",
                                  region_id=glob_wildcards(reference_recalibration_region_dir_path / "/intervals/region_{region_id}.list")[0],
                                  sample_id=[wildcards.sample_id])
    output:
        alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.recal.table"
    params:
        input_files=alignment_dir_path / "{sample_id}/baserecalibrator/*",
        recal_table_list=alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.recal.table.list"
    log:
        std=log_dir_path / "{sample_id}.gatherbsqrreports.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.gatherbsqrreports.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.gatherbsqrreports.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/gatherbsqrreports.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["gatherbsqrreports_threads"],
        time=config["gatherbsqrreports_time"],
        mem=config["gatherbsqrreports_mem_mb"],
    threads: config["gatherbsqrreports_threads"]
    shell:
        "ls {params.input_files} | sort -V > {params.recal_table_list}; "
        "gatk --java-options '-Xmx{resources.mem}m' GatherBQSRReports -I {params.recal_table_list} -O {output} > {log.std} 2>&1"

rule applybsqr:
    input:
        bam=rules.bwa_map.output.bam,
        reference=config["reference"],
        table=rules.gatherbsqrreports.output
    output:
        bam=alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.recalibrated.bam",
        bai=alignment_dir_path / "{sample_id}/{sample_id}.sorted.mkdup.recalibrated.bai"
    log:
        std=log_dir_path / "{sample_id}/applybsqr.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.applybsqr.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.applybsqr.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/applybsqr.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["applybsqr_threads"],
        time=config["applybsqr_time"],
        mem=config["applybsqr_mem_mb"]
    threads: config["applybsqr_threads"]
    shell:
        "gatk ApplyBQSR -R {input.reference} -I {input.bam} --bqsr-recal-file {input.table} -O {output} > {log.std} 2>&1"
