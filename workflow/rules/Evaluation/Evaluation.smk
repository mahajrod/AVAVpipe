rule variantcalling_metrics:
    input:
        vcf=rules.select_good_variants.output.vcf,
        dbsnp=known_variants_dbsnp_path,
        reference_dict=reference_dict_path
    output:
        summary=joint_snpcall_dir_path / "all_samples.recalibrated.good.variant_calling_summary_metrics",
        detailed=joint_snpcall_dir_path / "all_samples.recalibrated.good.variant_calling_detail_metrics"
    params:
        output_prefix=joint_snpcall_dir_path / "all_samples.recalibrated.good"
    log:
        std=log_dir_path / "variantcalling_metrics.log",
        cluster_log=cluster_log_dir_path / "variantcalling_metrics.cluster.log",
        cluster_err=cluster_log_dir_path / "variantcalling_metrics_.cluster.err"
    benchmark:
        benchmark_dir_path / "variantcalling_metrics.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["variantcalling_metrics_threads"],
        time=config["variantcalling_metrics_time"],
        mem=config["variantcalling_metrics_mem_mb"],
    threads: config["variantcalling_metrics_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' CollectVariantCallingMetrics"
        " -I {input.vcf} --DBSNP {input.dbsnp} -SD {input.reference_dict}"
        " -O {params.output_prefix} > {log.std} 2>&1"

"""
rule varianteval:
    input:
        vcf=rules.select_good_variants.output.vcf,
        dbsnp=known_variants_dbsnp_path,
        reference_=reference_path
    output:
        joint_snpcall_dir_path / "all_samples.recalibrated.varianteval.good"
    #params:
    #    metrics=" -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary "
    log:
        std=log_dir_path / "varianteval_metrics.log",
        cluster_log=cluster_log_dir_path / "varianteval_metrics.cluster.log",
        cluster_err=cluster_log_dir_path / "varianteval_metrics.cluster.err"
    benchmark:
        benchmark_dir_path / "varianteval_metrics.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["varianteval_metrics_threads"],
        time=config["varianteval_metrics_time"],
        mem=config["varianteval_metrics_mem_mb"],
    threads: config["varianteval_metrics_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' VariantEval"
        " -R {input.reference_} --eval {input.vcf} -D {input.dbsnp} " #  {params.metrics}"
        " -O {output} > {log.std} 2>&1"
"""