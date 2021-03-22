rule variantcalling_metrics:
    input:
        vcf=rules.applyvsqr_snp.output,
        dbsnp=known_variants_dbsnp_path,
        reference_dict=reference_dict_path
    output:
        log_dir /"variantcalling_metrics.log"
    params:
        output_prefix=joint_snpcall_dir / "all_samples.recalibrated.metrics"
    log:
        std="%s/variantcalling_metrics.log" % log_dir,
        cluster_log="%s/variantcalling_metrics.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/variantcalling_metrics_.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/variantcalling_metrics.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["variantcalling_metrics_threads"],
        time=config["variantcalling_metrics_time"],
        mem=config["variantcalling_metrics_mem_mb"],
    threads: config["variantcalling_metrics_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' CollectVariantCallingMetrics"
        " -I {input.vcf} --DBSNP {input.dbsnp} -SD {input.reference_dict}"
        " -O {params.output_prefix} > {log.std} 2>&1"

rule varianteval:
    input:
        vcf=rules.applyvsqr_snp.output,
        dbsnp=known_variants_dbsnp_path,
        reference_=reference_path
    output:
        joint_snpcall_dir / "all_samples.recalibrated.varianteval"
    params:
        metrics=" -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary "
    log:
        std="%s/varianteval_metrics.log" % log_dir,
        cluster_log="%s/varianteval_metrics.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/varianteval_metrics.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/varianteval_metrics.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["varianteval_metrics_threads"],
        time=config["varianteval_metrics_time"],
        mem=config["varianteval_metrics_mem_mb"],
    threads: config["varianteval_metrics_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' VariantEval"
        " -R {input.reference_} --eval {input.vcf} -D {input.dbsnp}  {params.metrics}"
        " -O {output} > {log.std} 2>&1"
