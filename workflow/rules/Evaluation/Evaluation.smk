rule variantcalling_metrics:
    input:
        vcf=joint_snpcall_dir / "all_samples.{variant_type}.recalibrated.vcf.gz",
        dbsnp=known_variants_dbsnp_path,
        reference_dict=reference_dict_path
    output:
        "%s/variantcalling_metrics_{variant_type}.log" % log_dir
    params:
        output_prefix=joint_snpcall_dir / "all_samples.{variant_type}.recalibrated.metrics"
    log:
        std="%s/variantcalling_metrics_{variant_type}.log" % log_dir,
        cluster_log="%s/variantcalling_metrics_{variant_type}cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/variantcalling_metrics_{variant_type}.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/variantcalling_metrics_{variant_type}.benchmark.txt" % benchmark_dir
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
        vcf=joint_snpcall_dir / "all_samples.{variant_type}.recalibrated.vcf.gz",
        dbsnp=known_variants_dbsnp_path,
        reference_=reference_path
    output:
        joint_snpcall_dir / "all_samples.{variant_type}.recalibrated.varianteval"
    params:
        metrics=" -noEV -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator -EV CountVariants -EV MultiallelicSummary "
    log:
        std="%s/varianteval_metrics_{variant_type}.log" % log_dir,
        cluster_log="%s/varianteval_metrics_{variant_type}cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/varianteval_metrics_{variant_type}.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/varianteval_metrics_{variant_type}.benchmark.txt" % benchmark_dir
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
