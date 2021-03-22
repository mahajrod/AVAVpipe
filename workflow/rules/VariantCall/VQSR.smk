
rule excess_filter:
    input:
        vcf=rules.genotypegvcfs.output
    output:
        joint_snpcall_dir / "all_samples.excesshet.vcf.gz"
    params:
        filter="'ExcessHet > 54.69'",
        filter_name="ExcessHet"
    log:
        std="%s/excess_filter.log" % log_dir,
        cluster_log="%s/excess_filter.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/excess_filter.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/excess_filter.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["excess_filter_threads"],
        time=config["excess_filter_time"],
        mem=config["excess_filter_mem_mb"],
    threads: config["excess_filter_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' VariantFiltration -V {input.vcf} "
        " --filter-expression {params.filter} --filter-name {params.filter_name}"
        " -O {output}> {log.std} 2>&1"

rule extract_sites:
    input:
        vcf=rules.excess_filter.output
    output:
        joint_snpcall_dir / "all_samples.excesshet.sites_only.vcf.gz"
    log:
        std="%s/extract_sites.log" % log_dir,
        cluster_log="%s/extract_sites.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/extract_sites.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/extract_sites.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["extract_sites_threads"],
        time=config["extract_sites_time"],
        mem=config["extract_sites_mem_mb"],
    threads: config["extract_sites_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' MakeSitesOnlyVcf -I {input.vcf} "
        " -O {output}> {log.std} 2>&1"

rule variantrecalibrator_indel:
    input:
        vcf=rules.extract_sites.output
    output:
        tranches=joint_snpcall_dir / "all_samples.indel.tranches",
        recal_table=joint_snpcall_dir / "all_samples.indel.recal.table",
    params:
        mode="INDEL",
        max_gaussians=4,
        tranche="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0",
        annotations="-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP",
        mills="mills,known=false,training=true,truth=true,prior=12",
        mills_path=known_variants_mills_path,
        axiompoly="axiomPoly,known=false,training=true,truth=false,prior=10",
        axiompoly_path=known_variants_axiompoly_path,
        dbsnp="dbsnp,known=true,training=false,truth=false,prior=2",
        dbsnp_path=known_variants_dbsnp_path
    log:
        std="%s/variantrecalibrator_indel.log" % log_dir,
        cluster_log="%s/variantrecalibrator_indel.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/variantrecalibrator_indel.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/variantrecalibrator_indel.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["variantrecalibrator_indel_threads"],
        time=config["variantrecalibrator_indel_time"],
        mem=config["variantrecalibrator_indel_mem_mb"],
    threads: config["variantrecalibrator_indel_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' VariantRecalibrator "
        " -V {input.vcf} --trust-all-polymorphic -mode {params.mode} --max-gaussians {params.max_gaussians}"
        " {params.tranche} {params.annotations}"
        " --resource {params.mills}:{params.mills_path}"
        " --resource {params.axiompoly}:{params.axiompoly_path}"
        " --resource {params.dbsnp}:{params.dbsnp_path}"
        " -O {output.recal_table} --tranches-file {output.tranches} > {log.std} 2>&1"

rule variantrecalibrator_snp:
    input:
        vcf=rules.extract_sites.output
    output:
        tranches=joint_snpcall_dir / "all_samples.snp.tranches",
        recal_table=joint_snpcall_dir / "all_samples.snp.recal.table",
    params:
        mode="SNP",
        max_gaussians=6,
        tranche="-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0",
        annotations="-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP",
        hapmap="hapmap,known=false,training=true,truth=true,prior=15",
        hapmap_path=known_variants_hapmap_path,
        omni="omni,known=false,training=true,truth=true,prior=12",
        omni_path=known_variants_omni_path,
        g1000="1000G,known=false,training=true,truth=false,prior=10",
        g1000_path=known_variants_1000g_path,
        dbsnp="dbsnp,known=true,training=false,truth=false,prior=2",
        dbsnp_path=known_variants_dbsnp_path
    log:
        std="%s/variantrecalibrator_snp.log" % log_dir,
        cluster_log="%s/variantrecalibrator_snp.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/variantrecalibrator_snp.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/variantrecalibrator_snp.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["variantrecalibrator_snp_threads"],
        time=config["variantrecalibrator_snp_time"],
        mem=config["variantrecalibrator_snp_mem_mb"],
    threads: config["variantrecalibrator_snp_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' VariantRecalibrator "
        " -V {input.vcf} --trust-all-polymorphic -mode {params.mode} --max-gaussians {params.max_gaussians}"
        " {params.tranche} {params.annotations}"
        " --resource {params.hapmap}:{params.hapmap_path}"
        " --resource {params.omni}:{params.omni_path}"
        " --resource {params.g1000}:{params.g1000_path}"
        " --resource {params.dbsnp}:{params.dbsnp_path}"
        " -O {output.recal_table} --tranches-file {output.tranches} > {log.std} 2>&1"

rule applyvsqr_indel:
    input:
        tranches=rules.variantrecalibrator_indel.output.tranches,
        recal_table=rules.variantrecalibrator_indel.output.recal_table,
        vcf=rules.excess_filter.output
    output:
        joint_snpcall_dir / "all_samples.indel.recalibrated.vcf.gz"
    params:
        mode="INDEL",
        truth_sensitivity_filter_level=99.7,
        create_output_variant_index="true"
    log:
        std="%s/applyvsqr_indel.log" % log_dir,
        cluster_log="%s/applyvsqr_indel.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/applyvsqr_indel.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/applyvsqr_indel.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["applyvsqr_indel_threads"],
        time=config["applyvsqr_indel_time"],
        mem=config["applyvsqr_indel_mem_mb"],
    threads: config["applyvsqr_indel_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' ApplyVQSR "
        " -V {input.vcf} --recal-file {input.recal_table} --tranches-file {input.tranches}"
        " --truth-sensitivity-filter-level {params.truth_sensitivity_filter_level} "
        " --create-output-variant-index {params.create_output_variant_index}"
        " -mode {params.mode} -O {output} > {log.std} 2>&1"

rule applyvsqr_snp:
    input:
        tranches=rules.variantrecalibrator_snp.output.tranches,
        recal_table=rules.variantrecalibrator_snp.output.recal_table,
        vcf=rules.excess_filter.output
    output:
        joint_snpcall_dir / "all_samples.snp.recalibrated.vcf.gz"
    params:
        mode="SNP",
        truth_sensitivity_filter_level=99.7,
        create_output_variant_index="true"
    log:
        std="%s/applyvsqr_snp.log" % log_dir,
        cluster_log="%s/applyvsqr_snp.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/applyvsqr_snp.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/applyvsqr_snp.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["applyvsqr_snp_threads"],
        time=config["applyvsqr_snp_time"],
        mem=config["applyvsqr_snp_mem_mb"],
    threads: config["applyvsqr_snp_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' ApplyVQSR "
        " -V {input.vcf} --recal-file {input.recal_table} --tranches-file {input.tranches}"
        " --truth-sensitivity-filter-level {params.truth_sensitivity_filter_level} "
        " --create-output-variant-index {params.create_output_variant_index}"
        " -mode {params.mode} -O {output} > {log.std} 2>&1"
