localrules: create_sample_file


rule haplotypecaller_gvcf:
    input:
        region=reference_region_dir_path / "intervals/region_{region_id}.list",
        bam=rules.applybsqr.output.bam,
        bai=rules.applybsqr.output.bai,
        reference=config["reference"],
    output:
        gvcf=temp(snpcall_dir_path / "{sample_id}/haplotypecaller_gvcf/{sample_id}.region_{region_id}.gvcf")
    log:
        std=log_dir_path / "{sample_id}/haplotypecaller_gvcf/haplotypecaller_gvcf.region_{region_id}.log",
        cluster_log=cluster_log_dir_path / "/haplotypecaller_gvcf/{sample_id}.haplotypecaller_gvcf.region_{region_id}.cluster.log",
        cluster_err=cluster_log_dir_path / "haplotypecaller_gvcf/{sample_id}.haplotypecaller_gvcf.region_{region_id}.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/haplotypecaller_gvcf/haplotypecaller_gvcf.region_{region_id}.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["haplotypecaller_threads"],
        time=config["haplotypecaller_time"],
        mem=config["haplotypecaller_mem_mb"],
    threads: config["haplotypecaller_threads"]
    shell:
        "gatk --java-options '-Xmx{resources.mem}m' HaplotypeCaller -R {input.reference} -L {input.region} "
        " -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation"
        " -I {input.bam} -O {output} > {log.std} 2>&1"

rule merge_splited_gvcf:
    input:
        lambda wildcards:  expand("%s/{sample_id}/haplotypecaller_gvcf/{sample_id}.region_{region_id}.gvcf" % snpcall_dir_path,
                                  region_id=glob_wildcards("%s/intervals/region_{region_id}.list" % reference_region_dir_path)[0],
                                  sample_id=[wildcards.sample_id])
    output:
        snpcall_dir_path / "{sample_id}/{sample_id}.gvcf"
    params:
        input_files=snpcall_dir_path / "{sample_id}/haplotypecaller_gvcf/*.gvcf",
        splited_gvcf_list=snpcall_dir_path / "{sample_id}/{sample_id}.splited_gvf_list",
        reference_dict=reference_dict_path
    log:
        std=log_dir_path / "{sample_id}.merge_splited_gvcf.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.merge_splited_gvcf.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.merge_splited_gvcf.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/merge_splited_gvcf.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["merge_splited_gvcf_threads"],
        time=config["merge_splited_gvcf_time"],
        mem=config["merge_splited_gvcf_mem_mb"],
    threads: config["merge_splited_gvcf_threads"]
    shell:
        "ls {params.input_files} | sort -V > {params.splited_gvcf_list}; "
        " workflow/scripts/combine_same_sample_vcf.py -f {params.splited_gvcf_list} -o {output} > {log.std} 2>&1"

rule index_merged_gvcf:
    input:
        rules.merge_splited_gvcf.output
    output:
        snpcall_dir_path / "{sample_id}/{sample_id}.gvcf.idx"
    log:
        std=log_dir_path / "{sample_id}.index_merged_gvcf.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.index_merged_gvcf.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.index_merged_gvcf.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/index_merged_gvcf.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["index_merged_gvcf_threads"],
        time=config["index_merged_gvcf_time"],
        mem=config["index_merged_gvcf_mem_mb"],
    threads: config["index_merged_gvcf_threads"]
    shell:
        "gatk --java-options '-Xmx{resources.mem}m' IndexFeatureFile -I ${input} > {log.std} 2>&1"

rule create_sample_file:
    input:
         expand("%s/{sample_id}/{sample_id}.gvcf" % snpcall_dir_path, sample_id=config["sample_list"])
    output:
        sample_file=joint_snpcall_dir_path / "sample_file.tsv"
    run:
        with open(output.sample_file, "w") as out_fd:
            for sample in config["sample_list"]:
                out_fd.write("{0}\t{1}/{0}/{0}.gvcf\n".format(sample, str(snpcall_dir_path)))

rule genomicsdbimport:
    input:
        gvcfs=expand("%s/{sample_id}/{sample_id}.gvcf" % snpcall_dir_path, sample_id=config["sample_list"]),
        gvcf_indexes=expand("%s/{sample_id}/{sample_id}.gvcf.idx" % snpcall_dir_path, sample_id=config["sample_list"]),
        sample_file=rules.create_sample_file.output,
        interval_file=rules.prepare_genotyping_whitelist_intervals.output
    output:
        directory(joint_snpcall_dir_path / "gvcf_database/callset.json")
    params:
        batch_size=50,
        reader_threads=config["genomicsdbimport_reader_threads"],
        interval_threads=config["genomicsdbimport_interval_threads"],
    log:
        std=log_dir_path / "genomicsdbimport.log",
        cluster_log=cluster_log_dir_path / "genomicsdbimport.cluster.log",
        cluster_err=cluster_log_dir_path / "genomicsdbimport.cluster.err"
    benchmark:
        benchmark_dir_path / "/genomicsdbimport.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["genomicsdbimport_reader_threads"] + config["genomicsdbimport_interval_threads"],
        time=config["genomicsdbimport_time"],
        mem=config["genomicsdbimport_mem_mb"],
    threads: config["genomicsdbimport_reader_threads"] + config["genomicsdbimport_interval_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' GenomicsDBImport --batch-size {params.batch_size} "
        " --sample-name-map {input.sample_file} "
        " --max-num-intervals-to-import-in-parallel {params.interval_threads}"
        " --reader-threads {params.reader_threads}"
        " --genomicsdb-workspace-path {output} "
        " -L {input.interval_file} > {log.std} 2>&1"

rule genotypegvcfs:
    input:
        database=rules.genomicsdbimport.output,
        reference=reference_path
    output:
        joint_snpcall_dir_path / "all_samples.vcf.gz"
    log:
        std=log_dir_path / "genotypegvcfs.log",
        cluster_log=cluster_log_dir_path / "genotypegvcfs.cluster.log",
        cluster_err=cluster_log_dir_path / "genotypegvcfs.cluster.err"
    benchmark:
        benchmark_dir_path / "genotypegvcfs.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["genotypegvcfs_threads"],
        time=config["genotypegvcfs_time"],
        mem=config["genotypegvcfs_mem_mb"],
    threads: config["genotypegvcfs_threads"]
    shell:
        " gatk --java-options '-Xmx{resources.mem}m' GenotypeGVCFs -R {input.reference} "
        " -G StandardAnnotation -G AS_StandardAnnotation"
        " -V gendb://{input.database}"
        " -O {output} > {log.std} 2>&1"
