localrules: create_sample_file

rule haplotypecaller_gvcf:
    input:
        region="%s/intervals/region_{region_id}.list" % reference_region_dir_path,
        bam=rules.applybsqr.output.bam,
        bai=rules.applybsqr.output.bai,
        reference=config["reference"],
    output:
        gvcf=temp("%s/{sample_id}/haplotypecaller_gvcf/{sample_id}.region_{region_id}.gvcf" % snpcall_dir)
    params:
        ""
    log:
        std="%s/{sample_id}/haplotypecaller_gvcf/haplotypecaller_gvcf.region_{region_id}.log" % log_dir,
        cluster_log="%s/haplotypecaller_gvcf/{sample_id}.haplotypecaller_gvcf.region_{region_id}.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/haplotypecaller_gvcf/{sample_id}.haplotypecaller_gvcf.region_{region_id}.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/haplotypecaller_gvcf/haplotypecaller_gvcf.region_{region_id}.benchmark.txt" % benchmark_dir
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
        lambda wildcards:  expand("%s/{sample_id}/haplotypecaller_gvcf/{sample_id}.region_{region_id}.gvcf" % snpcall_dir,
                                  region_id=glob_wildcards("%s/intervals/region_{region_id}.list" % reference_region_dir_path)[0],
                                  sample_id=[wildcards.sample_id])
    output:
        "%s/{sample_id}/{sample_id}.gvcf" % snpcall_dir
    params:
        input_files="%s/{sample_id}/haplotypecaller_gvcf/*.gvcf" % snpcall_dir,
        splited_gvcf_list="%s/{sample_id}/{sample_id}.sorted.mkdup.recal.table.list" % snpcall_dir,
        reference_dict=reference_dict_path
    log:
        std="%s/{sample_id}.merge_splited_gvcf.log" % log_dir,
        cluster_log="%s/{sample_id}.merge_splited_gvcf.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/{sample_id}.merge_splited_gvcf.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/{sample_id}/merge_splited_gvcf.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["merge_splited_gvcf_threads"],
        time=config["merge_splited_gvcf_time"],
        mem=config["merge_splited_gvcf_mem_mb"],
    threads: config["merge_splited_gvcf_threads"]
    shell:
        "ls {params.input_files} | sort -V > {params.splited_gvcf_list}; "
        " picard --java-options '-Xmx{resources.mem}m' MergeVcfs -I {params.splited_gvcf_list} -O {output} "
        " -D {params.reference_dict}> {log.std} 2>&1"

rule create_sample_file:
    input:
         expand("%s/{sample_id}/{sample_id}.gvcf" % snpcall_dir, sample_id=config["sample_list"])
    output:
        sample_file=joint_snpcall_dir / "sample_file.tsv"
    run:
        with open(output.sample_file, "w") as out_fd:
            for sample in config["sample_list"]:
                out_fd.write("{0}\t{1}/{0}/{0}.gvcf\n".format(sample, str(snpcall_dir)))

rule genomicsdbimport:
    input:
        gvcfs=expand("%s/{sample_id}/{sample_id}.gvcf" % snpcall_dir, sample_id=config["sample_list"]),
        sample_file=rules.create_sample_file.output
    output:
        directory(joint_snpcall_dir / "gvcf_database")
    params:
        batch_size=50,
        reader_threads=config["genomicsdbimport_reader_threads"],
        interval_threads=config["genomicsdbimport_interval_threads"],
        interval_list=reference_genotyping_whitelist_path
    log:
        std="%s/genomicsdbimport.log" % log_dir,
        cluster_log="%s/genomicsdbimport.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/genomicsdbimport.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/genomicsdbimport.benchmark.txt" % benchmark_dir
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
        " -L {params.interval_list}> {log.std} 2>&1"
