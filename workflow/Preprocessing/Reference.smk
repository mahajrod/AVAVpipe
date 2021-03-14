rule ref_faidx:
    input:
        config["reference"]
    output:
        "%s.fai" % config["reference"]
    log:
        std="%s/faidx.log" % log_dir,
        cluster_log="%s/faidx.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/faidx.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/faidx.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["faidx_threads"],
        time=config["faidx_time"],
        mem=config["faidx_mem_mb"],
    threads:
        config["faidx_threads"]
    shell:
         "samtools faidx {input} > {log.std} 2>&1"

rule ref_dict:
    input:
        config["reference"]
    output:
        "%s.dict" % (os.path.splitext(config["reference"])[0])
    log:
        std="%s/dict.log" % log_dir,
        cluster_log="%s/dict.cluster.log" % config["cluster_log_dir"],
        cluster_err="%s/dict.cluster.err" % config["cluster_log_dir"]
    benchmark:
        "%s/dict.benchmark.txt" % benchmark_dir
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["dict_threads"],
        time=config["dict_time"],
        mem=config["dict_mem_mb"],
    threads:
        config["dict_threads"]
    shell:
         "picard CreateSequenceDictionary -R {input} > {log.std} 2>&1"