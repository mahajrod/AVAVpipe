rule jellyfish_count:
    input:
        forward_reads=rules.trimmomatic.output.pe_forward, #"%s/{sample_id}/{sample_id}.trimmed_1.fastq.gz" % filtered_read_dir,
        reverse_reads=rules.trimmomatic.output.pe_reverse,#"%s/{sample_id}/{sample_id}.trimmed_2.fastq.gz" % filtered_read_dir
    output:
        kmer_dir_path / ("{sample_id}/{sample_id}.%i.jf" % config["jellyfish_kmer_length"])
    params:
        kmer_length=config["jellyfish_kmer_length"]
    log:
        std=log_dir_path / "{sample_id}/jellyfish_count.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.jellyfish_count.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.jellyfish_count.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/jellyfish_count.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["jellyfish_count_threads"],
        time=config["jellyfish_count_time"],
        mem=config["jellyfish_count_mem_mb"],
    threads:
        config["jellyfish_count_threads"]
    shell:
         "jellyfish count -C -m {params.kmer_length} -s 20G -t {threads} -o {output}  "
         "<(gunzip -c {input.forward_reads}) <(gunzip -c {input.reverse_reads}) 2>{log.std}"

rule jellyfish_histo:
    input:
        rules.jellyfish_count.output
    output:
        kmer_dir_path / ("{sample_id}/{sample_id}.%i.histo" % config["jellyfish_kmer_length"])
    params:
        min_coverage=1,
        max_coverage=100000000,
        increment=1
    log:
        std=log_dir_path / "{sample_id}/jellyfish_histo.log",
        cluster_log=cluster_log_dir_path / "{sample_id}.jellyfish_histo.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample_id}.jellyfish_histo.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample_id}/jellyfish_histo.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["jellyfish_histo_threads"],
        time=config["jellyfish_histo_time"],
        mem=config["jellyfish_histo_mem_mb"],
    threads:
        config["jellyfish_histo_threads"]
    shell:
         "jellyfish histo -o {output} -t {threads} -l {params.min_coverage} -h {params.max_coverage} -i {params.increment} {input}"