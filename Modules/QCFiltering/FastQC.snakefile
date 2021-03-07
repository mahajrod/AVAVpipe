rule fastqc:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "{out_dir}/{fastqc_dir}/{sample}/"
    params:
        kmer=10
    log:
        "{out_dir}/{log_dir}/{sample}/fastqc.log"
    benchmark:
        "{out_dir}/{benchmark_dir}/{sample}/fastqc.benchmark.txt"
    threads: 2
    shell:
        "fastqc --nogroup -k {params.kmer} -t {threads} -o {output} {input} 2> {log}"
