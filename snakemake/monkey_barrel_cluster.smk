# Monkey Barrel Demultiplexer - Clustering Workflow

### Setup ###
configfile: "monkey_barrel_cluster_config.yaml"
workdir: config["workdir"]

import pandas
samples = pandas.read_csv(config["barcodes_file"], header=None, delimiter="\t").loc[:,0].tolist()
#### # ###

### Output ###
rule all:
    input:
        [f"output/6_sorted/{sample}.fasta" for sample in samples]
### # ###

### Workflow ###
rule multimer_demux:
    input:
        config["basecalled_fastq"]
    output:
        expand("output/1_demux/{sample}.fastq", sample=samples)
    params:
        mamba_env=config["mamba_env"],
        barcodes=config["barcodes_file"],
        order=config["barcodes_order"],
        flush=config["flush_size"],
        min_len=config["monomer_min_length"],
        max_len=config["monomer_max_length"],
        e=config["max_edit_dist"],
        d=config["max_deletions"],
        i=config["max_insertions"],
        s=config["max_substitutions"],
    threads:
        config["demux_threads"]
    shell:
        "mamba run -n {params.mamba_env} python multimer_demux.py "
        "{input} "
        "{params.barcodes} "
        "--barcode_order {params.order} "
        "--out_dir output/1_demux "
        "--flush_size {params.flush} "
        "--workers {threads} "
        "--min_length {params.min_len} "
        "--max_length {params.max_len} "
        "--max_edit_dist {params.e} "
        "--max_deletions {params.d} "
        "--max_insertions {params.i} "
        "--max_substitutions {params.s} "
        "--trim"


rule derep:
    input:
        "output/1_demux/{sample}.fastq"
    output:
        "output/2_derep/{sample}.fastq"
    params:
        mamba_env=config["mamba_env"],
        q_max=config["q_max"],
        fasta_width=config["fasta_width"]
    threads:
        1
    shell:        
        "mamba run -n {params.mamba_env} vsearch "
        "--fastx_uniques {input} "
        "--fastqout {output} "
        "--fastq_qout_max "
        "--fastq_qmaxout {params.q_max} "
        "--sizeout "
        "--fasta_width {params.fasta_width}"

rule fastq_filter:
    input:
        "output/2_derep/{sample}.fastq"
    output: 
        "output/3_filter/{sample}.fasta"
    params:
        mamba_env=config["mamba_env"],
        max_ee=config["max_expected_error"],
        min_len=config["min_len"],
        max_len=config["max_len"],
        max_ns=config["max_ns"],
        q_max=config["q_max"],
        fasta_width=config["fasta_width"]
    threads:
        1
    shell:
        "mamba run -n {params.mamba_env} vsearch "
        "--fastq_filter {input} "
        "--fastq_maxee {params.max_ee} "
        "--fastq_minlen {params.min_len} "
        "--fastq_maxlen {params.max_len} "
        "--fastq_maxns {params.max_ns} "
        "--fastq_qmax {params.q_max} "
        "--fastaout {output} "
        "--fasta_width {params.fasta_width}"

rule cluster_or_denoise:
    input:
        "output/3_filter/{sample}.fasta"
    output:
        "output/4_cluster/{sample}.fasta"
    params:
        mamba_env=config["mamba_env"],
        cluster_id=config["cluster_id"],
        iddef=config["cluster_id_definition"],
        fasta_width=config["fasta_width"]
    threads:
        config["cluster_threads"]
    shell:
        "mamba run -n {params.mamba_env} vsearch " 
        "--cluster_size {input} "
        "--id {params.cluster_id} "
        "--iddef {params.iddef} "
        "--threads {threads} "
        "--sizein "
        "--sizeout "
        "--fasta_width {params.fasta_width} "
        "--consout {output}"

rule uchime3_denovo:
    input:
        "output/4_cluster/{sample}.fasta"
    output:
        "output/5_nonchimeras/{sample}.fasta"
    params:
        mamba_env=config["mamba_env"],
        abskew=config["abundance_skew"],
        fasta_width=config["fasta_width"]
    threads:
        1
    shell:
        "mamba run -n {params.mamba_env} vsearch " 
        "--uchime3_denovo {input} "
        "--abskew {params.abskew} "
        "--sizein "
        "--sizeout "
        "--fasta_width {params.fasta_width} "
        "--qmask none "
        "--nonchimeras {output}"

rule size_sort:
    input:
        "output/5_nonchimeras/{sample}.fasta"
    output:
        "output/6_sorted/{sample}.fasta"
    params:
        mamba_env=config["mamba_env"],
        fasta_width=config["fasta_width"],
        minsize=config["size_sort_minsize"]
    shell:
        "mamba run -n {params.mamba_env} vsearch " 
        "--sortbysize {input} "
        "--sizein "
        "--sizeout "
        "--fasta_width {params.fasta_width} "
        "--minsize {params.minsize} "
        "--output {output}"
