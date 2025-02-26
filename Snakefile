import os
from glob import glob

configfile: "config.yaml"

REF = config["ref_fasta"]

# filename check and extract sample name
SAMPLES = []
for filename in os.listdir(config["input_read_dir"]):
    # check extension
    if filename.endswith(".fastq.gz"):
        sample = filename[:-9]
    elif filename.endswith(".fastq"):
        sample = filename[:-6]
    elif filename.endswith(".fq.gz"):
        sample = filename[:-6]
    elif filename.endswith(".fq"):
        sample = filename[:-3]
    else:
        raise ValueError("Input files must be in fastq format :(")
    SAMPLES.append(sample)


#####################
##### ALIGNMENT #####
#####################

fastq_stats_dir = os.path.join(config["outdir"], "fastq_stats")
splice_aln_dir = os.path.join(config["outdir"], "splice_aln")
splice_aln_hap_partitioned_dir = os.path.join(splice_aln_dir, "haplotype_partitioned")

rule all:
    input:
        expand(os.path.join(fastq_stats_dir, "{sample}.fastq.stats"), sample=SAMPLES),
        expand(os.path.join(splice_aln_dir, "{sample}." + config["ref_name"] + ".bam"), sample=SAMPLES),
        expand(os.path.join(splice_aln_dir, "{sample}." + config["ref_name"] + ".bamstats"), sample=SAMPLES),
        expand(os.path.join(splice_aln_hap_partitioned_dir, "{sample}." + config["ref_name"] + ".hap{hap}.bam"), sample=SAMPLES, hap=["A","B"]),
        os.path.join(config["outdir"], "all_samples_read_alignment_stats_summary.tsv")


rule seqkit_fastq_quality_stats:
    input:
        fastq = os.path.join(config["input_read_dir"], "{sample}.fastq.gz")
    output:
        stats = os.path.join(fastq_stats_dir, "{sample}.fastq.stats")
    group: "group0"
    shell:
        "set +eu "
        " && . $(conda info --base)/etc/profile.d/conda.sh"
        " && conda activate /g/data/xf3/miniconda/envs/common-tools; "
        "seqkit stats --all -T {input.fastq} > {output.stats}"


rule minimap2_splice_mapping:
    input:
        ref = REF,
        fastq = os.path.join(config["input_read_dir"], "{sample}.fastq.gz")
    output:
        bam = os.path.join(splice_aln_dir, "{sample}." + config["ref_name"] + ".bam")
    envmodules:
        "minimap2/2.28",
        "samtools/1.12"
    threads: 48
    shell:
        "minimap2 -ax splice -ub -G 3000 --secondary=no -t {threads} {input.ref} {input.fastq} |  samtools sort -@ {threads} -O BAM -o {output.bam}; "
        "samtools index -@ {threads} {output.bam}"


rule generate_mapping_stats:
    input:
        bam = os.path.join(splice_aln_dir, "{sample}." + config["ref_name"] + ".bam")
    output:
        stats = os.path.join(splice_aln_dir, "{sample}." + config["ref_name"] + ".bamstats")
    envmodules:
        "samtools/1.12"
    group: "group0"
    shell:
        "samtools flagstat {input.bam} > {output.stats}"


rule aggregate_read_and_mapping_stats:
    input:
        expand(os.path.join(fastq_stats_dir, "{sample}.fastq.stats"), sample=SAMPLES),
        expand(os.path.join(splice_aln_dir, "{sample}." + config["ref_name"] + ".bamstats"), sample=SAMPLES)
    params:
        fastq_stats_dir,
        splice_aln_dir,
        config["ref_name"]
    output:
        os.path.join(config["outdir"], "all_samples_read_alignment_stats_summary.tsv")
    script:
        "scripts/aggregate_read_alignment_stats_reports.py"


rule partition_haplotype_bam:
    input:
        bam = os.path.join(splice_aln_dir, "{sample}." + config["ref_name"] + ".bam")
    output:
        hapA = os.path.join(splice_aln_hap_partitioned_dir, "{sample}." + config["ref_name"] + ".hapA.bam"),
        hapB = os.path.join(splice_aln_hap_partitioned_dir, "{sample}." + config["ref_name"] + ".hapB.bam")
    envmodules:
        "samtools/1.12"
    group: "group0"
    shell:
        "samtools view -bh -F4 -@4 {input.bam} $(for n in {{1..18}}; do echo chr${{n}}A; done) -o {output.hapA}; "
        "samtools view -bh -F4 -@4 {input.bam} $(for n in {{1..18}}; do echo chr${{n}}B; done) -o {output.hapB}; "
        "samtools index -@4 {output.hapA}; "
        "samtools index -@4 {output.hapB}"




#####################
##### ALIGNMENT #####
#####################

