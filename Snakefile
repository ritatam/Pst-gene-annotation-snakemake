import os
from glob import glob

configfile: "config.yaml"

REF = config["ref_fasta"]
HAPLOTYPES = ["A", "B"]

# filename check and extract sample name
SAMPLES = []
for filename in os.listdir(config["input_read_dir"]):
    # check extension
    extensions = {".fastq.gz": 9, ".fastq": 6, ".fq.gz": 6, ".fq": 3}
    for extension, trim_len in extensions.items():
        if filename.endswith(extension):
            sample = filename[:-trim_len]
            SAMPLES.append(sample)
            break
    else:
        raise ValueError("Input read files must be in .fastq/fastq.gz/fq/fq.gz format.")

# extract fields in the sample name 
SAMPLE_IDS, CONDITIONS, REPLICATES = set(), set(), set()
for sample in SAMPLES:
    sample_id, condition, replicate = sample.split("_")
    SAMPLE_IDS.add(sample_id)
    CONDITIONS.add(condition)
    REPLICATES.add(replicate)
assert len(SAMPLE_IDS) == 1, "All samples must have the same sample name. ([sample_id]_[condition]_rep[replicate_number].fastq.gz)"
SAMPLE_ID = SAMPLE_IDS.pop()
CONDITIONS = list(CONDITIONS)
REPLICATES = sorted([int(rep.split("rep")[-1]) for rep in list(REPLICATES)])
print("\nPlease check fields extracted from sample read fastq filename:\n")
print("SAMPLES:", SAMPLES, "\n")
print("SAMPLE_ID:", SAMPLE_ID, "\n")
print("CONDITIONS:", CONDITIONS, "\n")
print("REPLICATES:", REPLICATES, "\n\n")


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
        expand(os.path.join(splice_aln_hap_partitioned_dir, "{sample}." + config["ref_name"] + ".hap{hap}.bam"), sample=SAMPLES, hap=HAPLOTYPES),
        os.path.join(config["outdir"], "all_samples_read_alignment_stats_summary.tsv"),
        expand(os.path.join(splice_aln_hap_partitioned_dir, SAMPLE_ID+".{condition}."+config["ref_name"] + ".hap{hap}.merged.bam"), condition=CONDITIONS, hap=HAPLOTYPES)


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


rule merge_bams_per_condition:
    input:
        expand(os.path.join(splice_aln_hap_partitioned_dir, "{sample}." + config["ref_name"] + ".hap{haplotype}.bam"), sample=SAMPLES, haplotype=HAPLOTYPES)
    output:
        hapA_bam = os.path.join(splice_aln_hap_partitioned_dir, SAMPLE_ID+".{condition}."+config["ref_name"] + ".hapA.merged.bam"),
        hapB_bam = os.path.join(splice_aln_hap_partitioned_dir, SAMPLE_ID+".{condition}."+config["ref_name"] + ".hapB.merged.bam")
    params:
        hapA_bams = os.path.join(splice_aln_hap_partitioned_dir, "*_{condition}_rep*.hapA.bam"),
        hapB_bams = os.path.join(splice_aln_hap_partitioned_dir, "*_{condition}_rep*.hapB.bam")
    envmodules:
        "samtools/1.12"
    shell:
        "samtools merge -@4 {output.hapA_bam} $(ls {params.hapA_bams}); "
        "samtools merge -@4 {output.hapB_bam} $(ls {params.hapB_bams}); "
        "samtools index -@4 {output.hapA_bam}; "
        "samtools index -@4 {output.hapB_bam}"



# rule merge_replicates_per_condition:
#     input:
#         expand(os.path.join(splice_aln_hap_partitioned_dir, "{sample}." + config["ref_name"] + ".hap{haplotype}.bam"), sample=SAMPLES, haplotype=HAPLOTYPES)
#     output:
#         hapA_merged_bam = os.path.join(splice_aln_hap_partitioned_dir,  "{sample_name}_{condition}_rep{rep}."+config["ref_name"] + ".hapA.merged.bam"),
#         hapB_merged_bam = os.path.join(splice_aln_hap_partitioned_dir,  "{sample_name}_{condition}_rep{rep}."+config["ref_name"] + ".hapB.merged.bam"),
#     run:
#         for condition in CONDITIONS:
#             hapA_bams = " ".join([n for n in input if condition in os.path.basename(n) and "hapA" in os.path.basename(n)])
#             hapB_bams = " ".join([n for n in input if condition in os.path.basename(n) and "hapB" in os.path.basename(n)])
#             shell("samtools merge -@4 -o {output.hapA_merged_bam} {hapA_bams}")
#             shell("samtools merge -@4 -o {output.hapB_merged_bam} {hapB_bams}")

        

#################################
##### STRUCTURAL ANNOTATION #####
#################################

# stringtie3_dir = os.path.join(config["outdir"], "stringtie3")
# espresso_dir = os.path.join(config["outdir"], "espresso")

# rule stringtie3_assembly:
#     input:
#         expand(os.path.join(splice_aln_hap_partitioned_dir, "{sample}." + config["ref_name"] + ".hap{haplotype}}.bam"), sample=SAMPLES, haplotype=HAPLOTYPES)
#     output:
#         gtf = os.path.join(stringtie3_dir, "{sample}.stringtie3.gtf")
#     params:
#         stringtie3 = config["stringtie3_bin"]
#     shell:
#         "{params.stringtie3} {input.bam} -p {threads} -G {config['ref_gtf']} -o {output.gtf}"



# stringtie3


# espresso


# gene prediction: codingquarry