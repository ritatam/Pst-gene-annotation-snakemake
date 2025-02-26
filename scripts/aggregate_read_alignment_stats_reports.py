from glob import glob
import os
import pandas as pd


fastq_stats_dir = snakemake.params[0]
splice_aln_dir = snakemake.params[1]
ref_name = snakemake.params[2]
output_file = snakemake.output[0]

# for testing
# outdir = "/scratch/xe2/ht5438/198e_anno/output"
# fastq_stats_dir = os.path.join(outdir, "fastq_stats")
# splice_aln_dir = os.path.join(outdir, "splice_aln")
# ref_name = "Pst198E"
# output_file = "test.tsv"

stats_files = sorted(glob(fastq_stats_dir + "/*.fastq.stats"))
sample_names = [os.path.basename(n)[:-12] for n in stats_files]

# merge all read fastq stats into output file
with open(output_file, "w") as f:
    for i, file in enumerate(stats_files):
        with open(file, "r") as stats:
            lines = stats.readlines()
            if i == 0:
                f.write(lines[0])
            f.write(lines[1])

# get sample names for bamstats lookup later
readstats_df = pd.read_csv(output_file, sep="\t")
readstats_df["sample"] = readstats_df["file"].apply(lambda x: os.path.basename(x)[:-9])
readstats_df.rename(columns={"N50_num":"L50", "Q1":"q1_len", "Q2":"q2_len", "Q3":"q3_len"}, inplace=True)

# merge all bamstats into a dictionary
bamstats_df = {}
for sample in sample_names:
    with open(os.path.join(splice_aln_dir, f"{sample}.{ref_name}.bamstats"), "r") as bamstats:
        lines = bamstats.readlines()
        bamstats_df[sample] = {}
        total_read_count = int(lines[0].split(" ")[0])
        mapped_read_count = int(lines[4].split(" ")[0])
        mapped_read_perc = round((mapped_read_count/total_read_count)*100, 2)
        bamstats_df[sample]["sample"] = sample
        bamstats_df[sample]["#total_read"] = total_read_count
        bamstats_df[sample]["#mapped_read"] = mapped_read_count
        bamstats_df[sample]["%mapped_reads"] = mapped_read_perc
bamstats_df = pd.DataFrame.from_dict(bamstats_df, orient="index")

df = readstats_df[["sample", "num_seqs", "sum_len", "min_len", "max_len", "avg_len", "q1_len", "q2_len", "q3_len", "N50", "L50", "Q20(%)", "Q30(%)", "AvgQual", "GC(%)"]]
df["ref_name"] = ref_name

# merge bamstats with fastqstats
df = df.merge(bamstats_df, on="sample", how="inner")
df.to_csv(output_file, sep="\t", index=False)