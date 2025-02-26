#!/bin/bash
#PBS -P xf3
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=1
#PBS -l mem=2G
#PBS -l storage=scratch/xf3+gdata/xf3+scratch/xe2+gdata/xe2
#PBS -l wd
#PBS -j oe
#PBS -m abe

###change settings here###
script_dir=/scratch/xe2/ht5438/198e_anno/snakemake
conda_snakemake_env=/g/data/xe2/ht5438/conda_env/snakemake
storage=scratch/xf3+gdata/xf3+scratch/xe2+gdata/xe2
##########################

source ${script_dir}/gadimod.sh
conda activate ${conda_snakemake_env}

set -ueo pipefail
logdir=pbs_log
mkdir -p $logdir
export TMPDIR=${PBS_JOBFS:-$TMPDIR}
TARGET=${TARGET:-all}

QSUB="qsub -q {cluster.queue} -l ncpus={cluster.ncpus} -l ngpus={cluster.ngpus}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem} -N {cluster.name} -l storage={cluster.storage}"
QSUB="$QSUB -l wd -m abe -j oe -o $logdir -P {cluster.project}" 

snakemake																	    \
	-j 1000																	    \
	--max-jobs-per-second 2													    \
	--cluster-config ${script_dir}/cluster.yaml		\
	--local-cores ${PBS_NCPUS:-1}											    \
	--js ${script_dir}/jobscript.sh			    	\
	--nolock																    \
	--keep-going															    \
	--rerun-incomplete														    \
	--use-envmodules														    \
	--cluster "$QSUB"														    \
	--use-conda																    \
	"$TARGET"