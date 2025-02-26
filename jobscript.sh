#!/bin/bash

source /scratch/xe2/ht5438/198e_anno/snakemake/gadimod.sh

export TMPDIR=$PBS_JOBFS

set -ueo pipefail
{exec_job}