module purge
export TMPDIR=${PBS_JOBFS:-/tmp}

function useconda() {
	                eval "$(/g/data/xf3/miniconda/bin/conda shell.zsh hook)"
			        }


useconda

source /g/data/xf3/miniconda/etc/profile.d/conda.sh