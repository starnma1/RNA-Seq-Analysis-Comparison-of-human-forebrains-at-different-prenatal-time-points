#!/usr/bin/env bash
#SBATCH --job-name=star-index
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G              # 40 GB total; human genome needs ~32 GB
#SBATCH --time=04:00:00
#SBATCH --output=build-index-%j.log   # written in the folder you submit from

set -euo pipefail

# --- environment -------------------------------------------------------
export MAMBA_ROOT_PREFIX=/cluster/scratch/$USER/sys-genomics/micromamba
export MAMBA_NO_RC=true
export PATH="$HOME/.local/bin:$PATH"
ENV_NAME=sys-genomics-align           # check with: micromamba env list
eval "$(micromamba shell hook -s bash)"
micromamba activate "$ENV_NAME"

# --- build -------------------------------------------------------------
cd "$SLURM_SUBMIT_DIR"                # run from the index folder
rm -rf _STARtmp                       # STAR refuses to start if a leftover one exists

STAR --runThreadN "$SLURM_CPUS_PER_TASK" \
     --runMode genomeGenerate \
     --genomeDir star-index \
     --genomeFastaFiles ./reference_genome/hg38.fa