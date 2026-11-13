#!/usr/bin/env bash
#SBATCH --job-name=feature-counting
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/feature-counting-%j.log

set -euo pipefail

cd "$SLURM_SUBMIT_DIR/.."

if ! type module >/dev/null 2>&1; then
    source /cluster/software/profiles/lmod.sh
    source /cluster/software/profiles/module_path.sh
fi
module load stack/2025-06
module load r/4.4.1

Rscript alignment/processing/featureCounting.R
