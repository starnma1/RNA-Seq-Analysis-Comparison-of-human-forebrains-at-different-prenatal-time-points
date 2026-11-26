#!/usr/bin/env bash
#SBATCH --job-name=pixi-install
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --output=logs/pixi-install-%j.log

set -euo pipefail
: "${PROJECT_ROOT:?Submit using preprocessing/scripts/submit-install.sh}"
: "${PIXI_BIN:?Submit using preprocessing/scripts/submit-install.sh}"
: "${SLURM_CPUS_PER_TASK:?Slurm CPU allocation is required}"

export RAYON_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export TOKIO_WORKER_THREADS="$SLURM_CPUS_PER_TASK"
cd "$PROJECT_ROOT/preprocessing"
mkdir -p .pixi
rm -f .pixi/preprocessing-ready.sha256
./pixi.sh install --locked --concurrent-downloads "$SLURM_CPUS_PER_TASK" --concurrent-solves 1
./pixi.sh run --as-is --offline bash -c '
    set -e
    STAR --version
    samtools --version | head -n 1
    fastqc --version
    multiqc --version
    Rscript -e "library(Rsubread); library(rtracklayer)"
'
sha256sum pixi.toml pixi.lock > .pixi/preprocessing-ready.sha256.tmp
mv .pixi/preprocessing-ready.sha256.tmp .pixi/preprocessing-ready.sha256
