#!/usr/bin/env bash
#SBATCH --job-name=star-index
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G              # 40 GB total; human genome needs ~32 GB
#SBATCH --time=04:00:00
#SBATCH --output=build-index-%j.log   # written in the folder you submit from

set -euo pipefail

: "${PROJECT_ROOT:?Submit using the preprocessing index launcher}"
: "${PIXI_BIN:?Submit using the preprocessing index launcher}"
cd "$PROJECT_ROOT/index"
mkdir -p star-index
rm -rf _STARtmp                       # STAR refuses to start if a leftover one exists

"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline STAR \
     --runThreadN "$SLURM_CPUS_PER_TASK" \
     --runMode genomeGenerate \
     --genomeDir star-index \
     --genomeFastaFiles ./reference_genome/hg38.fa