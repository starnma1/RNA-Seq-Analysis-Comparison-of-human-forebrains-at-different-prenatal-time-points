#!/usr/bin/env bash
#SBATCH --job-name=star-align
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G # 64 GB
#SBATCH --time=04:00:00
#SBATCH --output=logs/align-%A_%a.log

set -euo pipefail

: "${PROJECT_ROOT:?Submit using the preprocessing alignment launcher}"
: "${PIXI_BIN:?Submit using the preprocessing alignment launcher}"
GENOME_DIR="$PROJECT_ROOT/index/star-index"
FASTQ_DIR="$PROJECT_ROOT/raw/fastq"
OUT_DIR="$PROJECT_ROOT/alignment/mapping"

# --- pick this task's file ---------------------------------------------
FILES=("$FASTQ_DIR"/*.fastq.gz)       # glob is sorted, so indices are stable
fq="${FILES[$SLURM_ARRAY_TASK_ID]:-}"
if [[ -z "$fq" ]]; then
    echo "No file for array index $SLURM_ARRAY_TASK_ID (found ${#FILES[@]} files)" >&2
    exit 1
fi
sample=$(basename "$fq" .fastq.gz)
out="$OUT_DIR/$sample"

if [[ -s "$out/Aligned.sortedByCoord.out.bam" ]]; then
    echo "Skipping $sample (already aligned)"
    exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: aligning $sample"
mkdir -p "$out"
rm -rf "$out/_STARtmp"                # leftover from a failed attempt

"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline STAR --genomeDir "$GENOME_DIR" \
     --runThreadN "$SLURM_CPUS_PER_TASK" \
     --readFilesIn "$fq" \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix "$out/"

"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline samtools index \
    -@ "$SLURM_CPUS_PER_TASK" "$out/Aligned.sortedByCoord.out.bam"