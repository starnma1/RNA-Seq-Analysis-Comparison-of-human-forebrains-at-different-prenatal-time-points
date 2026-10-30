#!/usr/bin/env bash
#SBATCH --job-name=star-align
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G # 64 GB
#SBATCH --time=04:00:00
#SBATCH --output=logs/align-%A_%a.log

set -euo pipefail

# --- paths (job is submitted from the alignment folder) ----------------
cd "$SLURM_SUBMIT_DIR"
GENOME_DIR=../index/star-index
FASTQ_DIR=../fastq_prenatal_forebrain/fastq
OUT_DIR=./mapping

# --- environment -------------------------------------------------------
export MAMBA_ROOT_PREFIX="$SLURM_SUBMIT_DIR/../micromamba"
export MAMBA_NO_RC=true
export PATH="$HOME/.local/bin:$PATH"
ENV_NAME=sys-genomics-align           # check with: micromamba env list
eval "$(micromamba shell hook -s bash)"
micromamba activate "$ENV_NAME"

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

STAR --genomeDir "$GENOME_DIR" \
     --runThreadN "$SLURM_CPUS_PER_TASK" \
     --readFilesIn "$fq" \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix "$out/"

samtools index -@ "$SLURM_CPUS_PER_TASK" "$out/Aligned.sortedByCoord.out.bam"