#!/usr/bin/env bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fastqc-%A_%a.log

# STEP=fastqc (default): one array task per FASTQ file
# STEP=multiqc         : combine all FastQC reports (submitted by the launcher)
set -euo pipefail

: "${PROJECT_ROOT:?Submit using the preprocessing QC launcher}"
: "${PIXI_BIN:?Submit using the preprocessing QC launcher}"
FASTQ_DIR="$PROJECT_ROOT/raw/fastq"
FASTQC_DIR="$PROJECT_ROOT/library_qc/fastqc"
MULTIQC_DIR="$PROJECT_ROOT/library_qc/multiqc"
STEP=${STEP:-fastqc}

if [[ "$STEP" == "multiqc" ]]; then
    "$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline multiqc \
        "$FASTQC_DIR" --outdir "$MULTIQC_DIR" --filename fastq_qc_report --force
    exit 0
fi

FILES=("$FASTQ_DIR"/*.fastq.gz)
fq="${FILES[$SLURM_ARRAY_TASK_ID]:-}"
if [[ -z "$fq" ]]; then
    echo "No file for array index $SLURM_ARRAY_TASK_ID (found ${#FILES[@]} files)" >&2
    exit 1
fi
sample=$(basename "$fq" .fastq.gz)

if [[ -s "$FASTQC_DIR/${sample}_fastqc.zip" ]]; then
    echo "Skipping $sample (already done)"
    exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: FastQC on $sample"
"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline fastqc \
    --threads 1 --outdir "$FASTQC_DIR" "$fq"