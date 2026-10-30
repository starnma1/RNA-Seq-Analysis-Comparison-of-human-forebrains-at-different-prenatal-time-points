#!/usr/bin/env bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fastqc-%A_%a.log

# STEP=fastqc (default): one array task per FASTQ file
# STEP=multiqc         : combine all FastQC reports (submitted by the launcher)
set -euo pipefail

cd "$SLURM_SUBMIT_DIR"                       # the quality_control folder
FASTQ_DIR=../raw/fastq
FASTQC_DIR=./fastqc
MULTIQC_DIR=./multiqc
STEP=${STEP:-fastqc}

export MAMBA_ROOT_PREFIX="$SLURM_SUBMIT_DIR/../micromamba"
export MAMBA_NO_RC=true
export PATH="$HOME/.local/bin:$PATH"
ENV_NAME=sys-genomics-align                  # check with: micromamba env list
eval "$(micromamba shell hook -s bash)"
micromamba activate "$ENV_NAME"

if [[ "$STEP" == "multiqc" ]]; then
    multiqc "$FASTQC_DIR" --outdir "$MULTIQC_DIR" --filename fastq_qc_report --force
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
fastqc --threads 1 --outdir "$FASTQC_DIR" "$fq"