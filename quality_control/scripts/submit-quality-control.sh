#!/usr/bin/env bash
# Usage (from anywhere): bash quality_control/scripts/submit-quality-control.sh
set -euo pipefail

cd "$(dirname "$0")/.."                      # the quality_control folder
mkdir -p logs fastqc multiqc

N=$(ls ../fastq_prenatal_forebrain/fastq/*.fastq.gz | wc -l)
echo "Submitting FastQC for $N files"

jid=$(sbatch --parsable --array=0-$((N-1))%20 ./scripts/slurm-quality-control.sh)

# MultiQC runs once the whole array has finished (even if a few tasks failed)
sbatch --dependency=afterany:"$jid" \
       --export=ALL,STEP=multiqc \
       --job-name=multiqc \
       --time=00:30:00 \
       --output=logs/multiqc-%j.log \
       ./scripts/slurm-quality-control.sh