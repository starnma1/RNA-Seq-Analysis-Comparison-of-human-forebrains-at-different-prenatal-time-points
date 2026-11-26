#!/usr/bin/env bash
# Usage (from anywhere): bash preprocessing/scripts/qc/submit-quality-control.sh
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/../common.sh"
shopt -s nullglob
files=("$PROJECT_ROOT"/raw/fastq/*.fastq.gz)
if (( ${#files[@]} == 0 )); then
    echo "Error: no FASTQs found in $PROJECT_ROOT/raw/fastq." >&2
    exit 1
fi
cd "$PROJECT_ROOT/library_qc"
mkdir -p logs fastqc multiqc

N=${#files[@]}
echo "Submitting FastQC for $N files"

jid=$(sbatch --parsable --array=0-$((N-1))%20 "$SCRIPT_DIR/slurm-quality-control.sh")

# MultiQC runs once the whole array has finished (even if a few tasks failed)
sbatch --dependency=afterany:"$jid" \
       --export=ALL,STEP=multiqc \
       --job-name=multiqc \
       --time=00:30:00 \
       --output=logs/multiqc-%j.log \
       "$SCRIPT_DIR/slurm-quality-control.sh"