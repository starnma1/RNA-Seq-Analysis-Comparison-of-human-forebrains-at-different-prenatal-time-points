#!/usr/bin/env bash
# Usage (from anywhere): bash preprocessing/scripts/qc/submit-fastqc.sh
# Submits FastQC only; it can be run independently of alignment and counting.
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
mkdir -p logs fastqc

N=${#files[@]}
echo "Submitting FastQC for $N files"

sbatch --array=0-$((N-1))%20 "$SCRIPT_DIR/slurm-fastqc.sh"