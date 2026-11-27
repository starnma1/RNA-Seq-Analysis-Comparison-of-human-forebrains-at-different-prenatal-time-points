#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/../common.sh"
shopt -s nullglob
files=("$PROJECT_ROOT"/raw/fastq/*.fastq.gz)
if (( ${#files[@]} == 0 )); then
    echo "Error: no FASTQs found in $PROJECT_ROOT/raw/fastq." >&2
    exit 1
fi
cd "$PROJECT_ROOT/alignment"
mkdir -p logs
export FORCE_REALIGN="${FORCE_REALIGN:-0}"
sbatch --export=ALL,FORCE_REALIGN="$FORCE_REALIGN" \
    --array=0-$((${#files[@]}-1))%5 "$SCRIPT_DIR/slurm-align-fastq.sh"
