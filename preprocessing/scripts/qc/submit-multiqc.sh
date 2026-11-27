#!/usr/bin/env bash
# Usage (from anywhere): bash preprocessing/scripts/qc/submit-multiqc.sh
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/../common.sh"

"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline python \
    "$SCRIPT_DIR/prepare-multiqc-inputs.py" \
    --project-root "$PROJECT_ROOT" --check-only

cd "$PROJECT_ROOT/library_qc"
mkdir -p logs multiqc
sbatch "$SCRIPT_DIR/slurm-multiqc.sh"
