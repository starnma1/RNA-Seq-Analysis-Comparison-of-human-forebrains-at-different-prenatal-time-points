#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/../common.sh"
cd "$PROJECT_ROOT/alignment"
mkdir -p logs
sbatch "$SCRIPT_DIR/slurm-package-outputs.sh"
