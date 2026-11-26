#!/usr/bin/env bash
#SBATCH --job-name=feature-counting
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/feature-counting-%j.log

set -euo pipefail

: "${PROJECT_ROOT:?Submit using the preprocessing counting launcher}"
: "${PIXI_BIN:?Submit using the preprocessing counting launcher}"
cd "$PROJECT_ROOT"
exec "$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline \
    Rscript "$PROJECT_ROOT/preprocessing/scripts/counting/featureCounting.R"
