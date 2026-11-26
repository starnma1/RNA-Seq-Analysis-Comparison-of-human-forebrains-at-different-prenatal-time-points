#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PREPROCESSING_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)
PROJECT_ROOT=$(cd -- "$PREPROCESSING_DIR/.." && pwd)
PIXI_BIN=${PIXI_BIN:-$(command -v pixi || true)}
if [[ -z "$PIXI_BIN" || ! -x "$PIXI_BIN" ]]; then
    echo "Error: Pixi not found. Put it on PATH or set PIXI_BIN to its executable." >&2
    exit 1
fi
PIXI_BIN="$(cd -- "$(dirname -- "$PIXI_BIN")" && pwd)/$(basename -- "$PIXI_BIN")"
if [[ ! -f "$PREPROCESSING_DIR/pixi.lock" ]]; then
    echo "Error: $PREPROCESSING_DIR/pixi.lock is missing." >&2
    exit 1
fi
export PROJECT_ROOT PIXI_BIN
mkdir -p "$PREPROCESSING_DIR/logs"
cd "$PREPROCESSING_DIR"
sbatch "$SCRIPT_DIR/slurm-install.sh"
