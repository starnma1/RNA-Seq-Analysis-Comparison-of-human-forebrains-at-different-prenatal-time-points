#!/usr/bin/env bash

PREPROCESSING_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PROJECT_ROOT=$(cd -- "$PREPROCESSING_DIR/.." && pwd)
PIXI_BIN=${PIXI_BIN:-$(command -v pixi || true)}
if [[ -z "$PIXI_BIN" || ! -x "$PIXI_BIN" ]]; then
    echo "Error: Pixi not found. Put it on PATH or set PIXI_BIN to its executable." >&2
    exit 1
fi
PIXI_BIN="$(cd -- "$(dirname -- "$PIXI_BIN")" && pwd)/$(basename -- "$PIXI_BIN")"
if [[ ! -d "$PREPROCESSING_DIR/.pixi/envs/default" ||
      ! -f "$PREPROCESSING_DIR/.pixi/preprocessing-ready.sha256" ]] ||
   ! (cd "$PREPROCESSING_DIR" && sha256sum --check --status .pixi/preprocessing-ready.sha256); then
    echo "Error: preprocessing Pixi environment is missing or out of date; submit preprocessing/scripts/submit-install.sh first." >&2
    exit 1
fi
export PROJECT_ROOT PIXI_BIN
