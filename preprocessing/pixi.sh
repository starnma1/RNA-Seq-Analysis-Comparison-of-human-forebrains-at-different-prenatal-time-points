#!/usr/bin/env bash
set -euo pipefail

PREPROCESSING_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PIXI_BIN=${PIXI_BIN:-$(command -v pixi || true)}
if [[ -z "$PIXI_BIN" || ! -x "$PIXI_BIN" ]]; then
    echo "Error: Pixi not found. Put it on PATH or set PIXI_BIN to its executable." >&2
    exit 1
fi
if (( $# == 0 )); then
    echo "Usage: $0 <pixi-command> [arguments...]" >&2
    exit 1
fi

RUNTIME_DIR="$PREPROCESSING_DIR/.runtime"
mkdir -p "$RUNTIME_DIR/home" "$RUNTIME_DIR/cache" "$RUNTIME_DIR/tmp" \
    "$RUNTIME_DIR/r-library" "$RUNTIME_DIR/config" "$RUNTIME_DIR/data" \
    "$RUNTIME_DIR/state"
export HOME="$RUNTIME_DIR/home"
export PIXI_HOME="$RUNTIME_DIR/pixi"
export PIXI_CACHE_DIR="$RUNTIME_DIR/cache/pixi"
export PIXI_NO_CONFIG=1 PIXI_NO_TELEMETRY=1
export XDG_CACHE_HOME="$RUNTIME_DIR/cache"
export XDG_CONFIG_HOME="$RUNTIME_DIR/config"
export XDG_DATA_HOME="$RUNTIME_DIR/data"
export XDG_STATE_HOME="$RUNTIME_DIR/state"
export TMPDIR="$RUNTIME_DIR/tmp"
export R_LIBS_USER="$RUNTIME_DIR/r-library"
export PYTHONNOUSERSITE=1
unset PYTHONPATH R_LIBS R_LIBS_SITE R_PROFILE_USER R_ENVIRON_USER

command=$1
shift
exec "$PIXI_BIN" "$command" --manifest-path "$PREPROCESSING_DIR/pixi.toml" "$@"
