#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/../common.sh"
shopt -s nullglob
BAMS=("$PROJECT_ROOT"/alignment/mapping/*/*.bam)

if (( ${#BAMS[@]} == 0 )); then
    echo "Error: no BAM files found in $PROJECT_ROOT/alignment/mapping." >&2
    exit 1
fi

echo "Checking ${#BAMS[@]} BAM files with samtools quickcheck..."
if "$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline samtools quickcheck -v "${BAMS[@]}"; then
    echo "All ${#BAMS[@]} BAM files passed quickcheck."
else
    echo "Error: BAM quickcheck failed; see the reported files above." >&2
    exit 1
fi
