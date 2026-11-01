#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
shopt -s nullglob
BAMS=("$SCRIPT_DIR"/../mapping/*/*.bam)

if (( ${#BAMS[@]} == 0 )); then
    echo "Error: no BAM files found in $SCRIPT_DIR/../mapping." >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "Error: samtools is not available. Run this script in the sys-genomics-align environment." >&2
    exit 1
fi

echo "Checking ${#BAMS[@]} BAM files with samtools quickcheck..."
if samtools quickcheck -v "${BAMS[@]}"; then
    echo "All ${#BAMS[@]} BAM files passed quickcheck."
else
    echo "Error: BAM quickcheck failed; see the reported files above." >&2
    exit 1
fi
