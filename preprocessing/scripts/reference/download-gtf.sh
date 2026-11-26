#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd -- "$SCRIPT_DIR/../../.." && pwd)
mkdir -p "$PROJECT_ROOT/index/reference_genome"
cd "$PROJECT_ROOT/index/reference_genome"
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz
gzip -d gencode.v50.basic.annotation.gtf.gz