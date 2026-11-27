#!/usr/bin/env bash
#SBATCH --job-name=package-outputs
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G # 64 GB (zstd --ultra -22 needs ~1-2 GB per thread)
#SBATCH --time=24:00:00
#SBATCH --output=logs/package-%j.log

set -euo pipefail

: "${PROJECT_ROOT:?Submit using the preprocessing packaging launcher}"
: "${PIXI_BIN:?Submit using the preprocessing packaging launcher}"
cd "$PROJECT_ROOT/alignment"
SRC=mapping
ARCHIVE=mapping.tar.zst
TMP="$ARCHIVE.partial"

# --- compress at the maximum zstd level ---------------------------------
# --ultra -22 is zstd's highest level; its 128 MB window is still within the
# default decoder limit, so the archive extracts with a plain `tar --zstd -xf`.
echo "Packaging $SRC -> $ARCHIVE ($(du -sh "$SRC" | cut -f1))"
"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline tar --exclude='_STARtmp' \
    -I "zstd --ultra -22 -T$SLURM_CPUS_PER_TASK" \
    -cvf "$TMP" "$SRC"

"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline zstd -t "$TMP"
mv "$TMP" "$ARCHIVE"
echo "Done: $(du -sh "$ARCHIVE" | cut -f1) $ARCHIVE"
