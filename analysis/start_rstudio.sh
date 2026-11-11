#!/usr/bin/env bash
# Open this project's RStudio Project using the R from the forebrain-rnaseq conda env.
# Lives in <project>/analysis/start_rstudio.sh and works from any directory:
#
#   ./analysis/start_rstudio.sh            # from the project root
#   ./start_rstudio.sh                     # from inside analysis/
#   ./analysis/start_rstudio.sh other-env  # use a different conda env

set -euo pipefail

ENV_NAME="${1:-forebrain-rnaseq}"
RSTUDIO_BIN="/Applications/RStudio.app/Contents/MacOS/RStudio"

# Project root = the folder above analysis/ (where this script lives)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Find the .Rproj file in the project root
shopt -s nullglob
RPROJ_FILES=("$PROJECT_ROOT"/*.Rproj)
shopt -u nullglob
if [[ ${#RPROJ_FILES[@]} -eq 0 ]]; then
  echo "Error: no .Rproj file found in $PROJECT_ROOT" >&2
  exit 1
fi
RPROJ="${RPROJ_FILES[0]}"

# Make "conda activate" work inside a script
CONDA_BASE="$(conda info --base 2>/dev/null)" || {
  echo "Error: conda not found. Is Anaconda on your PATH?" >&2
  exit 1
}
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"

# conda's activation scripts use unset variables, so relax "set -u" just for this step
set +u
conda activate "$ENV_NAME" || {
  echo "Error: could not activate conda environment '$ENV_NAME'." >&2
  echo "Check it exists with: conda info --envs" >&2
  echo "Create it with: conda env create -f \"$PROJECT_ROOT/environments/analysis.yml\"" >&2
  exit 1
}
set -u

R_PATH="$CONDA_PREFIX/bin/R"
[[ -x "$R_PATH" ]] || { echo "Error: no R in env '$ENV_NAME' ($R_PATH)." >&2; exit 1; }
[[ -x "$RSTUDIO_BIN" ]] || { echo "Error: RStudio not found at $RSTUDIO_BIN" >&2; exit 1; }

export RSTUDIO_WHICH_R="$R_PATH"
cd "$PROJECT_ROOT"

echo "R:       $R_PATH"
echo "Project: $RPROJ"

# Start RStudio detached from this terminal
nohup "$RSTUDIO_BIN" "$RPROJ" >/dev/null 2>&1 &
disown
