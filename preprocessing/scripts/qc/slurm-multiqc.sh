#!/usr/bin/env bash
#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/multiqc-%j.log

set -euo pipefail

: "${PROJECT_ROOT:?Submit using the MultiQC launcher}"
: "${PIXI_BIN:?Submit using the MultiQC launcher}"
FASTQC_DIR="$PROJECT_ROOT/library_qc/fastqc"
MULTIQC_DIR="$PROJECT_ROOT/library_qc/multiqc"
MAPPING_DIR="$PROJECT_ROOT/alignment/mapping"
INDEX_DIR="$PROJECT_ROOT/index/star-index"
COUNT_DIR="$PROJECT_ROOT/quantification/Output"

for required in \
    "$INDEX_DIR/genomeParameters.txt" \
    "$INDEX_DIR/Log.out" \
    "$INDEX_DIR/chrNameLength.txt" \
    "$COUNT_DIR/featurecounts_summary_unique.tsv" \
    "$COUNT_DIR/featurecounts_summary_multiple.tsv"; do
    if [[ ! -s "$required" ]]; then
        echo "Error: integrated MultiQC requires completed index creation, alignment, and feature counting; missing $required" >&2
        exit 1
    fi
done

mapfile -t FASTQC_REPORTS < <(find "$FASTQC_DIR" -type f -name '*_fastqc.zip' -print)
mapfile -t STAR_LOGS < <(find "$MAPPING_DIR" -type f -name 'Log.final.out' -print)
mapfile -t STAR_GENE_COUNTS < <(find "$MAPPING_DIR" -type f -name 'ReadsPerGene.out.tab' -print)
if (( ${#FASTQC_REPORTS[@]} == 0 || ${#STAR_LOGS[@]} == 0 || ${#STAR_GENE_COUNTS[@]} == 0 )); then
    echo "Error: no FastQC, STAR alignment, or STAR gene-count reports found; check that upstream jobs completed." >&2
    exit 1
fi

INPUT_DIR=$(mktemp -d "${TMPDIR:-/tmp}/multiqc-inputs.XXXXXX")
trap 'rm -rf -- "$INPUT_DIR"' EXIT
"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline python \
    "$PROJECT_ROOT/preprocessing/scripts/qc/prepare-multiqc-inputs.py" \
    --project-root "$PROJECT_ROOT" --output-dir "$INPUT_DIR"

INPUT_LIST="$INPUT_DIR/multiqc_file_list.txt"
{
    printf '%s\n' "${FASTQC_REPORTS[@]}"
    printf '%s\n' "${STAR_LOGS[@]}"
    printf '%s\n' "${STAR_GENE_COUNTS[@]}"
    printf '%s\n' \
        "$INPUT_DIR/featurecounts_unique.summary" \
        "$INPUT_DIR/featurecounts_multiple.summary" \
        "$INPUT_DIR/featurecounts_modes_mqc.yaml" \
        "$INPUT_DIR/star_index_mqc.html"
} > "$INPUT_LIST"

"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline multiqc \
    --file-list "$INPUT_LIST" --outdir "$MULTIQC_DIR" \
    --filename fastq_qc_report --title "RNA-seq quality report" --force
