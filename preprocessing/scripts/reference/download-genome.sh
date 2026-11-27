#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd -- "$SCRIPT_DIR/../../.." && pwd)
REFERENCE_DIR="$PROJECT_ROOT/index/reference_genome"
mkdir -p "$REFERENCE_DIR"
cd "$REFERENCE_DIR"

readonly URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
readonly RAW_FASTA="hg38.download.fa"
readonly OUTPUT_FASTA="hg38.fa"

# UCSC hg38 contains primary chromosomes plus alternative haplotypes and
# patches.  Keep only the chromosome names used by the GENCODE GTF so reads
# cannot be assigned equally to an unannotated alternative copy.
# GRCh38 chrY PAR1 and PAR2 coordinates are 1-based and inclusive; masking
# them leaves the PAR sequence on chrX as the single counting location.
curl --fail --location --retry 8 --retry-connrefused --connect-timeout 30 \
    --output hg38.fa.gz "$URL"
gzip -dc hg38.fa.gz | awk '
    function keep_chromosome(chromosome) {
        return chromosome ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/
    }

    function emit_sequence() {
        if (!keep || sequence == "") {
            return
        }
        if (chromosome == "chrY") {
            line_start = sequence_start
            line_end = sequence_start + length(sequence) - 1
            for (i = 1; i <= 2; i++) {
                mask_start = (i == 1 ? 10000 : 56887903)
                mask_end = (i == 1 ? 2781479 : 57217415)
                overlap_start = (line_start > mask_start ? line_start : mask_start)
                overlap_end = (line_end < mask_end ? line_end : mask_end)
                for (position = overlap_start; position <= overlap_end; position++) {
                    offset = position - line_start + 1
                    sequence = substr(sequence, 1, offset - 1) "N" substr(sequence, offset + 1)
                }
            }
        }
        print sequence
    }

    /^>/ {
        emit_sequence()
        chromosome = $1
        sub(/^>/, "", chromosome)
        keep = keep_chromosome(chromosome)
        sequence = ""
        sequence_start = 1
        if (keep) {
            print ">" chromosome
        }
        next
    }

    keep {
        sequence = sequence $0
        emit_sequence()
        sequence_start += length($0)
        sequence = ""
    }
' > "$RAW_FASTA"

mv "$RAW_FASTA" "$OUTPUT_FASTA"
rm -f hg38.fa.gz