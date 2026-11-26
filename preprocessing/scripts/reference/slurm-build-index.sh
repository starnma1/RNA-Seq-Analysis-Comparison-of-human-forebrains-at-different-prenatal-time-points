#!/usr/bin/env bash
#SBATCH --job-name=star-index
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G              # 40 GB total; human genome needs ~32 GB
#SBATCH --time=04:00:00
#SBATCH --output=build-index-%j.log   # written in the folder you submit from

set -euo pipefail

: "${PROJECT_ROOT:?Submit using the preprocessing index launcher}"
: "${PIXI_BIN:?Submit using the preprocessing index launcher}"
cd "$PROJECT_ROOT/index"
mkdir -p star-index
rm -rf _STARtmp                       # STAR refuses to start if a leftover one exists
FASTA_FILE="./reference_genome/hg38.fa"
GTF_FILE="./reference_genome/gencode.v50.basic.annotation.gtf"
if [[ ! -s "$FASTA_FILE" ]]; then
    echo "Error: filtered hg38 FASTA not found at $PWD/$FASTA_FILE; run download-genome.sh first." >&2
    exit 1
fi
if [[ ! -s "$GTF_FILE" ]]; then
    echo "Error: GENCODE GTF not found at $PWD/$GTF_FILE; run download-gtf.sh first." >&2
    exit 1
fi

"$PROJECT_ROOT/preprocessing/pixi.sh" run --as-is --offline STAR \
     --runThreadN "$SLURM_CPUS_PER_TASK" \
     --runMode genomeGenerate \
     --genomeDir star-index \
     --genomeFastaFiles "$FASTA_FILE" \
     --sjdbGTFfile "$GTF_FILE" \
     --sjdbOverhang 100