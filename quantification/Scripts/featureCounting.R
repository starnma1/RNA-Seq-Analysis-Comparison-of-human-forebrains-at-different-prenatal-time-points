#!/usr/bin/env Rscript
# Count reads per gene from STAR BAMs and write gene x sample count matrices,
# once with unique reads only and once including multimapping reads.

library(Rsubread)
library(rtracklayer)

bam_dir  <- "alignment/mapping"
gtf_file <- "reference/Homo_sapiens.GRCh38.116.gtf"   # same GTF used for the STAR index
out_dir  <- "counts"
threads  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

dir.create(out_dir, showWarnings = FALSE)

# 1. Find the BAMs (STAR names them <sample>_Aligned.sortedByCoord.out.bam)
bams <- list.files(bam_dir, pattern = "\\.bam$", recursive = TRUE, full.names = TRUE)
sample_names <- sub("_?Aligned\\.sortedByCoord\\.out\\.bam$|\\.bam$", "", basename(bams))

# 2. Count: exons are grouped by gene_id, so each row is one gene
feature_count <- function(bam_files, annotation, multiple_reads) {
  featureCounts(
    files                  = bam_files,
    annot.ext              = annotation,
    isGTFAnnotationFile    = TRUE,
    GTF.featureType        = "exon",
    GTF.attrType           = "gene_id",
    useMetaFeatures        = TRUE,            # gene-level (meta-feature) counts
    strandSpecific         = 2,               # reverse-stranded; check "Assigned" in the stat file
    countMultiMappingReads = multiple_reads,  # FALSE = unique only, TRUE = include multimappers
    nthreads               = threads
  )
}

results <- list(
  unique   = feature_count(bams, gtf_file, multiple_reads = FALSE),
  multiple = feature_count(bams, gtf_file, multiple_reads = TRUE)
)

# Gene names and biotypes from the GTF (read once, used for both)
gtf   <- import(gtf_file)
genes <- as.data.frame(unique(mcols(gtf)[gtf$type == "gene",
                                         c("gene_id", "gene_name", "gene_biotype")]))

# 3-5. For each counting mode: build the matrix, annotate, save
for (mode in names(results)) {
  fc <- results[[mode]]
  
  # Matrix: genes = rows, samples = columns
  counts <- fc$counts
  colnames(counts) <- sample_names
  
  stat <- fc$stat
  colnames(stat) <- c("Status", sample_names)
  
  # Annotated table
  g <- genes[match(rownames(counts), genes$gene_id), ]
  annotated <- data.frame(gene_id      = rownames(counts),
                          gene_name    = g$gene_name,
                          gene_biotype = g$gene_biotype,
                          gene_length  = fc$annotation$Length,
                          counts,
                          check.names  = FALSE)
  
  # Save, with the mode in each file name
  write.csv(counts,    file.path(out_dir, paste0("gene_counts_matrix_", mode, ".csv")))
  write.csv(annotated, file.path(out_dir, paste0("gene_counts_annotated_", mode, ".csv")),
            row.names = FALSE)
  write.table(stat,    file.path(out_dir, paste0("featurecounts_summary_", mode, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(fc,          file.path(out_dir, paste0("featurecounts_result_", mode, ".rds")))
  
  # Quick check: fraction of reads assigned per sample
  assigned <- unlist(stat[stat$Status == "Assigned", -1]) / colSums(stat[, -1])
  cat("\nFraction of reads assigned (", mode, "):\n", sep = "")
  print(round(assigned, 3))
}