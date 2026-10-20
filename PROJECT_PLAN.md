# Human Prenatal Forebrain Development (E-MTAB-6814)

**Question:** Comparison of human forebrains at different prenatal time points (4–20 WPC)
**Data:** Cardoso-Moreira *et al.* 2019, *Nature* — ArrayExpress E-MTAB-6814

Mark done by changing ⬜ to ✅. Replace the target date with the date it was actually done.

| Done | Date | Phase | Step | Tools | Output |
|:---:|---|---|---|---|---|
| ⬜ | 23 Oct | 1. Download & organize Euler | Euler access + conda environment for both | Euler, conda | Working environment, `environment.yml` |
| ⬜ | 23 Oct | 1. Download & organize Euler | Set up folder structure (home + scratch) and GitHub repo | Euler, git | Project folders, repo |
| ⬜ | 23 Oct | 1. Download & organize Euler | Download sample sheet (SDRF) | BioStudies | `E-MTAB-6814.sdrf.txt` |
| ⬜ | 30 Oct | 1. Download & organize Euler | Select prenatal forebrain samples | awk / R | `samples_forebrain_prenatal.tsv`, `run_list.txt` |
| ⬜ | 30 Oct | 1. Download & organize Euler | Download FASTQs | ENA, wget | Raw `.fastq.gz` files on scratch |
| ⬜ | 30 Oct | 1. Download & organize Euler | Check downloads are complete and intact | gzip, md5sum | Download check log |
| ⬜ | 30 Oct | 2. QC & alignment | Raw read quality control | FastQC, MultiQC | MultiQC report (raw reads) |
| ⬜ | 30 Oct | 2. QC & alignment | Decide on trimming | fastp (if needed) | Decision + trimmed FASTQs (if needed) |
| ⬜ | 6 Nov | 2. QC & alignment | Download reference genome + annotation | GENCODE | Genome FASTA, GTF |
| ⬜ | 6 Nov | 2. QC & alignment | Build genome index | STAR | STAR index |
| ⬜ | 6 Nov | 2. QC & alignment | Align all samples | STAR | Sorted BAM files |
| ⬜ | 6 Nov | 2. QC & alignment | Alignment QC (mapping rates, strandedness) | STAR logs, MultiQC | Mapping stats table |
| ⬜ | 13 Nov | 3. Assignment & functional characterization | Assign reads to genes | featureCounts | Gene count table |
| ⬜ | 13 Nov | 3. Assignment & functional characterization | Assignment QC (% assigned per sample) | featureCounts summary, MultiQC | Assignment stats table |
| ⬜ | 13 Nov | 3. Assignment & functional characterization | Library composition (gene types, mito, rRNA) | featureCounts | Composition plot |
| ⬜ | 13 Nov | 3. Assignment & functional characterization | Sex check | R | Sex check plot |
| ⬜ | 13 Nov | 3. Assignment & functional characterization | Gene annotation (symbols, gene types, markers) | GENCODE GTF, biomaRt | Gene annotation table |
| ⬜ | 13 Nov | 3. Assignment & functional characterization | Copy results off scratch | rsync | Counts + QC saved in home |
| ⬜ | 20 Nov | 4. Final data matrices | Build count matrix + sample metadata | R | `counts_raw`, `sample_metadata` |
| ⬜ | 20 Nov | 4. Final data matrices | Remove failed / outlier samples | R | List of excluded samples |
| ⬜ | 20 Nov | 4. Final data matrices | Filter lowly expressed genes | edgeR | Filtered count matrix |
| ⬜ | 20 Nov | 4. Final data matrices | Normalize | edgeR, DESeq2 | Normalized matrices (logCPM, VST) |
| ⬜ | 20 Nov | 4. Final data matrices | Exploration: sample correlation, clustering, PCA | R, DESeq2 | Correlation heatmap, PCA plot |
| ⬜ | 20 Nov | 4. Final data matrices | Marker gene trajectories over time | R | Marker plots |
| ⬜ | 27 Nov | 5. Differential expression | DE across time points | DESeq2 | DESeq2 results table |
| ⬜ | 27 Nov | 5. Differential expression | DE across time points | edgeR | edgeR results table |
| ⬜ | 27 Nov | 5. Differential expression | Compare DESeq2 vs edgeR | R | Overlap plot, final DEG list |
| ⬜ | 27 Nov | 5. Differential expression | Adjacent time point comparisons | DESeq2 / edgeR | DEG counts per transition, plot |
| ⬜ | 4 Dec | 5. Differential expression | Group DEGs by temporal pattern | R (clustering) | Cluster heatmap, cluster profiles |
| ⬜ | 4 Dec | 5. Differential expression | Functional enrichment of DEG clusters | clusterProfiler / DAVID | GO enrichment tables + plots |
| ⬜ | 4 Dec | 5. Differential expression | Compare with the paper's results | R | Comparison table |
| ⬜ | 11 Dec | Wrap-up | Final figures + write-up | — | Report / presentation |
