############################################################
# GEO: GSE282742
# Analysis: Differential Gene Expression using DESeq2
# Author: Midhun
# Project: Machine Learning-Based Early Alzheimer’s Detection Using Blood Transcriptomics
# Threshold: pvalue < 0.05 (no log2FC cutoff)
# Date: 24-02-2026
############################################################

# -----------------------------
# 1. Install & Load Packages
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

library(DESeq2)
library(biomaRt)

set.seed(123)

# -----------------------------
# 2. Load Counts Matrix
# -----------------------------
cat("Select counts matrix (CSV: rows = genes, columns = samples)\n")
counts_mat <- read.csv(file.choose(), row.names = 1, check.names = FALSE)

# Ensure numeric
counts_mat <- as.matrix(counts_mat)
mode(counts_mat) <- "numeric"

# Round counts (if small decimals exist)
counts_mat <- round(counts_mat)
storage.mode(counts_mat) <- "integer"

cat("Counts matrix dimensions:", dim(counts_mat), "\n")

# -----------------------------
# 3. Load Metadata
# -----------------------------
cat("Select metadata file (must contain SampleID & Condition)\n")
meta <- read.csv(file.choose(), row.names = 1)

meta$Condition <- factor(meta$Condition)

# Check sample matching
if (!all(colnames(counts_mat) %in% rownames(meta))) {
  stop("Sample names in counts do not match metadata!")
}

meta <- meta[colnames(counts_mat), , drop = FALSE]

cat("Metadata loaded. Condition levels:\n")
print(levels(meta$Condition))

# -----------------------------
# 4. Create DESeq2 Dataset
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = meta,
  design    = ~ Condition
)

# Prefilter low counts
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# -----------------------------
# 5. Run Contrast (Modify if needed)
# -----------------------------
# Example: AD vs Control
# Replace level names based on your metadata

res <- results(dds)

res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)

cat("Total genes tested:", nrow(res_df), "\n")

# -----------------------------
# 6. Apply DEG Threshold
# -----------------------------
# Final Threshold:
# pvalue < 0.05
# No log2FC cutoff

deg_df <- res_df[
  !is.na(res_df$pvalue) &
    res_df$pvalue < 0.05,
]

deg_df$direction <- ifelse(deg_df$log2FoldChange > 0, "UP", "DOWN")

cat("Significant DEGs (pvalue < 0.05):", nrow(deg_df), "\n")
cat("UP genes:", sum(deg_df$direction == "UP"), "\n")
cat("DOWN genes:", sum(deg_df$direction == "DOWN"), "\n")

# -----------------------------
# 7. Annotate Gene Symbols
# -----------------------------
cat("Mapping Ensembl IDs to HGNC symbols...\n")

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = deg_df$ensembl_id,
  mart       = mart
)

deg_annot <- merge(
  deg_df,
  gene_map,
  by.x = "ensembl_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)

# -----------------------------
# 8. Save Results
# -----------------------------
write.csv(
  deg_annot,
  "GEO4_DEG_DESeq2_pvalue05.csv",
  row.names = FALSE
)

writeLines(
  deg_annot$ensembl_id[deg_annot$direction == "UP"],
  "GEO4_UP_ENS.txt"
)

writeLines(
  deg_annot$ensembl_id[deg_annot$direction == "DOWN"],
  "GEO4_DOWN_ENS.txt"
)

cat("\nAnalysis Complete ✅\n")
cat("Files generated:\n")
cat("- GEO4_DEG_DESeq2_pvalue05.csv\n")
cat("- GEO4_UP_ENS.txt\n")
cat("- GEO4_DOWN_ENS.txt\n")