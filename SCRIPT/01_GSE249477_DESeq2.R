############################################################
# GEO: GSE249477
# Analysis: Differential Gene Expression using DESeq2
# Final threshold: pvalue <= 0.1 |log2FoldChange| >= 0.25
# Author: Midhun
# Project: Machine Learning-Based Early Alzheimer’s Detection Using Blood Transcriptomics
# Date: 24-02-2026
############################################################

# =========================
# 1. Load Libraries
# =========================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

library(DESeq2)

#Set Seed (Reproducibility)
set.seed(123)


# =========================
# 2. Load Input Files
# =========================
# Replace with your file paths
raw_counts  <- read.csv("raw_counts.csv", header = TRUE)
metadata    <- read.csv("metadata.csv", header = TRUE)

# =========================
# 3. Extract Count Matrix
# =========================

gene_names <- make.unique(raw_counts$Name)

count_cols <- grep("Total.counts", colnames(raw_counts), value = TRUE)

counts <- raw_counts[, count_cols]
rownames(counts) <- gene_names

counts <- as.matrix(counts)
mode(counts) <- "integer"

# =========================
# 4. Prepare Metadata
# =========================

metadata$Condition <- factor(metadata$Condition,
                             levels = c("Control", "MCI", "AD"))

rownames(metadata) <- metadata$SampleID

# Ensure alignment
counts <- counts[, metadata$SampleID]

# =========================
# 5. DESeq2 Pipeline
# =========================

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ Condition
)

# Remove very low count genes
dds <- dds[rowSums(counts(dds)) > 5, ]

dds <- DESeq(dds)

# =========================
# 6. Extract Results
# =========================

res_ADvsCtrl  <- results(dds, contrast = c("Condition", "AD", "Control"))
res_MCIvsCtrl <- results(dds, contrast = c("Condition", "MCI", "Control"))

res_ADvsCtrl_df  <- as.data.frame(res_ADvsCtrl)
res_MCIvsCtrl_df <- as.data.frame(res_MCIvsCtrl)

res_ADvsCtrl_df$gene  <- rownames(res_ADvsCtrl_df)
res_MCIvsCtrl_df$gene <- rownames(res_MCIvsCtrl_df)

# =========================
# 7. Apply FINAL Threshold
# =========================

p_cut  <- 0.1
lfc_cut <- 0.25

extract_deg <- function(res_df, name) {
  
  sig <- res_df[
    !is.na(res_df$pvalue) &
      res_df$pvalue <= p_cut &
      abs(res_df$log2FoldChange) >= lfc_cut,
  ]
  
  up   <- sig[sig$log2FoldChange > 0, ]
  down <- sig[sig$log2FoldChange < 0, ]
  
  cat(name, ": total =", nrow(sig),
      "| UP =", nrow(up),
      "| DOWN =", nrow(down), "\n")
  
  return(list(sig = sig, up = up, down = down))
}

deg_AD  <- extract_deg(res_ADvsCtrl_df,  "AD_vs_Control")
deg_MCI <- extract_deg(res_MCIvsCtrl_df, "MCI_vs_Control")

# =========================
# 8. Save DEG Outputs
# =========================

outdir <- "GEO1_DEG_outputs"
dir.create(outdir, showWarnings = FALSE)

write.csv(deg_AD$sig,
          file.path(outdir, "AD_vs_Control_DEG_full.csv"),
          row.names = FALSE)

write.csv(deg_MCI$sig,
          file.path(outdir, "MCI_vs_Control_DEG_full.csv"),
          row.names = FALSE)

# Gene-only lists
write.csv(data.frame(Gene = deg_AD$up$gene),
          file.path(outdir, "AD_UP_genes.csv"),
          row.names = FALSE)

write.csv(data.frame(Gene = deg_AD$down$gene),
          file.path(outdir, "AD_DOWN_genes.csv"),
          row.names = FALSE)

write.csv(data.frame(Gene = deg_MCI$up$gene),
          file.path(outdir, "MCI_UP_genes.csv"),
          row.names = FALSE)

write.csv(data.frame(Gene = deg_MCI$down$gene),
          file.path(outdir, "MCI_DOWN_genes.csv"),
          row.names = FALSE)

# =========================
# 9. Create FULL DEG File for ML
# =========================

make_full_deg_ml <- function(res_df, comparison_name) {
  
  df <- res_df[, c(
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "gene"
  )]
  
  df$comparison <- comparison_name
  df$direction  <- ifelse(df$log2FoldChange > 0, "UP", "DOWN")
  
  return(df)
}

full_AD  <- make_full_deg_ml(res_ADvsCtrl_df,  "AD_vs_Control")
full_MCI <- make_full_deg_ml(res_MCIvsCtrl_df, "MCI_vs_Control")

full_all <- rbind(full_AD, full_MCI)

write.csv(full_all,
          file.path(outdir, "GEO1_ALL_FULL_DEG_for_ML.csv"),
          row.names = FALSE)

cat("\n✔ GEO1 pipeline completed successfully\n")
cat("Outputs saved in:", normalizePath(outdir), "\n")