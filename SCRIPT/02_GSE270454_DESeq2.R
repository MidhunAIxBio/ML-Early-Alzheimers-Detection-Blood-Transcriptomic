############################################################
# GEO : GSE27045
# Analysis: Differential Gene Expression using DESeq2
# Final threshold: pvalue < 0.1 & |log2FC| > 0.25
# Author: Midhun
# Project: Machine Learning-Based Identification of Blood Transcriptomic Biomarkers for Early Alzheimer’s Disease
# Date: 24-02-2026
############################################################

# ===============================
# 1. Install / Load Packages
# ===============================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("DESeq2", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE)
  library(p, character.only = TRUE)
}

# Set Seed (Reproducibility)
set.seed(123)

# ===============================
# 2. User Parameters
# ===============================
min_count_sum <- 10
p_cut  <- 0.1
lfc_cut <- 0.25

counts_path <- "GSE270454_RNAseq-combined-counts-matrix.csv"
meta_path   <- "metadata_filled.csv"

outdir <- "outputs_geo2"
if (!dir.exists(outdir)) dir.create(outdir)

# ===============================
# 3. Load Data
# ===============================
counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_path, row.names = 1, check.names = FALSE)

# Clean names
clean_names <- function(x) trimws(gsub('"', '', x))
colnames(counts) <- clean_names(colnames(counts))
rownames(meta)   <- clean_names(rownames(meta))

# Align samples
meta <- meta[colnames(counts), ]
stopifnot(all(colnames(counts) == rownames(meta)))

# Merge ASM/ASO → Control
meta$condition <- as.character(meta$condition)
meta$condition[meta$condition %in% c("ASM","ASO")] <- "Control"
meta$condition <- factor(meta$condition)
meta$condition <- relevel(meta$condition, ref = "Control")

# ===============================
# 4. Filter Low Counts
# ===============================
counts_f <- counts[rowSums(counts) >= min_count_sum, ]

# ===============================
# 5. Run DESeq2
# ===============================
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_f)),
  colData = meta,
  design = ~ condition
)

dds <- DESeq(dds)

# ===============================
# 6. Define Contrasts
# ===============================
contrasts <- list(
  AD_vs_Control  = c("condition", "AD", "Control"),
  MCI_vs_AD      = c("condition", "MCI", "AD"),
  MCI_vs_Control = c("condition", "MCI", "Control")
)

combined_up   <- list()
combined_down <- list()

# ===============================
# 7. Extract Results
# ===============================
for (name in names(contrasts)) {
  
  message("Processing: ", name)
  
  res <- results(dds, contrast = contrasts[[name]])
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Save FULL DEG table (for ML)
  full_out <- file.path(outdir, paste0(name, "_FULL_DEG_for_ML.csv"))
  write.csv(res_df, full_out, row.names = FALSE, quote = FALSE)
  
  # Apply FINAL threshold
  df_relaxed <- res_df %>%
    filter(!is.na(pvalue) & pvalue < p_cut &
             !is.na(log2FoldChange) &
             abs(log2FoldChange) > lfc_cut)
  
  # Save full relaxed table
  relaxed_out <- file.path(outdir,
                           paste0(name, "_relaxed_p01_lfc025_full.csv"))
  write.csv(df_relaxed, relaxed_out, row.names = FALSE, quote = FALSE)
  
  # Split UP and DOWN
  up_genes <- df_relaxed %>%
    filter(log2FoldChange > lfc_cut) %>%
    pull(gene) %>% unique()
  
  down_genes <- df_relaxed %>%
    filter(log2FoldChange < -lfc_cut) %>%
    pull(gene) %>% unique()
  
  write.csv(data.frame(Gene = up_genes),
            file.path(outdir,
                      paste0(name, "_relaxed_UP_p01_lfc025.csv")),
            row.names = FALSE, quote = FALSE)
  
  write.csv(data.frame(Gene = down_genes),
            file.path(outdir,
                      paste0(name, "_relaxed_DOWN_p01_lfc025.csv")),
            row.names = FALSE, quote = FALSE)
  
  combined_up[[name]]   <- up_genes
  combined_down[[name]] <- down_genes
  
  message("  -> total: ", nrow(df_relaxed),
          " | UP: ", length(up_genes),
          " | DOWN: ", length(down_genes))
}

# ===============================
# 8. Create Master Combined Lists
# ===============================
all_up   <- sort(unique(unlist(combined_up)))
all_down <- sort(unique(unlist(combined_down)))

write.csv(data.frame(Gene = all_up),
          file.path(outdir, "ALL_combined_UP_p01_lfc025.csv"),
          row.names = FALSE, quote = FALSE)

write.csv(data.frame(Gene = all_down),
          file.path(outdir, "ALL_combined_DOWN_p01_lfc025.csv"),
          row.names = FALSE, quote = FALSE)

message("✅ GEO 2 analysis complete.")
message("Results saved in folder: ", outdir)