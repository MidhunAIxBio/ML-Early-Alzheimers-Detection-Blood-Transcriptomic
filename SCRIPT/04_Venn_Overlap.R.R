############################################################
# Venn Diagram + Overlap Analysis for 3 GEO Datasets
############################################################

# Install if needed
if (!require(VennDiagram)) install.packages("VennDiagram")
library(VennDiagram)

# ===============================
# 1️⃣ Choose Gene List Files
# ===============================

cat("Select GSE249477 gene list file\n")
g1 <- scan(file.choose(), what="character")

cat("Select GSE270454 gene list file\n")
g2 <- scan(file.choose(), what="character")

cat("Select GSE282742 gene list file\n")
g3 <- scan(file.choose(), what="character")


# Remove duplicates
g1 <- unique(g1)
g2 <- unique(g2)
g3 <- unique(g3)

# ===============================
# 2️⃣ Compute Overlaps
# ===============================

common_all <- Reduce(intersect, list(g1, g2, g3))

cat("Total overlap across all three:", length(common_all), "\n")

# Pairwise overlaps
cat("GSE249477 ∩ GSE270454:", length(intersect(g1, g2)), "\n")
cat("GSE249477 ∩ GSE282742:", length(intersect(g1, g3)), "\n")
cat("GSE270454 ∩ GSE282742:", length(intersect(g2, g3)), "\n")

# ===============================
# 3️⃣ Save 73 Common Genes
# ===============================

writeLines(common_all, "Common_3GEO_Genes.txt")

# ===============================
# 4️⃣ Generate Venn Diagram
# ===============================

venn.diagram(
  x = list(
    GSE249477 = g1,
    GSE270454 = g2,
    GSE282742 = g3
  ),
  filename = "3GEO_Venn.png",
  fill = c("blue", "red", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Overlap of Significant Genes Across 3 GEO Datasets"
)

cat("Venn diagram saved as 3GEO_Venn.png\n")
