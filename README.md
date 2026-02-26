# ML-Based Early Alzheimer’s Detection Using Blood Transcriptomics
## 🧠 Project Workflow

![Workflow](docs/workflow.png)

## 📖 Project Overview
This project integrates multi-GEO bulk RNA-seq datasets to identify early biomarkers of Alzheimer’s disease using differential expression analysis, protein-protein interaction (PPI) network analysis, and machine learning.

## 📂 Datasets
- GSE249477
- GSE270454
- GSE282742

All datasets were obtained from the Gene Expression Omnibus (GEO).

---

## 🧬 Analysis Workflow

### 1️⃣ Differential Expression (DEG)
- Tool: DESeq2 (R)
- Threshold:
  - Adjusted p-value < 0.05
  - |log2FoldChange| ≥ 1

### 2️⃣ Overlap Analysis
- 3-GEO Venn analysis
- Common genes identified for downstream analysis

### 3️⃣ PPI Network Construction
- STRING database
- Confidence score ≥ 0.4
- Hub gene identification using Cytoscape (CytoHubba - MCC method)

### 4️⃣ Machine Learning
Machine learning models were applied to prioritize and evaluate candidate biomarkers:

- Support Vector Regression (SVR)
- Random Forest
- Linear Regression
- Ridge / Lasso Regression

Model evaluation included:
- R² Score
- Mean Squared Error (MSE)
- Residual analysis
- Predicted vs Actual comparison
- Feature importance ranking

---

## 🛠 Tools Used
- R (DESeq2)
- Python (Scikit-learn, Pandas, NumPy, Matplotlib, Seaborn)
- STRING
- Cytoscape (CytoHubba – MCC algorithm)

---

## 🎯 Objective
To identify robust blood-based transcriptomic biomarkers for early Alzheimer’s disease detection using integrated computational and machine learning approaches.

---

## 👨‍🔬 Author
Midhun AI x Bio
