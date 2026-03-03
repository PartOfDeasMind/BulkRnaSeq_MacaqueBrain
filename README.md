# Bulk RNA-seq Analysis of Developing vs Adult Macaque Brain

## 📌 Project Overview

This project analyzes bulk RNA-seq data from macaque brain samples across developmental time points.  
The dataset includes prenatal and postnatal samples, allowing investigation of gene expression changes during brain development.

The goal was to compare two analysis strategies:

1. **Grouped Analysis** – Samples grouped into:
   - Prenatal
   - Postnatal

2. **Time-Resolved Analysis** – Each developmental age treated as a distinct class.

This project was completed as part of ETH Systems Genomics Fall 2025 course.

---

## 🧠 Biological Question

How does gene expression change during macaque brain development?

Specifically:
- What genes differ between prenatal and postnatal stages?
- How does the expression of those keys genes change across time?
- - What additional structure appears when modeling age as separate developmental classes?

---

## 🧪 Analysis Approaches

### 1️⃣ Grouped Analysis (Prenatal vs Postnatal)

Samples were collapsed into two biological conditions:
- Prenatal
- Postnatal

Differential expression analysis was performed between these two groups.

📂 Script:
`SysGen_Grouping_FINAL.R`

This approach increases statistical power but ignores finer temporal resolution.

---

### 2️⃣ Time-Resolved Analysis (Age as Distinct Classes)

Each age was treated as its own condition.

📂 Script:
`SysGen_NoGrouping_FINAL.R`

This approach captures more detailed developmental transitions but may reduce power due to smaller group sizes.

---

## ⚙️ Methods

Pipeline steps included:

- Quality control
- Normalization
- PCA and visualization
- Differential expression analysis
- K-means clustering (in time-resolved analysis only)
- Functional enrichment 

Main tools used:

Differential Expression
- **DESeq2** – normalization and differential expression testing  
- **tximport** – transcript-level import and gene-level summarization  

Data Manipulation
- **dplyr** – data wrangling and preprocessing  

Annotation
- **biomaRt** – gene annotation and identifier mapping  
- **org.Mmu.eg.db** – Macaca mulatta gene annotation database   

Visualization
- **ggplot2** – data visualization  
- **pheatmap** – heatmap generation  
- **viridis** – perceptually uniform color palettes  

Functional Enrichment
- **clusterProfiler** – Gene Ontology and pathway enrichment analysis  

---

## 📁 Repository Structure
```
.
├── SysGen_Grouping_FINAL.R        # Differential expression: Prenatal vs Postnatal grouping
├── SysGen_NoGrouping_FINAL.R      # Differential expression: Each age as distinct class
├── Presentation.pdf               # Final project presentation (exported from PowerPoint)
└── README.md                      # Project description and documentation
```

## 📂 Data Availability

The raw RNA-seq data are not included in this repository due to file size and licensing restrictions.

The dataset originates from:

Zhu Y, Sousa AMM, Gao T, Skarica M, Li M, Santpere G, Esteller-Cucala P, Juan D, Ferrández-Peral L, Gulden FO, Yang M, Miller DJ, Marques-Bonet T, Imamura Kawasawa Y, Zhao H, Sestan N. Spatiotemporal transcriptomic divergence across human and macaque brain development. Science. 2018 Dec 14;362(6420):eaat8077. doi: 10.1126/science.aat8077. Epub 2018 Dec 13. PMID: 30545855; PMCID: PMC6900982.

To reproduce the analysis:
1. Download the dataset from the original publication or GEO.
2. Update file paths in the scripts.
3. Run the R scripts.

## 👤 Author

Dea Karam
MSc Computational Biology & Bioinformatics  
ETH Zürich
