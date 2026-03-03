############################################################
## 0. Setup: Working Directory & Libraries
############################################################

cat("=== 0. SETUP: Working directory and libraries ===\n")

setwd("C:/Users/deaka/OneDrive/Desktop/ETHZ/Sem1/SystemsGenomics")

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(biomaRt)
  library(ggplot2)
  library(pheatmap)
  library(tximport)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Mmu.eg.db)  
  library(org.Hs.eg.db)
  library(viridis)
})

cat("Libraries loaded successfully.\n\n")


############################################################
## 1. Load RSEM TPM Expression & Metadata
############################################################

cat("=== 1. LOADING TPM EXPRESSION AND META DATA ===\n")

# List all sample files
samples <- list.files("rsem_genes_results")
cat("Number of sample files found:", length(samples), "\n")

# Build expression matrix (TPM)
expr <- sapply(samples, function(sample){
  file <- paste0("rsem_genes_results/", sample)
  quant <- read.csv(file, sep = "\t", header = TRUE)
  tpm   <- setNames(quant$TPM, quant$gene_id)
  return(tpm)
})

# Clean column names: remove ".genes.results"
colnames(expr) <- sub("\\.genes\\.results$", "", colnames(expr))

# Load metadata
tmp1 <- read.csv("SraRunTable_new.csv", sep = ",", header = TRUE)
meta <- dplyr::select(tmp1, c("Run", "Age", "sex", "Sample_Name", "tissue"))

# Align columns of expr with meta$Run
expr <- expr[, meta$Run]

cat("Expression matrix dimensions (genes x samples): ", 
    paste(dim(expr), collapse = " x "), "\n")
cat("Meta data dimensions (samples x variables): ", 
    paste(dim(meta), collapse = " x "), "\n\n")


############################################################
## 2. Quick QC: Average Expression Histograms
############################################################

cat("=== 2. QC: Average expression distribution ===\n")

avg_expr <- rowMeans(expr)

layout(matrix(1:2, nrow = 1))
hist(avg_expr, main = "Mean TPM per gene", xlab = "TPM")
hist(log10(avg_expr + 1), main = "Log10(mean TPM + 1)", xlab = "log10(TPM + 1)")

cat("Plotted raw and log10-transformed mean TPM distributions.\n\n")


############################################################
## 3. Gene Annotation with biomaRt (Macaque Ensembl)
############################################################

cat("=== 3. GENE ANNOTATION WITH BIOMART ===\n")

ensembl <- useEnsembl(biomart = "genes", dataset = "mmulatta_gene_ensembl")

meta_genes_raw <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "gene_biotype"
  ),
  filters = "external_gene_name",
  values  = rownames(expr),   # gene symbols in expr
  mart    = ensembl
)

expr_ids <- data.frame(
  external_gene_name = rownames(expr),
  stringsAsFactors = FALSE
)

meta_genes <- expr_ids %>%
  left_join(meta_genes_raw, by = "external_gene_name") %>%
  distinct(external_gene_name, .keep_all = TRUE)

cat("Annotation table dimensions: ", 
    paste(dim(meta_genes), collapse = " x "), "\n")
cat("Rownames(expr) == meta_genes$external_gene_name ? ",
    all(rownames(expr) == meta_genes$external_gene_name), "\n\n")


############################################################
## 4. LOC Genes & Biotype Distribution
############################################################

cat("=== 4. LOC GENES & BIOTYPE DISTRIBUTION ===\n")

loc_genes <- grepl("^LOC", rownames(expr))
cat("Number of LOC* genes: ", sum(loc_genes), "\n")
cat("Number of non-LOC genes: ", sum(!loc_genes), "\n")

cat("Example of annotated genes:\n")
print(head(meta_genes))

cat("Unique gene biotypes:\n")
print(unique(meta_genes$gene_biotype))

pc_count <- sum(meta_genes$gene_biotype == "protein_coding", na.rm = TRUE)
cat("Number of protein_coding genes (with annotation): ", pc_count, "\n")

biotype_counts <- meta_genes %>%
  dplyr::count(gene_biotype) %>%       # <- force dplyr version
  dplyr::arrange(desc(n)) %>%
  dplyr::mutate(percent = 100 * n / sum(n))


cat("Biotype counts and percentages:\n")
print(biotype_counts)

ggplot(biotype_counts, aes(x = reorder(gene_biotype, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(n, " (", sprintf("%.1f", percent), "%)")), 
            vjust = -0.4, size = 3.5) +
  labs(
    x = "Gene Biotype",
    y = "Number of Genes",
    title = "Gene Biotype Distribution (Counts + Percentages)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cat("Plotted biotype distribution.\n\n")


############################################################
## 5. Filter to Protein-Coding & Expression Filters
############################################################

cat("=== 5. FILTERING TO PROTEIN-CODING AND EXPRESSED GENES ===\n")

keep_pc <- !is.na(meta_genes$gene_biotype) &
  meta_genes$gene_biotype == "protein_coding"

expr_pc <- expr[keep_pc, ]

cat("Protein-coding expression dimensions: ", 
    paste(dim(expr_pc), collapse = " x "), "\n")

avg_expr_pc <- rowMeans(expr_pc)
layout(matrix(1:2, nrow = 1))
hist(avg_expr_pc, main = "Mean TPM (protein-coding)", xlab = "TPM")
hist(log10(avg_expr_pc + 1), main = "Log10(mean TPM + 1)", xlab = "log10(TPM + 1)")

ggplot(data.frame(avg_expr = avg_expr_pc), aes(x = avg_expr)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(
    breaks = c(0, 1, 10, 100, 1000, 10000, 20000),
    trans = "log1p",
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(0, 1),
    expand = c(0, 0),
    trans = "log1p"
  ) +
  theme_minimal()

num_det <- rowSums(expr_pc > 0)
hist(num_det, main = "Number of samples with TPM > 0", xlab = "Samples detected")

expressed_pc <- rowMeans(expr_pc > 0) >= 0.5 | rowMeans(expr_pc) >= 1
expr_pc <- expr_pc[which(expressed_pc), ]

cat("Filtered expressed protein-coding matrix dimensions: ", 
    paste(dim(expr_pc), collapse = " x "), "\n")
cat("Number of genes kept after expression filter: ", nrow(expr_pc), "\n\n")


############################################################
## 6. Sample Correlations & Clustering
############################################################

cat("=== 6. SAMPLE CORRELATIONS & HIERARCHICAL CLUSTERING ===\n")

corr_pearson  <- cor(log1p(expr_pc))
corr_spearman <- cor(expr_pc, method = "spearman")

hcl_pearson   <- hclust(as.dist(1 - corr_pearson))
hcl_spearman  <- hclust(as.dist(1 - corr_spearman))

layout(matrix(1:2, nrow = 1))
plot(hcl_spearman, labels = meta$Age, main = "Spearman clustering by Age")
plot(hcl_pearson,  labels = meta$Age, main = "Pearson clustering by Age")

layout(matrix(1:2, nrow = 1))
plot(hcl_spearman, labels = meta$sex, main = "Spearman clustering by Sex")
plot(hcl_pearson,  labels = meta$sex, main = "Pearson clustering by Sex")

cat("Plotted clustering trees for Age and Sex.\n\n")


############################################################
## 7. PCA & Age Grouping
############################################################

cat("=== 7. PCA ANALYSIS & AGE GROUPING ===\n")

pca <- prcomp(log1p(t(expr_pc)), 
              center = TRUE, 
              scale.  = TRUE)

eigs <- pca$sdev^2
prop <- eigs / sum(eigs)
plot(1:length(prop), prop, xlab = "PC", ylab = "Proportion of variance",
     main = "PCA variance explained")

cat("Top 5 PCs variance explained (%):\n")
print(round(prop[1:5] * 100, 2))

# Order ages
age_levels <- c("E60","E81","E82","E110",
                "P0","P2",
                "7M",
                "1Y","2Y","3.5Y","4Y","5Y","7Y","11Y")

meta$Age <- factor(meta$Age, levels = age_levels)

# PCA by Age (earlier)
ggplot(data.frame(pca$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = Age), size = 5) +
  theme_bw() +
  ggtitle("PCA colored by age")

# PCA by sex (earlier)
ggplot(data.frame(pca$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = sex), size = 5) +
  theme_bw() +
  ggtitle("PCA colored by Sex")

# Define prenatal vs postnatal
meta$AgeGroup <- ifelse(grepl("^E", meta$Age), "Prenatal", "Postnatal")
meta$AgeGroup <- factor(meta$AgeGroup, levels = c("Prenatal", "Postnatal"))

cat("Counts per AgeGroup:\n")
print(table(meta$AgeGroup))

# PCA colored by AgeGroup
pca <- prcomp(log1p(t(expr_pc)), center = TRUE, scale. = TRUE)
var_expl <- pca$sdev^2 / sum(pca$sdev^2)

pc1_lab <- paste0("PC1 (", round(var_expl[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")

ggplot(data.frame(pca$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = AgeGroup), size = 5) +
  xlab(pc1_lab) +
  ylab(pc2_lab) +
  theme_bw() +
  ggtitle("PCA colored by Age Group")

cat("PCA plots generated (Sex and AgeGroup).\n\n")


############################################################
## 8. Correlation Heatmap with Age Annotation
############################################################

cat("=== 8. SAMPLE CORRELATION HEATMAP ===\n")

stopifnot(identical(colnames(expr_pc), colnames(corr_pearson)))

annotation_col <- data.frame(Age = meta$Age)
rownames(annotation_col) <- colnames(expr_pc)

age_levels <- c("E60","E81","E82","E110","P0","P2","7M",
                "1Y","2Y","3.5Y","4Y","5Y","7Y","11Y")
annotation_col$Age <- factor(annotation_col$Age, levels = age_levels)

age_colors <- setNames(
  viridis(length(age_levels)),
  age_levels
)

pheatmap(
  corr_pearson,
  annotation_col    = annotation_col,
  annotation_colors = list(Age = age_colors),
  color             = colorRampPalette(c("skyblue","yellow","red"))(100),
  border_color      = NA,
  main              = "Sample Correlation Matrix (Pearson)",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete"
)

cat("Correlation heatmap plotted.\n\n")


############################################################
## 9. ANOVA Across Ages & AgeGroup
############################################################

cat("=== 9. ANOVA ACROSS AGES AND PRENATAL/POSTNATAL ===\n")

# Across all ages
pvals_age <- apply(expr_pc, 1, function(g){
  summary(aov(g ~ meta$Age))[[1]][["Pr(>F)"]][1]
})
pvals_age_adj <- p.adjust(pvals_age, method = "BH")
sig_genes_age <- rownames(expr_pc)[pvals_age_adj < 0.05]

cat("Number of genes significant across all ages (BH FDR < 0.05): ",
    length(sig_genes_age), "\n")

# Prenatal vs Postnatal
pvals_group <- apply(expr_pc, 1, function(g){
  summary(aov(g ~ meta$AgeGroup))[[1]][["Pr(>F)"]][1]
})
pvals_group_adj <- p.adjust(pvals_group, method = "BH")
sig_genes_group <- rownames(expr_pc)[pvals_group_adj < 0.05]

cat("Number of genes significant for Prenatal vs Postnatal (BH FDR < 0.05): ",
    length(sig_genes_group), "\n\n")


############################################################
## 10. DESeq2 Differential Expression (AgeGroup factor)
############################################################

cat("=== 10. DESEQ2 DIFFERENTIAL EXPRESSION (Wald Test ~ AgeGroup) ===\n")

files <- file.path("rsem_genes_results", samples)
names(files) <- samples

txi <- tximport(files, type = "rsem", geneIdCol = "gene_id")
txi$length[txi$length == 0] <- 1

# Align meta to samples
samples_no_ext <- sub("\\.genes\\.results$", "", samples)
meta <- meta[match(samples_no_ext, meta$Run), ]

cat("After alignment, meta rows: ", nrow(meta), "\n")

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ AgeGroup)

# keep only genes that are in expr_pc as well
dds_filtered <- dds[
  intersect(rownames(expr_pc), rownames(dds)),
]

cat("DESeq2 filtered dataset dimensions: ",
    paste(dim(dds_filtered), collapse = " x "), "\n")

dds_filtered <- DESeq(
  dds_filtered
)

res_DESeq2 <- results(dds_filtered)
cat("DESeq2 Wald Test results summary:\n")
print(summary(res_DESeq2))
cat("\n")


############################################################
## 11. Volcano Plot (Postnatal vs Prenatal, from Wald Test)
############################################################

cat("=== 11. VOLCANO PLOT ===\n")

res_df <- as.data.frame(res_DESeq2)
res_df <- res_df[!is.na(res_df$padj), ]


sum(res_df$padj < 0.05)

# Define significance groups: here using |log2FC| > 2 as threshold
res_df$group <- "NS"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange >  2] <- "Up_Postnatal"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < -2] <- "Up_Prenatal"

res_df$group <- factor(res_df$group,
                       levels = c("NS", "Up_Postnatal", "Up_Prenatal"))

cat("Counts per group in volcano plot:\n")
print(table(res_df$group))

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("grey70", "red3", "dodgerblue3"),
                     name   = "",
                     labels = c("NS",
                                "Up_Postnatal",
                                "Up_Prenatal")) +
  labs(
    x = "log2 fold change",
    y = "-log10(FDR)"
  ) +
  theme_minimal()

cat("Volcano plot generated.\n\n")


############################################################
## 12. Heatmap of Top 50 DE Genes (VST-scaled)
############################################################

cat("=== 12. HEATMAP OF TOP 100 DE GENES ===\n")

vsd <- vst(dds_filtered, blind = FALSE)
mat <- assay(vsd)

# Top 50 genes by adjusted p-value
top_genes <- rownames(res_df)[order(res_df$padj)][1:1000]
cat("Top 100 genes (by padj):\n")
print(head(top_genes))

mat_top <- mat[top_genes, ]
mat_top_scaled <- t(scale(t(mat_top)))   # row-wise Z-score

annotation_col <- data.frame(Age = meta$Age)
rownames(annotation_col) <- colnames(mat_top_scaled)

age_colors <- setNames(
  viridis(length(age_levels)),
  age_levels
)


pheatmap(
  mat_top_scaled,
  annotation_col = annotation_col,
  annotation_colors = list(Age = age_colors),
  show_rownames  = FALSE,
  fontsize_col   = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  main = "Top 1000 Significant Genes"
)

cat("Top 50 gene heatmap generated.\n\n")


############################################################
## 14. GO Enrichment 
############################################################


res_df$gene <- rownames(res_df)
up_postnatal <- res_df$gene[res_df$group == "Up_Postnatal"]
up_prenatal  <- res_df$gene[res_df$group == "Up_Prenatal"]

universe <- rownames(expr_pc)

length(up_postnatal)
length(up_prenatal)

ego <- enrichGO(
  gene          = up_postnatal,
  universe      = universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)


dotplot(ego, showCategory = 15) +
  ggtitle(paste0("GO BP enrichment: Postnatal "))


ego <- enrichGO(
  gene          = up_prenatal,
  universe      = universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)


dotplot(ego, showCategory = 15) +
  ggtitle(paste0("GO BP enrichment: Prenatal "))
