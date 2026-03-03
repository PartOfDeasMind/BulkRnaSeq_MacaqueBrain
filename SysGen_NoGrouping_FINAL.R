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
meta <- dplyr::select(tmp1, c("Run", "Age", "sex", "Sample_Name", "tissue", "Instrument" ))

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
  geom_bar(stat = "identity", fill = "#215caf") +
  geom_text(aes(label = paste0(n, " (", sprintf("%.1f", percent), "%)")), 
            vjust = -0.4, size = 3.5) +
  labs(
    x = "Gene Biotype",
    y = "Number of Genes",
    title = "Gene Biotype Distribution (Counts + Percentages)"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

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


pc1_var <- round(prop[1] * 100, 1)
pc2_var <- round(prop[2] * 100, 1)

# PCA by Age (earlier)
ggplot(data.frame(pca$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = Age), size = 5) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    title = "PCA by Age"
  )


# PCA by sex (earlier)
ggplot(data.frame(pca$x, meta)) +
  geom_point(aes(x = PC1, y = PC2, color = sex), size = 5) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    title = "PCA by Sex"
  )


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
  color             = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  border_color      = NA,
  main              = "Sample Correlation Matrix",
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
## 10. DESeq2 Differential Expression (Age factor)
############################################################

cat("=== 10. DESEQ2 DIFFERENTIAL EXPRESSION (LRT ~ Age) ===\n")

files <- file.path("rsem_genes_results", samples)
names(files) <- samples

txi <- tximport(files, type = "rsem", geneIdCol = "gene_id")
txi$length[txi$length == 0] <- 1

# Align meta to samples
samples_no_ext <- sub("\\.genes\\.results$", "", samples)
meta <- meta[match(samples_no_ext, meta$Run), ]

cat("After alignment, meta rows: ", nrow(meta), "\n")

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ Age)

# keep only genes that are in expr_pc as well
dds_filtered <- dds[
  intersect(rownames(expr_pc), rownames(dds)),
]

cat("DESeq2 filtered dataset dimensions: ",
    paste(dim(dds_filtered), collapse = " x "), "\n")

dds_filtered <- DESeq(
  dds_filtered,
  test    = "LRT",
  reduced = ~ 1
)

res_DESeq2 <- results(dds_filtered)
cat("DESeq2 LRT results summary:\n")
print(summary(res_DESeq2))
cat("\n")


############################################################
## 11. Volcano Plot LRT)
############################################################

cat("=== 11. VOLCANO PLOT ===\n")

res_df <- as.data.frame(res_DESeq2)
res_df <- res_df[!is.na(res_df$padj), ]

sum(res_df$padj < 0.01)

# Define significance groups: here using |log2FC| > 2 as threshold
res_df$group <- "NS"
res_df$group[res_df$padj < 0.01 & res_df$log2FoldChange >  2] <- "Up"
res_df$group[res_df$padj < 0.01 & res_df$log2FoldChange < -2] <- "Down"

res_df$group <- factor(res_df$group,
                       levels = c("NS", "Up", "Down"))

cat("Counts per group in volcano plot:\n")
print(table(res_df$group))

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  scale_color_manual(values = c("grey70", "red3", "dodgerblue3"),
                     name   = "",
                     labels = c("NS",
                                "Up",
                                "Down")) +
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
  colorRampPalette(c("#440154", "#31688E", "#35B779", "#FDE725"))(length(age_levels)),
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
  main = "Top 50 DE Genes (VST, Z-score per gene)"
)

cat("Top 50 gene heatmap generated.\n\n")















############################################################
## HEATMAP from expr_pc (TPM), log-transform + Z-score
## Select top 1000 significant genes from res_df
############################################################

cat("=== 12. HEATMAP OF TOP 1000 SIGNIFICANT GENES (expr_pc TPM) ===\n")

# 1) Get top 1000 significant genes from res_df (DESeq2 LRT)
sig_res <- res_df %>%
  filter(!is.na(padj), padj < 0.01) %>%     # choose threshold (0.05 or 0.01)
  arrange(padj)

topN <- 100
top_genes <- rownames(sig_res)[1:min(topN, nrow(sig_res))]

cat("Top genes selected: ", length(top_genes), "\n")

# 2) Keep only genes that exist in expr_pc
top_genes <- intersect(top_genes, rownames(expr_pc))
cat("Top genes present in expr_pc: ", length(top_genes), "\n")

# 3) Subset expr_pc and transform
mat_top <- expr_pc[top_genes, , drop = FALSE]     # genes x samples

# log-transform TPM
mat_log <- log2(mat_top + 1)                      # or log1p(mat_top)

# 4) Row-wise Z-score (per gene)
mat_scaled <- t(scale(t(mat_log)))

# Optional: handle genes with zero variance (can become NA after scaling)
mat_scaled <- mat_scaled[complete.cases(mat_scaled), , drop = FALSE]

cat("Matrix for heatmap: ", paste(dim(mat_scaled), collapse = " x "), "\n")

# 5) Annotation
annotation_col <- data.frame(Age = meta$Age)
rownames(annotation_col) <- colnames(mat_scaled)

# Make sure Age is a factor with the correct order
annotation_col$Age <- factor(annotation_col$Age, levels = age_levels)

age_colors <- setNames(
  viridis(length(age_levels)),
  age_levels
)

# 6) Plot heatmap
pheatmap(
  mat_scaled,
  annotation_col    = annotation_col,
  annotation_colors = list(Age = age_colors),
  show_rownames     = FALSE,
  fontsize_col      = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  border_color      = NA,
  main              = paste0("Top ", nrow(mat_scaled), " significant genes")
)

cat("Heatmap generated from expr_pc.\n\n")





















































############################################################
## 13. K-means Clustering on Top 1000 DE Genes
############################################################

cat("=== 13. K-MEANS CLUSTERING ON TOP 1000 DE GENES ===\n")

# Filter for significant genes
sig_genes_km <- res_df[!is.na(res_df$padj) & res_df$padj < 0.01, ]
cat("Number of genes with padj < 0.01: ", nrow(sig_genes_km), "\n")

sig_genes_km <- sig_genes_km[order(sig_genes_km$padj), ]
top1000_genes <- rownames(sig_genes_km)[1:min(1000, nrow(sig_genes_km))]

cat("Number of genes taken for k-means (up to 1000): ", length(top1000_genes), "\n")

# Subset genes
expr_sub <- expr_pc[top1000_genes, , drop = FALSE]   # genes x samples

# 1) Log-transform TPM
#expr_log <- log2(expr_sub + 1)    # or use log1p(expr_sub)

# 2) Z-score per gene (row-wise)
expr_scaled <- t(scale(t(expr_sub)))

set.seed(123)
k <- 5
km <- kmeans(expr_scaled, centers = k, nstart = 50, iter.max = 10)
clusters <- km$cluster

cat("Cluster sizes (genes per cluster):\n")
print(table(clusters))

# Age factor
age <- factor(meta$Age)
age_levels <- levels(age)

cluster_genes <- vector("list", k)

par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))

for (cl in 1:k) {
  genes_cl <- names(clusters[clusters == cl])
  cluster_genes[[cl]] <- genes_cl
  
  expr_cl <- expr_scaled[genes_cl, , drop = FALSE]
  
  box_dat <- lapply(age_levels, function(a) {
    as.vector(expr_cl[, age == a, drop = FALSE])
  })
  
  boxplot(
    box_dat,
    names = age_levels,
    main  = sprintf("Cluster %d (%d genes)", cl, length(genes_cl)),
    xlab  = "Age",
    ylab  = "Scaled expression"
  )
}

cat("K-means clustering and cluster expression boxplots done.\n\n")


############################################################
## 14. GO Enrichment for Each Cluster
############################################################

cat("=== 14. GO ENRICHMENT PER CLUSTER ===\n")

# Define universe for GO analysis
universe <- rownames(expr_pc)
cat("GO universe size (protein-coding expressed genes): ", length(universe), "\n")

ego_list <- vector("list", k)

for (cl in 1:k) {
  cat(sprintf("Running GO enrichment for cluster %d (%d genes)...\n",
              cl, length(cluster_genes[[cl]])))
  
  ego <- enrichGO(
    gene          = cluster_genes[[cl]],
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  
  ego_list[[cl]] <- ego
}


ego <- ego_list[[1]]
dotplot(ego, showCategory = 15) +
  ggtitle(paste0("GO BP enrichment: Cluster ", 1))

ego <- ego_list[[2]]
dotplot(ego, showCategory = 15) +
  ggtitle(paste0("GO BP enrichment: Cluster ", 2))

ego <- ego_list[[3]]
dotplot(ego, showCategory = 15) +
  ggtitle(paste0("GO BP enrichment: Cluster ", 3))

ego <- ego_list[[4]]
dotplot(ego, showCategory = 15) +
  ggtitle(paste0("GO BP enrichment: Cluster ", 4))

ego <- ego_list[[5]]
dotplot(ego, showCategory = 15) +
  ggtitle(paste0("GO BP enrichment: Cluster ", 5))

cat("GO enrichment analysis finished for all clusters.\n\n")


############################################################
## 15. Extra: Gene Sets by Peak Age
############################################################

cat("=== 15. EXTRA: GENE SUBSETS BY PEAK AGE ===\n")

DEG <- rownames(res_df)[
  !is.na(res_df$padj) &
    res_df$padj < 0.01
]

cat("Number of DE genes used for peak-age analysis: ", length(DEG), "\n")

avg_expr_age <- sapply(age_levels, function(age) {
  rowMeans(expr_pc[, meta$Age == age, drop = FALSE])
})

avg_expr_DEG <- avg_expr_age[DEG, , drop = FALSE]

max_age_DEG <- setNames(
  colnames(avg_expr_DEG)[apply(avg_expr_DEG, 1, which.max)],
  rownames(avg_expr_DEG)
)

cat("Number of genes peaking at each age:\n")
print(table(max_age_DEG))

avg_expr_DEG_list <- tapply(names(max_age_DEG), max_age_DEG, function(genes) {
  avg_expr_age[genes, , drop = FALSE]
})

scaled_expr_DEG_list <- lapply(avg_expr_DEG_list, function(x) t(scale(t(x))))

par(mfrow = c(3, 5), mar = c(3, 3, 3, 3))

for (age in age_levels) {
  if (!is.null(scaled_expr_DEG_list[[age]])) {
    boxplot(
      scaled_expr_DEG_list[[age]],
      names = colnames(scaled_expr_DEG_list[[age]]),
      main  = paste0(age, " (", nrow(scaled_expr_DEG_list[[age]]), " genes)"),
      xlab  = "Age",
      ylab  = "Scaled expression"
    )
  }
}

cat("Peak-age expression boxplots generated.\n")
cat("=== PIPELINE COMPLETE ===\n")



library(clusterProfiler)
library(enrichplot)  # for emapplot, cnetplot
library(DOSE)        # required by enrichplot for similarity



for (cl in 1:k) {
  ego <- ego_list[[cl]]
  
  # skip clusters with no significant GO terms
  if (is.null(ego) || nrow(ego) == 0) {
    message("Cluster ", cl, ": no significant GO terms, skipping.")
    next
  }
  
  message("Cluster ", cl, ": plotting GO network...")
  
  # 1) compute pairwise semantic similarity between GO terms
  ego_sim <- pairwise_termsim(ego)
  
  # 2a) GO-term network (like the slide: nodes = GO BP terms)
  p_go <- emapplot(
    ego_sim,
    showCategory = 20,   # number of GO terms to display
    layout       = "kk"  # nice force-directed layout
  ) + ggtitle(paste("Network of GO BP – Cluster", cl))
  
  print(p_go)
  
  # 2b) Gene–GO network (nodes = genes + GO terms), colored by fold change
  #     Comment out 'foldChange' if you don't have gene_fc.
  p_gene <- cnetplot(
    ego_sim,
    showCategory = 15 # optional
  ) + ggtitle(paste("Gene–GO network – Cluster", cl))
  
  print(p_gene)
}
