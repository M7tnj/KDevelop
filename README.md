## Kidney Transcriptomic Analysis — DEG Pipeline

## Overview

Reproducible R pipeline for differential gene expression analysis of **GSE8611** (Affymetrix Human Genome U133A Array), comparing **renal tubular progenitor cells** (n=3) vs **adult renal proximal tubule epithelial cells** (n=3).

This README covers the pipeline from raw data download through DEG identification and visualization.

---

## Dataset

| Field | Value |
|-------|-------|
| GEO Accession | GSE8611 |
| Platform | Affymetrix Human Genome U133A Array |
| Organism | Homo sapiens |
| Group 1 | Renal tubular progenitor cells (n=3) |
| Group 2 | Adult renal proximal tubule epithelial cells (n=3) |
| CEL files | 6 total |

---

## Pipeline Summary

```
┌─────────────────────────┐
│  1. Data Download       │  GEOquery → raw CEL files
├─────────────────────────┤
│  2. RMA Normalization   │  affy::rma() → expr_mat
├─────────────────────────┤
│  3. Quality Control     │  PCA, RLE, NUSE, MA plots
├─────────────────────────┤
│  4. Probe Annotation    │  hgu133a.db → Symbol, EntrezID
│     - Remove cross-hyb  │  (probes with "///" in Symbol)
│     - Remove no-symbol  │  (probes with NA mapping)
├─────────────────────────┤
│  5. DEG Analysis        │  limma + eBayes at probe level
│     - Threshold         │  |log2FC| > 1, raw P < 0.05
│     - FDR computed      │  BH adjustment (reported, not gated)
├─────────────────────────┤
│  6. Collapse to Genes   │  Best probe per gene (lowest P)
├─────────────────────────┤
│  7. Visualization       │  Volcano, MA, Scatter, Heatmap,
│                         │  Forest plot, Boxplots
└─────────────────────────┘
```

---

## Required R Packages

### Bioconductor
```r
BiocManager::install(c(
  "affy",          # RMA normalization
  "affyPLM",       # RLE/NUSE quality metrics
  "limma",         # Differential expression
  "hgu133a.db",    # Probe annotation (HG-U133A)
  "hgu133acdf",    # CDF environment
  "GEOquery",      # GEO data download
  "AnnotationDbi"  # Annotation infrastructure
))
```

### CRAN
```r
install.packages(c(
  "ggplot2", "pheatmap", "ggrepel", "dplyr", "tidyr",
  "ggpubr", "patchwork", "RColorBrewer", "scales"
))
```

---

## Step-by-Step Code

### 1. Load Packages

```r
library(affy)
library(affyPLM)
library(limma)
library(hgu133a.db)
library(hgu133acdf)
library(GEOquery)
library(AnnotationDbi)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

cat("hgu133a.db version:", as.character(packageVersion("hgu133a.db")), "\n")
```

---

### 2. Download Raw CEL Files

```r
gse_id  <- "GSE8611"
cel_dir <- "raw_data"
dir.create(cel_dir, showWarnings = FALSE, recursive = TRUE)

getGEOSuppFiles(GEO = gse_id, makeDirectory = TRUE, baseDir = cel_dir)

tar_file <- list.files(file.path(cel_dir, gse_id), pattern = "\\.tar$", full.names = TRUE)
untar(tar_file, exdir = file.path(cel_dir, gse_id, "CEL"))
cel_path <- file.path(cel_dir, gse_id, "CEL")

cel_files <- list.files(cel_path, pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE)
cat("Found", length(cel_files), "CEL files.\n")
```

> **Alternative:** Download manually from  
> https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8611/suppl/

---

### 3. RMA Normalization

```r
raw_affy <- ReadAffy(filenames = cel_files, cdfname = "hgu133acdf")
rma_eset <- rma(raw_affy)

expr_mat <- exprs(rma_eset)
cat("Dimensions:", dim(expr_mat), "\n")
```

RMA performs three steps:
- **Background correction** — adjusts for non-specific binding
- **Quantile normalization** — makes distributions identical across arrays
- **Median polish summarization** — combines probe pairs into one probe-set value

---

### 4. Assign Sample Groups

```r
sample_names <- sampleNames(rma_eset)

group <- factor(
  rep(c("Progenitor", "Adult_RPTE"), each = 3),
  levels = c("Progenitor", "Adult_RPTE")
)

pheno_data <- data.frame(Sample = sample_names, Group = group, stringsAsFactors = FALSE)
pData(rma_eset) <- pheno_data
```

---

### 5. Quality Control

#### 5a. RLE Plot (Relative Log Expression)

```r
plm_fit <- fitPLM(raw_affy)

png("QC_RLE_plot.png", width = 1000, height = 600, res = 150)
par(mar = c(7, 5, 3, 2))
RLE(plm_fit, main = "RLE Plot - GSE8611", las = 2, cex.axis = 0.7)
dev.off()
```

**Interpretation:** All samples should have medians near 0 and similar spread. Outliers have shifted medians or wider boxes.

#### 5b. NUSE Plot (Normalized Unscaled Standard Errors)

```r
png("QC_NUSE_plot.png", width = 1000, height = 600, res = 150)
par(mar = c(7, 5, 3, 2))
NUSE(plm_fit, main = "NUSE Plot - GSE8611", las = 2, cex.axis = 0.7)
dev.off()
```

**Interpretation:** Medians should be near 1.0. Values > 1.2 indicate problematic samples.

#### 5c. PCA Plot

```r
pca_result <- prcomp(t(expr_mat), scale. = TRUE, center = TRUE)
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = group,
  Sample = sample_names
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 5, shape = 21, stroke = 1.5) +
  geom_text_repel(aes(label = Sample), size = 3) +
  scale_color_manual(values = c(Progenitor = "#E64B35", Adult_RPTE = "#4DBBD5")) +
  scale_fill_manual(values = c(Progenitor = "#E64B3580", Adult_RPTE = "#4DBBD580")) +
  labs(title = "PCA of GSE8611 Samples",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_bw()
```

**Expected result:** The two groups should separate along PC1, confirming biological differences dominate over technical variation.

#### 5d. Sample Clustering

```r
png("QC_Sample_Clustering.png", width = 800, height = 600, res = 150)
dist_mat <- dist(t(expr_mat), method = "euclidean")
hc <- hclust(dist_mat, method = "average")
plot(hc, main = "Hierarchical Clustering (RMA)", xlab = "", sub = "", cex = 0.8)
dev.off()
```

---

### 6. Probe-to-Gene Annotation

```r
probe_ids <- rownames(expr_mat)

symbol_map  <- mapIds(hgu133a.db, keys = probe_ids, column = "SYMBOL",
                       keytype = "PROBEID", multiVals = "first")
entrez_map  <- mapIds(hgu133a.db, keys = probe_ids, column = "ENTREZID",
                       keytype = "PROBEID", multiVals = "first")
```

#### Remove cross-hybridizing probes

```r
# Probes mapping to multiple genes contain "///" in Symbol
symbol_map[grepl("///", symbol_map)] <- NA
entrez_map[grepl("///", symbol_map)] <- NA

# Remove probes with no gene symbol
valid <- !is.na(symbol_map)
cat("Total probes:", length(probe_ids), "\n")
cat("After removing unmapped & cross-hybridizing:", sum(valid), "\n")
```

---

### 7. Differential Expression Analysis

#### 7a. Design Matrix & Linear Model

```r
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_mat, design)

contrast_mat <- makeContrasts(
  Adult_vs_Prog = Adult_RPTE - Progenitor,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)
```

**Why eBayes?** With only n=3 per group, individual gene variance estimates are unstable. eBayes borrows information across all genes to shrink variances toward a common value, producing more stable and reliable test statistics — this is standard best practice for small-sample microarray analysis.

#### 7b. Extract Probe-Level Results

```r
results_probes <- data.frame(
  ProbeID   = rownames(fit2),
  logFC     = as.numeric(fit2$coefficients[, "Adult_vs_Prog"]),
  AveExpr   = as.numeric(fit2$Amean),
  t         = as.numeric(fit2$t[, "Adult_vs_Prog"]),
  P.Value   = as.numeric(fit2$p.value[, "Adult_vs_Prog"]),
  SE        = as.numeric(fit2$stdev.unscaled[, "Adult_vs_Prog"] * sqrt(fit2$s2.post)),
  df        = as.numeric(fit2$df.total),
  adj.P.Val = p.adjust(as.numeric(fit2$p.value[, "Adult_vs_Prog"]), method = "BH"),
  Symbol    = symbol_map[rownames(fit2)],
  EntrezID  = entrez_map[rownames(fit2)],
  stringsAsFactors = FALSE
)

# Remove unmapped and cross-hybridizing probes
results_probes <- results_probes[!is.na(results_probes$Symbol), ]
results_probes <- results_probes[!grepl("///", results_probes$Symbol), ]
```

#### 7c. Apply DEG Thresholds

```r
deg_probes <- results_probes %>%
  filter(abs(logFC) > 1 & P.Value < 0.05) %>%
  mutate(Direction = ifelse(logFC > 0, "Upregulated", "Downregulated")) %>%
  arrange(P.Value)

cat("Total DEG probes:   ", nrow(deg_probes), "\n")
cat("  Upregulated:      ", sum(deg_probes$Direction == "Upregulated"), "\n")
cat("  Downregulated:    ", sum(deg_probes$Direction == "Downregulated"), "\n")
```

**Threshold justification:**
- `|log2FC| > 1` corresponds to linear fold change > 2 (matching the article's stated |FC| > 2)
- `P.Value < 0.05` (raw p-value, matching the original TAC analysis)
- FDR (adj.P.Val) is computed and reported for reviewer requirements but not used as the primary gate

#### 7d. Collapse to Gene-Level (Best Probe per Gene)

```r
deg_genes <- deg_probes %>%
  group_by(Symbol) %>%
  arrange(P.Value) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    Gene = Symbol,
    Direction = ifelse(logFC > 0, "Upregulated", "Downregulated")
  )
```

#### 7e. All-Gene Reference Table (for plots & downstream)

```r
results_all_gene <- results_probes %>%
  group_by(Symbol) %>%
  arrange(P.Value) %>%
  slice(1) %>%
  ungroup()

# Build expression matrix with gene symbols as rownames
expr_collapsed <- as.matrix(expr_mat[results_all_gene$ProbeID, , drop = FALSE])
rownames(expr_collapsed) <- results_all_gene$Symbol
```

#### 7f. Effect Sizes with 95% Confidence Intervals

```r
df_common <- fit2$df.total[1]
t_crit <- qt(0.975, df = df_common)

results_all_df <- deg_genes %>%
  as.data.frame() %>%
  mutate(
    EntrezID = results_all_gene$EntrezID[match(Gene, results_all_gene$Symbol)],
    SE       = results_all_gene$SE[match(Gene, results_all_gene$Symbol)],
    CI_lower = logFC - t_crit * SE,
    CI_upper = logFC + t_crit * SE,
    Direction = case_when(
      abs(logFC) > 1 & P.Value < 0.05 & logFC > 0 ~ "Upregulated",
      abs(logFC) > 1 & P.Value < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

write.csv(results_all_df, "DEGs.csv", row.names = FALSE)
```

---

### 8. Visualization

#### 8a. Volcano Plot

```r
results_all_plot <- results_all_gene %>%
  as.data.frame() %>%
  mutate(
    Gene = Symbol,
    Direction = case_when(
      abs(logFC) > 1 & P.Value < 0.05 & logFC > 0 ~ "Upregulated",
      abs(logFC) > 1 & P.Value < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

p_volcano <- ggplot(results_all_plot, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = filter(results_all_plot, Direction == "Not Significant"),
             color = "grey70", size = 0.8, alpha = 0.4) +
  geom_point(data = filter(results_all_plot, Direction == "Upregulated"),
             color = "#E64B35", size = 1.2, alpha = 0.7) +
  geom_point(data = filter(results_all_plot, Direction == "Downregulated"),
             color = "#4DBBD5", size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = filter(results_all_plot,
                  (abs(logFC) > 3 & P.Value < 1e-4) |
                    Gene %in% c("CDC42","CYCS","CAT","PIK3R1","FOXO1","NRAS","PPARGC1A","APOE")),
    aes(label = Gene), size = 2.8, max.overlaps = 25
  ) +
  labs(title = "Volcano Plot: Adult RPTE vs Tubular Progenitor Cells",
       subtitle = "Threshold: |log2FC| > 1, P < 0.05",
       x = "log2 Fold Change", y = "-log10(P-value)") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Volcano_Plot.png", p_volcano, width = 10, height = 8, dpi = 300)
```

**Reading the volcano plot:**
- X-axis: magnitude of expression change (log2 fold change)
- Y-axis: statistical significance (-log10 p-value)
- Red dots: upregulated in adult RPTE (logFC > 1)
- Blue dots: downregulated in adult RPTE (logFC < -1)
- Grey dots: not significant
- Dashed lines mark the thresholds (|log2FC| = 1, P = 0.05)
- Labeled points: hubs-of-hubs genes + extreme outliers

#### 8b. MA Plot

```r
p_ma <- ggplot(results_all_plot, aes(x = AveExpr, y = logFC)) +
  geom_point(data = filter(results_all_plot, Direction == "Not Significant"),
             color = "grey70", size = 0.6, alpha = 0.4) +
  geom_point(data = filter(results_all_plot, Direction != "Not Significant"),
             aes(color = Direction), size = 1, alpha = 0.7) +
  scale_color_manual(values = c(Upregulated = "#E64B35", Downregulated = "#4DBBD5")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(title = "MA Plot: Adult RPTE vs Tubular Progenitor Cells",
       x = "Average log2 Expression", y = "log2 Fold Change") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")

ggsave("MA_Plot.png", p_ma, width = 10, height = 8, dpi = 300)
```

**Reading the MA plot:**
- X-axis: average expression level across all samples
- Y-axis: log fold change between groups
- Reveals intensity-dependent bias (should be symmetric around 0)
- Most DEGs are upregulated (red) at higher expression levels

#### 8c. Scatter Plot

```r
avg_prog  <- rowMeans(expr_collapsed[, group == "Progenitor"])
avg_adult <- rowMeans(expr_collapsed[, group == "Adult_RPTE"])

scatter_df <- data.frame(
  Gene = names(avg_prog),
  Progenitor = avg_prog,
  Adult_RPTE = avg_adult,
  Direction = results_all_plot$Direction[match(names(avg_prog), results_all_plot$Gene)]
)

p_scatter <- ggplot(scatter_df, aes(x = Progenitor, y = Adult_RPTE)) +
  geom_point(data = filter(scatter_df, Direction == "Not Significant" | is.na(Direction)),
             color = "grey70", size = 0.6, alpha = 0.4) +
  geom_point(data = filter(scatter_df, !is.na(Direction) & Direction != "Not Significant"),
             aes(color = Direction), size = 1, alpha = 0.7) +
  scale_color_manual(values = c(Upregulated = "#E64B35", Downregulated = "#4DBBD5")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "Scatter Plot: Mean Expression by Group",
       x = "Mean Expression (Tubular Progenitor)",
       y = "Mean Expression (Adult RPTE)") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")

ggsave("Scatter_Plot.png", p_scatter, width = 8, height = 8, dpi = 300)
```

**Reading the scatter plot:**
- Each point is a gene
- Diagonal line = no change between groups
- Points above diagonal: upregulated in adult RPTE
- Points below diagonal: downregulated in adult RPTE

#### 8d. Hierarchical Clustering Heatmap

```r
deg_gene_symbols <- results_all_df %>%
  filter(Direction != "Not Significant") %>%
  pull(Gene)
deg_gene_symbols <- intersect(deg_gene_symbols, rownames(expr_collapsed))

deg_expr <- expr_collapsed[deg_gene_symbols, , drop = FALSE]
deg_z <- t(scale(t(deg_expr)))  # Z-score per gene

annotation_col <- data.frame(Group = group, row.names = colnames(deg_z))
ann_colors <- list(Group = c(Progenitor = "#E64B35", Adult_RPTE = "#4DBBD5"))

pheatmap(
  deg_z,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "none",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "DEG Heatmap: Adult RPTE vs Progenitor\n(|log2FC| > 1, P < 0.05)",
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  fontsize = 10,
  filename = "Heatmap_DEG.png",
  width = 10, height = 12, dpi = 300
)
```

**Reading the heatmap:**
- Rows = DEGs (Z-score normalized), clustered by similarity
- Columns = samples, color-coded by group
- Red = high expression relative to mean, Blue = low expression
- Clear separation between Progenitor and Adult_RPTE samples confirms DEGs discriminate the groups

#### 8e. Forest Plot (Hubs of Hubs with 95% CIs)

```r
hub_of_hubs_genes <- c("CDC42", "CYCS", "CAT", "PIK3R1", "FOXO1", "NRAS", "PPARGC1A", "APOE")

hoh_results <- results_all_df %>%
  filter(Gene %in% hub_of_hubs_genes) %>%
  arrange(desc(logFC)) %>%
  mutate(Gene = factor(Gene, levels = Gene))

p_forest <- ggplot(hoh_results, aes(x = logFC, y = Gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, linewidth = 0.8) +
  geom_point(size = 4, color = ifelse(hoh_results$logFC > 0, "#E64B35", "#4DBBD5")) +
  labs(title = "Hubs of Hubs: Effect Sizes with 95% CIs",
       subtitle = "Adult RPTE vs Tubular Progenitor Cells",
       x = "log2 Fold Change (95% CI)", y = "") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Forest_Plot_HubsOfHubs.png", p_forest, width = 10, height = 6, dpi = 300)
```

**Reading the forest plot:**
- Points = estimated log2 fold change for each hub gene
- Horizontal bars = 95% confidence interval
- Bars that do not cross 0 indicate statistically significant effects
- Red = upregulated, Blue = downregulated

#### 8f. Expression Boxplots (Hubs of Hubs)

```r
hoh_genes_found <- intersect(hub_of_hubs_genes, rownames(expr_collapsed))
hoh_plot_list <- list()

for (gene in hoh_genes_found) {
  plot_df <- data.frame(
    Expression = expr_collapsed[gene, ],
    Group = group,
    Sample = sample_names
  )
  p <- ggboxplot(plot_df, x = "Group", y = "Expression",
                 color = "Group", palette = c("#E64B35", "#4DBBD5"),
                 add = "jitter") +
    stat_compare_means(method = "t.test", label = "p.signif") +
    labs(title = gene, y = "log2 Expression") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
  hoh_plot_list[[gene]] <- p
}

p_combined <- wrap_plots(hoh_plot_list, ncol = 4)
ggsave("Boxplots_HubsOfHubs.png", p_combined, width = 16, height = 10, dpi = 300)
```

---

## Output Files

| File | Description |
|------|-------------|
| `DEGs.csv` | Gene-level DEGs with logFC, P, FDR, SE, 95% CI, Direction |
| `DEG_Gene_Level2.1.csv` | Probe-level DEGs (multiple probes per gene) |
| `Volcano_Plot.png` | Volcano plot of all genes |
| `MA_Plot.png` | MA plot (mean expression vs fold change) |
| `Scatter_Plot.png` | Group mean scatter plot |
| `Heatmap_DEG.png` | Hierarchical clustering heatmap of DEGs |
| `Forest_Plot_HubsOfHubs.png` | Effect sizes with 95% CIs for 8 hub genes |
| `Boxplots_HubsOfHubs.png` | Per-gene expression boxplots |
| `QC_RLE_plot.png` | RLE quality control plot |
| `QC_NUSE_plot.png` | NUSE quality control plot |
| `QC_Sample_Clustering.png` | Sample dendrogram |

---

## Key Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Normalization | RMA | Standard for Affymetrix; background correct + quantile normalize + summarize |
| DE method | limma + eBayes | Best practice for small-n microarray; moderates variances |
| Analysis level | Probe-level (then collapsed to gene) | Matches TAC approach; each probe tested independently |
| logFC threshold | \|log2FC\| > 1 | Equivalent to linear FC > 2 (article's stated threshold) |
| P-value threshold | Raw P < 0.05 | Matches original TAC analysis |
| FDR method | Benjamini-Hochberg | Computed and reported for reviewer; not used as primary gate |
| Probe collapse | Lowest P per gene | Most significant probe retained for each gene symbol |
| Cross-hybridizing probes | Removed (Symbol contains "///") | Prevents ambiguous gene assignments |
| Annotation | hgu133a.db | Standard Bioconductor annotation for HG-U133A |

---

## Notes on TAC vs R Comparison

The original article used **Transcriptome Analysis Console (TAC) v4.0** with the same thresholds and reported **1439 DEGs** (1060 up, 379 down). The R pipeline may produce a different count due to:

| Factor | TAC | R (this pipeline) |
|--------|-----|-------------------|
| RMA implementation | Proprietary C/C++ | Bioconductor `affy` package |
| Fold-change estimator | Tukey biweight on probe-pair differences | Least-squares coefficient (group mean diff) |
| P-value computation | ANOVA F-test with internal variance stabilization | limma empirical Bayes moderated t-test |
| Probe filtering | May include all probe sets | Removes unmapped & cross-hybridizing |

Small numerical differences in expression values cascade into threshold-dependent DEG counts. The **key biological findings** (hub genes, top pathways) should be consistent across both approaches.

---

## Reproducibility

```r
sessionInfo()
# Save at the end of your analysis to record all package versions
sink("Session_Info.txt")
sessionInfo()
sink()
```
