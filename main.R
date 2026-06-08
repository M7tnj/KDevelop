###########################
#### Packages Required ###
#########################

library(BiocManager)
library(affy)
library(affyPLM)
library(limma)
library(hgu133a.db)
library(hgu133acdf)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ReactomePA)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(data.table)
library(igraph)
library(ggraph)
library(ggsci)
library(corrplot)
library(patchwork)
library(ggpubr)
library(cluster)
library(dendextend)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(igraph)
library(ggraph)
library(GEOquery)
library(ggplot2)
library(shadowtext)
library(tibble)
library(tidyverse)
library(ggridges)
library(ReactomePA)
library(stringr)
library(clusterProfiler)
library(tabbycat)
library(igraph)
library(STRINGdb)
library(patchwork)


######################################
##### Download & Import Samples #####
####################################

dirs <- c("results", "results/QC", "results/DEG", "results/PPI",
          "results/Enrichment", "results/Modules", "results/GSEA",
          "results/HubsOfHubs")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

gse_id  <- "GSE8611"


getGEOSuppFiles(GEO = gse_id, makeDirectory = TRUE)

tar_file <- list.files(file.path( gse_id), pattern = "\\.tar$", full.names = TRUE)

if (length(tar_file) > 0) {
  untar(tar_file, exdir = file.path(gse_id, "CEL"))
  cel_path <- file.path( gse_id, "CEL")
} else {
  cel_path <- file.path(gse_id)
}

all <- list.files("GSE8611/CEL", pattern="\\.CEL\\.gz$", full.names=TRUE, ignore.case=TRUE)
cat("Total files found:", length(all), "\n")

# Files of interest
keep <- c("69","59","63","19","23","25")
pattern <- paste0("GSM\\d*(",paste(keep,collapse="|"),")\\.CEL\\.gz$")
del <- all[!grepl(pattern, basename(all), ignore.case=TRUE)]
file.remove(del)

########################
### Data Processing ###
######################

raw_affy <- ReadAffy(celfile.path = cel_path, cdfname = "hgu133acdf")
rma_eset <- rma(raw_affy) 

gsm_info <- getGEO(GEO = gse_id, GSEMatrix = TRUE)
sample_names <- sampleNames(rma_eset)
group <- factor(
  rep(c("Progenitor", "Adult_RPTE"), each = 3),
  levels = c("Progenitor", "Adult_RPTE")
)


pheno_data <- data.frame(
  Sample = sample_names,
  Group = group,
  stringsAsFactors = FALSE
)

pData(rma_eset) <- pheno_data

print(table(group))


#########################
#### Quality Control ###
#######################

plm_fit <- fitPLM(raw_affy)
png("results/QC/01_RLE_plot.png",
    width = 1080, 
    height = 720, 
    units = "px", 
    res = 300)
par(mar = c(7, 5, 3, 2))
RLE(plm_fit, main = "RLE Plot - GSE8611", las = 2, cex.axis = 0.7)
dev.off()

png("results/QC/02_NUSE_plot.png",
    width = 1080, 
    height = 720, 
    res = 300)
par(mar = c(7, 5, 3, 2))
NUSE(plm_fit, main = "NUSE Plot - GSE8611", las = 2, cex.axis = 0.7)
dev.off()

nuse_medians <- NUSE(plm_fit, type = "stats")["median", ]
rle_medians  <- RLE(plm_fit, type = "stats")["median", ]
cat("\nNUSE medians: should be close to 1.0 >>>\n")
print(round(nuse_medians, 4))
cat("\nRLE medians: should be close to 0 >>>:\n")
print(round(rle_medians, 4))

# PCA
expr_mat <- exprs(rma_eset)
pca_result <- prcomp(t(expr_mat), scale. = TRUE, center = TRUE)
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = group,
  Sample = sample_names
)
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 5, shape = 21, stroke = 1.5) +
  
  geom_shadowtext(aes(label = Sample), 
                  size = 3,
                  bg.colour = "white",    
                  bg.r = 0.2,              
                  check_overlap = TRUE,     
                  fontface = "plain") +     
  
  scale_color_manual(values = c("Progenitor" = "#E64B35", 
                                "Adult_RPTE" = "#4DBBD5")) +
  scale_fill_manual(values = c("Progenitor" = "#E64B3580", 
                               "Adult_RPTE" = "#4DBBD580")) +
  labs(
    title = "PCA of GSE8611 Samples (RMA-normalized)",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )
ggsave("results/QC/03_PCA_plot.png", p_pca, width = 8, height = 6, dpi = 300)


## Hierarchical Clustering Heat map
png("results/QC/04_Sample_Clustering.png", 
    width = 800, 
    height = 600, 
    res = 150)
dist_mat <- dist(t(expr_mat), method = "euclidean")
hc <- hclust(dist_mat, method = "average")
plot(hc, main = "Hierarchical Clustering of Samples (RMA)", xlab = "", sub = "",
     cex = 0.8, labels = sample_names)
dev.off()

## MA plot of all samples & refrence
raw_expr <- log2(exprs(raw_affy))

A <- rowMeans(raw_expr, na.rm = TRUE)
array_idx <- 1              
M <- raw_expr[, array_idx] - A

valid_idx <- is.finite(A) & is.finite(M)
A <- A[valid_idx]
M <- M[valid_idx]

png("results/QC/05_MA_plot_raw.png",
    width  = 1000,
    height = 600,
    res    = 150,
    units  = "px")
par(mfrow = c(1, 2))
plot(A, M, cex = 0.3, pch = ".",
     xlab = "A = mean expression (log2)",
     ylab = "M = log2 ratio",
     main = "MA Plot (Raw)")
abline(h = 0, col = "red", lty = 2)
norm_expr <- log2(exprs(rma_eset))
A_norm <- rowMeans(norm_expr, na.rm = TRUE)
M_norm <- norm_expr[, array_idx] - A_norm

valid_idx_norm <- is.finite(A_norm) & is.finite(M_norm)
A_norm <- A_norm[valid_idx_norm]
M_norm <- M_norm[valid_idx_norm]

plot(A_norm, M_norm, cex = 0.3, pch = ".",
     xlab = "A = mean expression (log2)",
     ylab = "M = log2 ratio",
     main = "MA Plot (RMA Normalised)")
abline(h = 0, col = "red", lty = 2)
dev.off()

###############################
#### Annotation & mapping ####
#############################

probe_ids <- rownames(expr_mat)

symbol_map <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first" 
)

entrez_map <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  column = "ENTREZID",
  keytype = "PROBEID",
  multiVals = "first"
)

genename_map <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  column = "GENENAME",
  keytype = "PROBEID",
  multiVals = "first"
)

multi_symbol <- mapIds(
  hgu133a.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "list"
)

cross_hybrid <- sapply(multi_symbol, function(x) {
  if (is.na(x[1])) return(FALSE)
  length(x) > 1
})

multi_map <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "PROBEID"
)

annot_df <- data.frame(
  ProbeID = probe_ids,
  Symbol  = symbol_map,
  EntrezID = entrez_map,
  GeneName = genename_map,
  stringsAsFactors = FALSE
)
annot_df <- annot_df[!is.na(annot_df$Symbol), ]
annot_df <- annot_df[!grepl("///", annot_df$Symbol), ]

expr_annot <- expr_mat[annot_df$ProbeID, , drop = FALSE]

probe_iqr <- apply(expr_annot, 1, IQR)
annot_df$IQR <- probe_iqr[match(annot_df$ProbeID, names(probe_iqr))]

annot_df_collapsed <- annot_df %>%
  group_by(Symbol) %>%
  dplyr::slice_max(order_by = IQR, n = 1, with_ties = FALSE) %>%
  ungroup()
cat("Unique genes after collapsing:", nrow(annot_df_collapsed), "\n")

expr_collapsed <- expr_mat[annot_df_collapsed$ProbeID, ]
rownames(expr_collapsed) <- annot_df_collapsed$Symbol
write.csv(annot_df_collapsed, "results/QC/Probe_to_Gene_Mapping.csv", row.names = FALSE)

#######################
#### DEG Analysis ####
#####################

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
print(design)

fit <- lmFit(expr_mat, design)
contrast_mat <- makeContrasts(
  Adult_vs_Prog = Adult_RPTE - Progenitor,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

se_probe    <- fit2$stdev.unscaled[, "Adult_vs_Prog"] * fit2$sigma
df_residual <- fit2$df.residual   
tstat_probe <- fit2$coefficients[, "Adult_vs_Prog"] / se_probe
pval_std    <- 2 * pt(-abs(tstat_probe), df = df_residual)

results_probes <- data.frame(
  ProbeID  = rownames(fit2),
  logFC    = as.numeric(fit2$coefficients[, "Adult_vs_Prog"]),
  AveExpr  = as.numeric(fit2$Amean),
  t        = as.numeric(fit2$t[, "Adult_vs_Prog"]),
  P.Value  = as.numeric(fit2$p.value[, "Adult_vs_Prog"]),
  SE       = as.numeric(fit2$stdev.unscaled[, "Adult_vs_Prog"] * sqrt(fit2$s2.post)),
  df       = as.numeric(fit2$df.total),
  adj.P.Val = as.numeric(fit2$p.value[, "Adult_vs_Prog"]),
  stringsAsFactors = FALSE
)
results_probes$adj.P.Val <- p.adjust(results_probes$P.Value, method = "BH")

results_probes$Symbol   <- symbol_map[results_probes$ProbeID]
results_probes$EntrezID <- entrez_map[results_probes$ProbeID]

## Remove probes without gene symbols and cross-hybridizing probes
results_probes <- results_probes[!is.na(results_probes$Symbol), ]
results_probes <- results_probes[!grepl("///", results_probes$Symbol), ]

deg_probes <- results_probes %>%
  filter(abs(logFC) > 1 & P.Value < 0.05) %>%
  mutate(Direction = ifelse(logFC > 0, "Upregulated", "Downregulated"))

cat("  Total probes:      ", nrow(deg_probes), "\n")
cat("  Upregulated:       ", sum(deg_probes$Direction == "Upregulated"), "\n")
cat("  Downregulated:     ", sum(deg_probes$Direction == "Downregulated"), "\n")

deg_genes <- deg_probes %>%
  dplyr::select(ProbeID, Symbol, EntrezID, logFC, AveExpr, t, P.Value, adj.P.Val, SE, Direction)
write.csv(deg_genes, "results/DEG/DEG Gene Level.csv", row.names = FALSE)

deg_genes <- deg_probes %>%
  group_by(Symbol) %>%
  arrange(P.Value) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(
    Gene = Symbol,
    Direction = ifelse(logFC > 0, "Upregulated", "Downregulated")
  )
df_common <- df_residual[1]
t_crit <- qt(0.975, df = df_common)


results_all_gene <- results_probes %>%
  group_by(Symbol) %>%
  arrange(P.Value) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(ProbeID, Symbol, everything())

results_all_gene <- results_all_gene %>%
  dplyr::select(Symbol, EntrezID, SE)

# 1. Convert to tibble and add row names as ProbeID
results_all_df <- deg_genes %>%
  as_tibble() %>%
  rownames_to_column("RowName")

results_all_df <- results_all_df %>%
  left_join(
    results_all_gene %>% dplyr::select(Symbol, EntrezID, SE),
    by = c("RowName" = "Symbol")
  )

results_all_df <- results_all_df %>%
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
write.csv(results_all_df, "results/DEG/DEGs.csv", row.names = FALSE)


############
############
up_genes <- results_all_df %>%
  filter(Direction == "Upregulated") %>%
  pull(Gene)

down_genes <- results_all_df %>%
  filter(Direction == "Downregulated") %>%
  pull(Gene)

############
############
# plotting #
results_all_gene <- results_probes %>%
  group_by(Symbol) %>%
  arrange(P.Value) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::select(Symbol, EntrezID, logFC, P.Value, AveExpr, t, adj.P.Val, SE)

results_all_plot <- results_all_gene %>%
  as_tibble() %>%
  mutate(
    Gene = Symbol,
    CI_lower = logFC - t_crit * SE,
    CI_upper = logFC + t_crit * SE,
    Direction = case_when(
      abs(logFC) > 1 & P.Value < 0.05 & logFC > 0 ~ "Upregulated",
      abs(logFC) > 1 & P.Value < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )
results_all_df <- results_all_plot
p_volcano <- ggplot(results_all_plot, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(
    data = filter(results_all_plot, Direction == "Not Significant"),
    color = "grey70", size = 0.8, alpha = 0.4
  ) +
  geom_point(
    data = filter(results_all_plot, Direction == "Upregulated"),
    color = "#E64B35", size = 1.2, alpha = 0.7
  ) +
  geom_point(
    data = filter(results_all_plot, Direction == "Downregulated"),
    color = "#4DBBD5", size = 1.2, alpha = 0.7
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = filter(results_all_plot,
                  (abs(logFC) > 3 & P.Value < 1e-4) |
                    Gene %in% c("CDC42","CYCS","CAT","PIK3R1","FOXO1","NRAS","PPARGC1A","APOE")),
    aes(label = Gene), size = 2.8, max.overlaps = 25
  ) +
  labs(
    title = "Volcano Plot: Adult RPTE vs Tubular Progenitor Cells",
    subtitle = "Threshold: |log2FC| > 1, P < 0.05",
    x = "log2 Fold Change",
    y = "-log10(P-value)"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11))
ggsave("results/DEG/06_Volcano_Plot.png", p_volcano, width = 10, height = 8, dpi = 300)

p_ma <- ggplot(results_all_plot, aes(x = AveExpr, y = logFC)) +
  geom_point(
    data = filter(results_all_plot, Direction == "Not Significant"),
    color = "grey70", size = 0.6, alpha = 0.4
  ) +
  geom_point(
    data = filter(results_all_plot, Direction != "Not Significant"),
    aes(color = Direction), size = 1, alpha = 0.7
  ) +
  scale_color_manual(values = c("Upregulated" = "#E64B35", "Downregulated" = "#4DBBD5")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(
    title = "MA Plot: Adult RPTE vs Tubular Progenitor Cells",
    subtitle = "Threshold: |log2FC| > 1, P < 0.05",
    x = "Average log2 Expression",
    y = "log2 Fold Change"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")
ggsave("results/DEG/07_MA_Plot.png", p_ma, width = 10, height = 8, dpi = 300)


avg_prog  <- rowMeans(expr_collapsed[, group == "Progenitor"])
avg_adult <- rowMeans(expr_collapsed[, group == "Adult_RPTE"])
scatter_df <- data.frame(
  Gene = names(avg_prog),
  Progenitor = avg_prog,
  Adult_RPTE = avg_adult,
  Direction = results_all_df$Direction[match(names(avg_prog), results_all_plot$Gene)]
)
p_scatter <- ggplot(scatter_df, aes(x = Progenitor, y = Adult_RPTE)) +
  geom_point(
    data = filter(scatter_df, Direction == "Not Significant" | is.na(Direction)),
    color = "grey70", size = 0.6, alpha = 0.4
  ) +
  geom_point(
    data = filter(scatter_df, !is.na(Direction) & Direction != "Not Significant"),
    aes(color = Direction), size = 1, alpha = 0.7
  ) +
  scale_color_manual(values = c("Upregulated" = "#E64B35", "Downregulated" = "#4DBBD5")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title = "Scatter Plot: Mean Expression by Group",
    subtitle = "Threshold: |log2FC| > 1, P < 0.05",
    x = "Mean Expression (Tubular Progenitor)",
    y = "Mean Expression (Adult RPTE)"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

ggsave("results/DEG/08_Scatter_Plot.png", p_scatter, width = 8, height = 8, dpi = 300)


expr_collapsed <- expr_mat[results_all_gene$ProbeID, , drop = FALSE]
rownames(expr_collapsed) <- results_all_gene$Symbol

deg_gene_symbols <- results_all_df %>%
  filter(Direction != "Not Significant") %>%
  pull(Gene)

deg_gene_symbols <- intersect(deg_gene_symbols, rownames(expr_collapsed))

cat("DEG genes found in expression matrix:", length(deg_gene_symbols), "\n")

if (length(deg_gene_symbols) >= 2) {
  deg_expr <- expr_collapsed[deg_gene_symbols, , drop = FALSE]
  deg_z <- t(scale(t(deg_expr)))
  annotation_col <- data.frame(
    Group = group,
    row.names = colnames(deg_z)
  )
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
    filename = "results/DEG/Heatmap_DEG.png",
    width = 10, height = 12, dpi = 300
  )
  cat("Heatmap saved: Heatmap_DEG.png\n")
} else {
  cat("First few Gene values:", head(results_all_df$Gene), "\n")
}
hub_of_hubs_genes <- c("CDC42", "CYCS", "CAT", "PIK3R1", "FOXO1", "NRAS", "PPARGC1A", "APOE")

hoh_results <- results_all_df %>%
  filter(Gene %in% hub_of_hubs_genes) %>%
  arrange(desc(logFC)) %>%
  mutate(Gene = factor(Gene, levels = Gene))

if (nrow(hoh_results) > 0) {
  p_forest <- ggplot(hoh_results, aes(x = logFC, y = Gene)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, linewidth = 0.8) +
    geom_point(size = 4, color = ifelse(hoh_results$logFC > 0, "#E64B35", "#4DBBD5")) +
    labs(
      title = "Hubs of Hubs: Effect Sizes with 95% CIs",
      subtitle = "Adult RPTE vs Tubular Progenitor Cells",
      x = "log2 Fold Change (95% CI)",
      y = ""
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("results/DEG/Forest_Plot_HubsOfHubs.png", p_forest, width = 10, height = 6, dpi = 300)
}
