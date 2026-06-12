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

