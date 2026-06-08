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




