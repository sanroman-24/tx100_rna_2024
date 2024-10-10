# Check clustering of mets by expression and how it relates to anatomical site / patient of origin
rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(here)

# PATHS
BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")

# LOAD DATA
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
vst <- readRDS(VST_PATH)

# FUNCTIONS
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "plotting_theme.R"))


# CLUSTERING PER PATIENT
met_patients <- annotation$Patient[annotation$type_collapsed != "PRIMARY"]
vst <- vst[,vst$Patient %in% met_patients]

# get I-TED 500, for example
top500_genes <- readRDS(file.path(BASE, "analysis", "outputs", "top500_variable_genes.rds"))
vst <- vst[rownames(vst) %in% top500_genes,]
vst$type_collapsed <- ifelse(
    vst$type_collapsed == "METASTASIS",
    str_remove(vst$type, "\\|.*$") %>% str_remove("_.*$"), 
    vst$type_collapsed
)
vst <- vst[,vst$type_collapsed != "#N/A"]

annotation_data <- data.frame(
    anatomical_site = vst$type_collapsed, 
    patient = vst$Patient
)

rownames(annotation_data) <- colnames(assay(vst))

palette <- brewer.pal(n = length(unique(vst$type_collapsed)), name = "Set3")
names(palette) <- unique(vst$type_collapsed)
annotation_colors <- list(anatomical_site = palette)

p <- pheatmap(
    as.matrix(t(assay(vst))), annotation_row = annotation_data, 
    annotation_colors = annotation_colors, show_rownames = FALSE, show_colnames = FALSE)

png(file.path(FIG_DIR, "heatmap_clustering_mets.png"), width = 800, height = 1000)
print(p)
dev.off()

