# Run PCA analysis in TRACERx Renal samples

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(here)
library(factoextra)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")
OUT_DIR <- file.path(BASE, "analysis", "outputs")

# LOAD DATA ---------------------------------------------------------------

vst <- readRDS(VST_PATH)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")

# As in TRACERx Lung (Martínez-Ruiz et al, Nature 2023),
# PCA performed with VST counts centering the data but not scaling
# as VST already scales expression
tx_pca <- prcomp(t(assay(vst)), center = TRUE, scale = FALSE)
saveRDS(tx_pca, file.path(OUT_DIR, "tx_transcriptional_pca.rds"))
fviz_eig(tx_pca)
# Scree plot suggests that 5 first components explain a high proportion of transcriptional variation
# Let's use those for downstream analyses
res.pca <- summary(tx_pca)
res.pca$importance[3, 5]
# In particular, 43% of the total variance is explained by these 5 principal components
