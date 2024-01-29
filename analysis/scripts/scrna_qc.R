# QC of the single-cell RNA-Seq metadata


# PACKAGES ----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(here)
library(nlme)
library(ggpubr)
library(ggbeeswarm)

# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis")
PLOT_DIR = file.path(OUT_DIR, "figures", "scrna") # TODO tmp
SEURAT_PATH = file.path(BASE, "data", "processed", "Chr3p_Subset_Merged_Clustered.rds")


# FUNCTIONS ---------------------------------------------------------------

filter_adaptive_thresh = function(srat, v, up = T){
  if (up){s = 1}else{s = -1}
  v = as_vector(srat[[v]])
  # Adaptive threshold: 3 MAD away from median
  thrsh = median(v) +  s*3*mad(v)
  if (up){idx = which(v <= thrsh)} else {idx = which(v > thrsh)}
  return (srat[,idx])
}

# LOAD DATA ---------------------------------------------------------------
srat = readRDS(SEURAT_PATH) # 50711 cells, 13830 genes
srat = filter_adaptive_thresh(srat, "percent.ribo") # 50344 cells
srat = filter_adaptive_thresh(srat, "percent.mito") # 50344 cells
# For the number of RNA features we know from Geoffrey that at least 1000 is OK
srat = srat[,srat$nFeature_RNA > 1000] # 38411 cells

# remove poorly expressed genes
# keep only genes that are detected in at least 0.5% of the cells
keep_genes = rowSums(srat@assays$RNA$counts > 1) >= 0.005*ncol(srat)
srat = srat[which(keep_genes), ]

saveRDS(srat, file.path(OUT_DIR, "outputs", "scrna",  "chr3p_subset_qc_exp_filt_merged_clustered.rds"))


