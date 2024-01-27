# Analysis of co-occurrence of 9p and 14q loss at the single-cell level

rm(list = ls(all = T))

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


# LOAD DATA ---------------------------------------------------------------
srat = readRDS(SEURAT_PATH)
meta = srat@meta.data
# Filter to only sporadic ccRCC
meta = meta[!is.na(meta$Disease_type) & meta$Disease_type == "ccRCC_sporadic", ]

# 9p loss is rather infrequent
table(meta$Chr14q23_31_Status)
table(meta$Chr9p21.3_Status)

# 
meta %>% filter(Chr14q23_31_Status != "Gain", Chr9p21.3_Status != "Gain") %>% 
  {table(.$Chr14q23_31_Status, .$Chr9p21.3_Status)}
