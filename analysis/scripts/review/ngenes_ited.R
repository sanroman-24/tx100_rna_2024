# Calculate I-TED scores obtained using top500 most variable genes or more
# Run on HPC -> 132GB RAM, 6 cores, +24h walltime

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(DESeq2)
library(energy)
library(pheatmap)
library(ggthemes)
library(doMC)

# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "outputs", "review")
PLOT_DIR = file.path(BASE, "analysis", "figures", "review")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")
ITED500_PATH <- file.path(BASE, "analysis", "outputs", "ITED_matrix_primary.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))


# LOAD DATA ---------------------------------------------------------------
annotation = read_delim(ANNOTATION_PATH, delim = "\t")
vst = read_rds(VST_PATH)
annotation$ITH <- as.numeric(annotation$ITH)
annotation$wgii <- as.numeric(annotation$wgii)
annotation$Regions <- as.numeric(annotation$Regions)

annotation <- annotation[annotation$sample != "K207-R4", ]
vst <- vst[,colnames(vst) != "K207-R4"]

# check that it is good to go
all(annotation$sample == colnames(vst))

# CALCULATE TRANSCRIPTIONAL DISTANCE --------------------------------------
d_mat = get_dist_matrix_pl(assay(vst), rownames(vst), f = "dcor", cores = 6)

write_rds(d_mat, file.path(OUT_DIR, "ITED_matrix_allgenes.rds"))

d_mat = read_rds(file.path(OUT_DIR, "ITED_matrix_allgenes.rds"))

annotation <- annotation %>% 
  filter(type_collapsed == "PRIMARY")

# get samples with more than one tumour primary regions
keep_patients <- annotation %>% 
  dplyr::group_by(Patient) %>% dplyr::summarise(n = n()) %>% # count number of different samples of each patient
  filter(n > 1) %>% dplyr::select(Patient) %>% as_vector() # patients with +1 primary region

keep_samples <- annotation[annotation$Patient %in% keep_patients, "sample"] %>%
  as_vector()

annotation <- annotation[annotation$sample %in% keep_samples, ]
vst <- vst[ ,colnames(vst) %in% keep_samples]

d_mat = get_dist_matrix_pl(assay(vst), rownames(vst), f = "dcor", cores = 6)
write_rds(d_mat, file.path(OUT_DIR, "ITED_matrix_primary.rds"))
