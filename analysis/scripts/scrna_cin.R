# Association between CIN and other signatures in scRNA-seq

# Inspiration from Nature 2018 Sam Bakhoum Fig 5
# Obtain a Z-score for each cell in the study for Bakhoum sigs
# Perform correlation between CIN signature and the rest sigs
# Be aware of ICB or treatment-naive setting

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(Seurat)
library(here)


# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "outputs", "scrna")
FIG_DIR = file.path(BASE, "analysis", "figures", "scrna")
# SEURAT_PATH = file.path(BASE, "data", "processed", "Chr3p_Subset_Merged_Clustered.rds")
SEURAT_PATH =  file.path(OUT_DIR,  "chr3p_subset_qc_exp_filt_merged_clustered.rds")
META_DIR = file.path(BASE, "data", "meta")
HALLMARK_GROUPS_PATH <- file.path(META_DIR, "martinez_ruiz_2023_hallmark_gs_groups.txt")
BAKHOUM_GENES_DIR <- file.path(META_DIR, "bakhoum_genes")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

get_sign_score = function(srat, path_gs, path_name){
  idx = which(rownames(srat) %in% path_gs)
  if (length(idx) < 2){
    return (srat)
  }
  subsrat = srat[idx]
  sign_score = colMeans(subsrat@assays$RNA$scale.data)
  srat[[path_name]] = sign_score
  return(srat)
}

# LOAD DATA ---------------------------------------------------------------

srat = readRDS(SEURAT_PATH)
srat = srat[,!is.na(srat$Disease_type) & srat$Disease_type == "ccRCC_sporadic" | 
              srat$Publication == "Obradovic-etal", ]

STING_genes <- c("TMEM173")
NFKB1_genes <- c("NFKB1", "IKBKG", "TRAF6", "RELA")
NFKB2_genes <- c("NFKB2", "RELB", "MAP3K14")
STING_act_genes <- c("CCL5", "CXCL9", "CXCL10", "CXCL11")
STING_supp_genes <- c("MB21D1", "ENPP1", "TREX1", "IFI16", "IRF3", "ATM", "TBK1", "PARP1", "NT5E", "SLC19A1")

# get them into a list format
cGAS_STING_genes <- list(
  STING = STING_genes,
  NFKB1 = NFKB1_genes,
  NFKB2 = NFKB2_genes,
  STING_act = STING_act_genes,
  STING_supp = STING_supp_genes
)

# Add the genes used by Bakhoum's group in one of their CIN papers
bak_signatures_files <- list.files(BAKHOUM_GENES_DIR, pattern = "CIN_.*Bakhoum.txt", full = T)
bak_signatures <- gsub("\\.txt", "", basename(bak_signatures_files))

bak_genes <- lapply(bak_signatures_files, function(filePath) as.vector(read.table(filePath, header = TRUE)[, 1]))
names(bak_genes) <- bak_signatures
cGAS_STING_genes <- c(cGAS_STING_genes, bak_genes)


# ASSOCIATION BETWEEN CIN SIGNATURE AND GENES -----------------------------
srat = NormalizeData(srat, assay = "RNA")
srat = ScaleData(srat, assay = "RNA")

for (i in seq_len(length(cGAS_STING_genes))){
  srat = get_sign_score(srat, cGAS_STING_genes[[i]], names(cGAS_STING_genes)[i])  
}

random_genes = sample(rownames(srat), 50)
srat = get_sign_score(srat, random_genes, "random")

srat@meta.data %>% 
  ggplot(aes(x = CIN_Signature_Bakhoum, y = STING_supp)) + geom_point() + 
  ggpubr::stat_cor(method = "spearman")

srat@meta.data %>% 
  ggplot(aes(x = CIN_Signature_Bakhoum, y = random)) + geom_point() + 
  ggpubr::stat_cor(method = "spearman")

srat$high_CIN = 
  ifelse(srat$CIN_Signature_Bakhoum > quantile(srat$CIN_Signature_Bakhoum, 0.8), "yes", 
         ifelse(srat$CIN_Signature_Bakhoum < quantile(srat$CIN_Signature_Bakhoum, 0.6), "no", NA))

srat@meta.data %>% filter(!is.na(high_CIN)) %>% 
  ggplot(aes(x = high_CIN, y = CIN_NNKFB_Regulators_Bakhoum)) + 
  geom_boxplot() + ggpubr::stat_compare_means()

srat@meta.data %>% filter(!is.na(high_CIN)) %>% 
  ggplot(aes(x = high_CIN, y = random)) + 
  geom_boxplot() + ggpubr::stat_compare_means()

srat@meta.data %>% filter(!is.na(high_CIN)) %>% 
  ggplot(aes(x = high_CIN, y = ENPP1)) + 
  geom_boxplot() + ggpubr::stat_compare_means()

srat@meta.data %>% 
  ggplot(aes(x = Patient_id, y = random, col = high_CIN)) + 
  geom_boxplot() + ggpubr::stat_compare_means()

srat@meta.data %>% 
  ggplot(aes(x = Patient_id, y = CIN_InflammationGenes_Bakhoum, col = high_CIN)) + 
  geom_boxplot() + ggpubr::stat_compare_means()