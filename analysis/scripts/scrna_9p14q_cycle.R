# Analysis of differences in cell cycle state of 9p and 14q cells
# compared to wild-type cells across different studies


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
SEURAT_PATH = file.path(OUT_DIR, "outputs", "scrna",  "chr3p_subset_qc_exp_filt_merged_clustered.rds")


# FUNCTIONS ---------------------------------------------------------------
plot_violin = function(df, x, y, col){
  p = ggplot(df, aes_string(x = x, y = y)) + 
    ggbeeswarm::geom_beeswarm(alpha = 0.2, cex = 0.1, 
                              aes_string(col = col), 
                              position = position_dodge2(width = .2)) + 
    geom_violin(alpha = .7) + 
    coord_capped_cart(bottom = "none", left = "none") + 
    stat_summary(
      geom = "point",
      fun.y = "mean",
      col = "black",
      size = 3,
      shape = 21,
      fill = "red"
    ) + 
    ggpubr::stat_compare_means() + 
    theme(legend.position = "none")
  return(p)
}

plot_violin_by_pat = function(df, x, y, col, patient_var){
  p = plot_violin(df, x, y, col)
  p = p + facet_wrap(as.formula(paste0("~", patient_var)))
}

source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
srat = readRDS(SEURAT_PATH)
meta = srat@meta.data
# Filter to only sporadic ccRCC
meta = meta[!is.na(meta$Disease_type) & meta$Disease_type == "ccRCC_sporadic", ]
# Binary label combining 9p and 14q loss status
meta = meta %>% 
  mutate(status = case_when(
    Chr9p21.3_Status == "Loss" & Chr14q23_31_Status == "Loss" ~ "9p_14q", 
    Chr9p21.3_Status != "Loss" & Chr14q23_31_Status == "Loss" ~ "14q", 
    Chr9p21.3_Status == "Loss" & Chr14q23_31_Status != "Loss" ~ "9p", 
    TRUE ~ "WT", 
  ))
meta$prolif.score = meta$G2M.Score + meta$S.Score

# ANALYSIS PROLIFERATION BY 9P AND 14Q ------------------------------------

# 9p
p = plot_violin(meta[meta$Chr9p21.3_Status %in% c("Neutral", "Loss"),], 
                "Chr9p21.3_Status", "S.Score", "Patient_id")
save_ggplot(p, file.path(PLOT_DIR, "S_phase_9p_WT"), w = 100, h = 100)
p = plot_violin_by_pat(meta[meta$Chr9p21.3_Status %in% c("Neutral", "Loss"),], 
                "Chr9p21.3_Status", "S.Score", "Patient_id", "Patient_id")
save_ggplot(p, file.path(PLOT_DIR, "S_phase_9p_WT_by_pat"), w = 100, h = 100)

summary(run_lme("prolif.score", "Chr9p21.3_Status", "Patient_id", meta)) # 0.0182 p-value and 2.36 t-value


# 14q
p = plot_violin(meta[meta$Chr14q23_31_Status %in% c("Neutral", "Loss"),], 
                "Chr14q23_31_Status", "S.Score", "Patient_id")
save_ggplot(p, file.path(PLOT_DIR, "S_phase_14q_WT"), w = 100, h = 100)
p = plot_violin_by_pat(meta[meta$Chr14q23_31_Status %in% c("Neutral", "Loss"),], 
                       "Chr14q23_31_Status", "S.Score", "Patient_id", "Patient_id")
save_ggplot(p, file.path(PLOT_DIR, "S_phase_14q_WT_by_pat"), w = 100, h = 100)

summary(run_lme("prolif.score", "Chr14q23_31_Status", "Patient_id", meta)) # 0.0054 p-value and -2.78 t-value
