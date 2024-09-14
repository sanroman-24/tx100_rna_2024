# Compare ssGSEA scores between matched primary and metastasis samples in TRACERx Renal



rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(nlme)
library(harmonicmeanp)
library(ggthemes)
library(patchwork)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
SSGSEA_NORMALS_PATH <- file.path(BASE, "data", "processed", "tx_normal_ssGSEA.tsv")
SSGSEA_TUMOURS_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
HALLMARK_GROUPS_PATH <- file.path(BASE, "data", "meta", "martinez_ruiz_2023_hallmark_gs_groups.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))


# LOAD DATA ---------------------------------------------------------------
gene_groups <- read_delim(HALLMARK_GROUPS_PATH, delim = "\t")
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
ssgsea_tumour <- read_delim(SSGSEA_TUMOURS_PATH)
ssgsea_tumour$patient <- str_remove(ssgsea_tumour$sample, "[-_].*$")

# subset ssgsea_tumour to only patients with matched primary normal sample
prim_smps <- annotation[annotation$type_collapsed == "PRIMARY", ]$sample
ssgsea_primary <- ssgsea_tumour[ssgsea_tumour$sample %in% prim_smps, ]
met_smps <- annotation[annotation$type_collapsed == "METASTASIS", ]$sample
ssgsea_met <- ssgsea_tumour[ssgsea_tumour$sample %in% met_smps, ]

keep_pats <- dplyr::intersect(ssgsea_primary$patient, ssgsea_met$patient)
ssgsea_primary <- ssgsea_primary[ssgsea_primary$patient %in% keep_pats, ]
ssgsea_met <- ssgsea_met[ssgsea_met$patient %in% keep_pats, ]

ssgsea_primary$sample_type <- "primary_tumour"
ssgsea_met$sample_type <- "met"
ssgsea_df <- rbind(ssgsea_met, ssgsea_primary)

# LME COMPARING SAMPLE TYPES ----------------------------------------------
hm_sigs <- colnames(ssgsea_df)[str_detect(colnames(ssgsea_df), "HALLMARK")]
ssgsea_difs_df <- data.frame(sig = c(), t_value = c(), p_value = c())
for (sig in hm_sigs) {
  res <- coef(summary(run_lme(sig, "sample_type", "patient", ssgsea_df)))
  ssgsea_difs_df <- rbind(ssgsea_difs_df, data.frame(sig = sig, t_value = res[2, 4], p_value = res[2, 5]))
}

ssgsea_difs_df$p_adj <- p.adjust(ssgsea_difs_df$p_value)

# collapse by gene groups
ssgsea_difs_df$Hallmark <- str_remove(str_to_lower(ssgsea_difs_df$sig), "hallmark_")
ssgsea_difs_df <- merge(ssgsea_difs_df, gene_groups, by = "Hallmark")

write_delim(ssgsea_difs_df, file.path(OUT_DIR, "ST2_ssgsea_difs_paired_primary_met.tsv"), delim = "\t")
write_delim(ssgsea_difs_df, file.path(OUT_DIR, "ST2_ssgsea_difs_paired_primary_met.csv"), delim = ",")
