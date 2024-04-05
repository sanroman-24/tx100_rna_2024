# Association between clonal distance and transcriptional distance

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(ggthemes)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
MET_PRIM_TD_PATH <- file.path(OUT_DIR, "td_primprim_primmet.tsv")
CLONES_PATH <- file.path(META_DIR, "clones_annotation.txt")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_clonal_dist.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH)
td_df <- read_delim(MET_PRIM_TD_PATH, delim = "\t")
rownames(td_df) <- colnames(td_df)
clones_annotation <- read_delim(CLONES_PATH, delim = "\t")
clones_annotation <- clones_annotation[!is.na(clones_annotation$parent_clone), ]

# subset association to only primary-met pairs
mets <- annotation %>%
  filter(type_collapsed %in% c("RENAL_MET", "THROMBUS", "LYMPH_NODE", "METASTASIS")) %>%
  dplyr::select(sample, Patient, type_collapsed)

# make sure that always one of the sample is a met and the other is primary
td_df <- td_df[td_df$sample1 %in% mets$sample |
  td_df$sample2 %in% mets$sample, ]

td_df <- td_df[!td_df$sample1 %in% mets$sample |
  !td_df$sample2 %in% mets$sample, ]

td_df <- td_df %>%
  mutate(
    prim_sample = ifelse(sample1 %in% mets$sample, sample2, sample1),
    met_sample = ifelse(sample1 %in% mets$sample, sample1, sample2)
  ) %>%
  dplyr::select(met_sample, prim_sample, t_d)
td_df <- merge(td_df, annotation, by.x = "met_sample", by.y = "sample")



# CALCULATE CLONAL DISTANCE BETWEEN SAMPLES -------------------------------

clonal_dists <- sapply(1:nrow(td_df), function(i) {
  prim_smp <- str_replace(td_df$prim_sample[i], "-", "_")
  met_smp <- str_replace(td_df$met_sample[i], "-", "_")
  get_met_prim_dist(met_smp, prim_smp, clones_annotation)
})


# PLOT ASSOCIATION TRANSCRIPTIONAL AND CLONAL DISTANCE --------------------
td_df$clonal_dists <- clonal_dists
td_df$metastasising_clone <- ifelse(clonal_dists == 0, "yes_met_clone", "no_met_clone")

td_df <- td_df[!is.na(td_df$metastasising_clone), ]

p <- violin_plot(
  df = td_df, x_str = "metastasising_clone",
  y_str = "t_d", x_title = "Primary region",
  y_title = "Transcriptional distance",
  labs_x = c("No metastasising clone", "Metastasising clone"),
  fun.y = "mean",
  ylim = c(0, 1)
)

save_ggplot(p, file.path(FIG_DIR, "Fig2D_td_metclone"), w = 40, h = 40)

# Run LME to control for inclusion mulitple pairs from same patient
summary(run_lme(
  "t_d", "metastasising_clone", "Patient",
  td_df
))
# Met clone more similar with p-val = 7e-04

# LME with continuous clonal distance
summary(run_lme(
  "t_d", "clonal_dists", "Patient",
  td_df
))

# general higher distance with higher clonal distance
# p-value = 0

td_df$clonal_dists <- as.character(td_df$clonal_dists)

p <- violin_plot(
  df = td_df,
  x_str = "clonal_dists",
  y_str = "t_d", x_title = "Clonal distance",
  labs_x = seq(
    min(td_df$clonal_dists, na.rm = T),
    max(td_df$clonal_dists, na.rm = T)
  ),
  y_title = "Transcriptional distance",
  fun.y = "mean",
  ylim = c(0, 1)
)

save_ggplot(p, file.path(FIG_DIR, "SupFig10_td_prim_met"),
  w = 75, h = 45
)
