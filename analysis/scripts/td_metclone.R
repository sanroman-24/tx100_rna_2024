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

purity_dists <- sapply(1:nrow(td_df), function(i) {
  prim_smp <- str_replace(td_df$prim_sample[i], "-", "_")
  met_smp <- str_replace(td_df$met_sample[i], "-", "_")
  annotation$sample <- clean_ids(annotation$sample)
  abs(as.numeric(annotation$purity[annotation$sample == prim_smp]) -
    as.numeric(annotation$purity[annotation$sample == met_smp]))[1]
})

# PLOT ASSOCIATION TRANSCRIPTIONAL AND CLONAL DISTANCE --------------------
td_df$clonal_dists <- clonal_dists
td_df$pd <- unlist(purity_dists)
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

# Check after correcting by purity
summary(lme(
  t_d ~ pd + metastasising_clone,
  random = ~ 1 | Patient,
  data = td_df[!is.na(td_df$pd), ]
))

summary(lme(
  t_d ~ pd + clonal_dists,
  random = ~ 1 | Patient,
  data = td_df[!is.na(td_df$pd), ]
))

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

# Analysis of whether closer clone best reflects ssGSEA
gene_groups <- read_delim("data/meta/martinez_ruiz_2023_hallmark_gs_groups.txt")
ssgsea <- read_delim("data/processed/tx_ssGSEA.tsv")
ssgsea$sample <- clean_ids(ssgsea$sample)
td_df$met_sample <- clean_ids(td_df$met_sample)
td_df$prim_sample <- clean_ids(td_df$prim_sample)

m_met <- match(td_df$met_sample, ssgsea$sample)
m_prim <- match(td_df$prim_sample, ssgsea$sample)

hallmark_signatures <- paste0("HALLMARK_", str_to_upper(hallmark_signatures))

signatures <- colnames(ssgsea) %>%
  {
    .[str_detect(., "motzer")]
  }

df <- data.table::rbindlist(
  lapply(signatures, function(signature) {
    td_df$ssgsea_prim <- ssgsea[[signature]][m_prim]
    td_df$ssgsea_met <- ssgsea[[signature]][m_met]
    td_df$dif_ssgsea <- (td_df$ssgsea_met - td_df$ssgsea_prim)
    return(data.frame(
      signature = signature,
      dif_ssgsea = abs(td_df$dif_ssgsea),
      is_met_clone = td_df$metastasising_clone,
      clonal_dist = td_df$clonal_dists,
      patient = td_df$Patient
    ))
  })
)

df$is_met_clone <- ifelse(df$is_met_clone == "yes_met_clone",
  "Seeding clone", "No seeding clone"
)

df$signature <- str_remove(df$signature, "motzer_") %>%
  str_replace_all("_", " ") %>%
  str_to_title() %>% 
  str_replace("Fas", "FAS") %>% str_replace("Fao Ampk", "FAO-AMPK")

p <- ggplot(df, aes(x = is_met_clone, y = dif_ssgsea)) +
  geom_violin(width = 0.5) +
  geom_point(position = position_dodge2(width = 0.2), alpha = .2) + 
  geom_boxplot(width = 0.1, outlier.shape =  NA) +
  ggpubr::stat_compare_means(size = 2) +
  facet_wrap(~signature) +
  theme(
    strip.text = element_text(size = 7) # Change the facet title size
  ) + 
  labs(y = "Difference in ssGSEA score", x = "Seeding primary")

save_ggplot(p,
  file.path(FIG_DIR, "ssgsea_scores_seeding_vs_no"),
  w = 125, h = 125
)


p <- ggplot(df, aes(x = is_met_clone, y = dif_ssgsea)) +
 geom_violin(width = 0.5) +
  geom_point(position = position_dodge2(width = 0.2), alpha = .2) + 
  geom_boxplot(width = 0.1, outlier.shape =  NA) +
  ggpubr::stat_compare_means(size = 2) +
  labs(y = "Difference in ssGSEA score", x = "Seeding primary")

save_ggplot(p,
  file.path(FIG_DIR, "ssgsea_scores_seeding_vs_no_all"),
  w = 50, h = 50
)

df <- as.data.frame(df)
df$met_clone <- as.factor(df$is_met_clone)
summary(run_lme("dif_ssgsea", "is_met_clone", "patient", df))
