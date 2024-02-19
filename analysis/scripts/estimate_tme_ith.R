# Estimate transcriptional intratumour heterogeneity in primary ccRCC

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(DESeq2)
library(energy)
library(pheatmap)
library(ggthemes)
library(lsa)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
PLOT_DIR <- file.path(BASE, "analysis", "figures")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
TME_PATH <- file.path(BASE, "data", "processed", "tum_consensustme.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))


# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
tme <- read_rds(TME_PATH)
annotation$ITH <- as.numeric(annotation$ITH)
annotation$wgii <- as.numeric(annotation$wgii)
annotation$Regions <- as.numeric(annotation$Regions)

annotation <- annotation[annotation$sample != "K207-R4", ]
tme <- tme[, colnames(tme) != "K207-R4"]

# check that it is good to go
all(annotation$sample == colnames(tme))

# CALCULATE TRANSCRIPTIONAL DISTANCE --------------------------------------
d_mat <- get_dist_matrix(tme, rownames(tme), f = "cosine")
write_rds(d_mat, file.path(OUT_DIR, "TME_dist_matrix.rds"))

d_mat = read_rds(file.path(OUT_DIR, "TME_dist_matrix.rds"))
# FILTER TO PATIENTS WITH MORE THAN 1 PRIMARY REGION SAMPLED

annotation <- annotation %>%
    filter(type_collapsed == "PRIMARY")

# get samples with more than one tumour primary regions
keep_patients <- annotation %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise(n = n()) %>% # count number of different samples of each patient
    filter(n > 1) %>%
    dplyr::select(Patient) %>%
    as_vector() # patients with +1 primary region

keep_samples <- annotation[annotation$Patient %in% keep_patients, "sample"] %>%
    as_vector()

annotation <- annotation[annotation$sample %in% keep_samples, ]
tme <- tme[, colnames(tme) %in% keep_samples]

d_mat <- get_dist_matrix(tme, rownames(tme), f = "cosine")
write_rds(d_mat, file.path(OUT_DIR, "TME_dist_primary.rds"))

d_mat = read_rds(file.path(OUT_DIR, "TME_dist_primary.rds"))
# p = pheatmap(d_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
# png(file.path(OUTPUT_PLOTS_DIR, "transcriptional_distance_matrix.png"), height = 100, width = 100, units = "mm", res = 300)
# p
# dev.off()


# CALCULATE I-TED ---------------------------------------------------------
summary_tme_ith <- data.frame(patient = keep_patients)

summary_tme_ith$min_tme_ith <- summarise_ited(d_mat, summary_tme_ith$patient, "min")
summary_tme_ith$max_tme_ith <- summarise_ited(d_mat, summary_tme_ith$patient, "max")
summary_tme_ith$median_tme_ith <- summarise_ited(d_mat, summary_tme_ith$patient, "median")

tme_ith_pairs <- unlist(summarise_ited(d_mat, unname(keep_patients), "c"))
tme_ith_pairs_df <- data.frame(pair_id = names(tme_ith_pairs), t_d = unname(tme_ith_pairs))
tme_ith_pairs_df$patient <- str_sub(tme_ith_pairs_df$pair_id, 1, 4)


# PLOT I-TED --------------------------------------------------------------
tme_ith_pairs_df$patient <- factor(tme_ith_pairs_df$patient)
tme_ith_pairs_df <- tme_ith_pairs_df %>%
    mutate(patient_ord = fct_reorder(patient, t_d))

p <- ggplot(tme_ith_pairs_df, aes(x = patient_ord, y = t_d)) +
    geom_point(alpha = 0.3, col = "lightblue") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p <- p +
    geom_point(
        data = summary_tme_ith,
        aes(x = patient, y = max_tme_ith), shape = 15, col = "orchid", alpha = 0.5
    ) +
    geom_point(
        data = summary_tme_ith,
        aes(x = patient, y = min_tme_ith), shape = 15, col = "orchid", alpha = 0.5
    ) +
    geom_point(
        data = summary_tme_ith,
        aes(x = patient, y = median_tme_ith), shape = 23, col = "black", fill = "blueviolet"
    ) +
    labs(x = "", y = "Transcriptional distance") +
    theme(axis.text.x = element_text(size = 8))

p <- change_axes(p)

save_ggplot(p, file.path(PLOT_DIR, "SuppFigX_tme_ith_primary"), w = 180, h = 70)
