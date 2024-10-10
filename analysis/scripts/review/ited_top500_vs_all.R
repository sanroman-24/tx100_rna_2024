# Compare I-TED with all genes vs to top500


rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(ggthemes)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs", "review")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
ITED500_PATH <- file.path(BASE, "analysis", "outputs", "ITED_matrix_primary.rds")
ITEDA_PATH <- file.path(BASE, "analysis", "outputs", "review", "ITED_matrix_ngenes", "ITED_matrix_primary_allgenes.rds")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))

load_ited_ngenes <- function(ngenes, pts) {
    d_mat <- readRDS(
        file.path(
            BASE, "analysis", "outputs", "review", "ITED_matrix_ngenes",
            glue::glue("ITED_matrix_primary_ngenes_{ngenes}.rds")
        )
    )
    return(tdmat2ited(d_mat, pts))
}

tdmat2ited <- function(d_mat, pts) {
    summary_ited <- data.frame(patient = keep_patients)

    summary_ited$min_ited <- summarise_ited(d_mat, summary_ited$patient, "min")
    summary_ited$max_ited <- summarise_ited(d_mat, summary_ited$patient, "max")
    summary_ited$median_ited <- summarise_ited(d_mat, summary_ited$patient, "median")

    ited_pairs <- unlist(summarise_ited(d_mat, unname(keep_patients), "c"))
    ited_pairs_df <- data.frame(pair_id = names(ited_pairs), t_d = unname(ited_pairs))
    ited_pairs_df$patient <- str_sub(ited_pairs_df$pair_id, 1, 4)

    return(summary_ited)
}

# LOAD DATA ---------------------------------------------------------------
tdmat500 <- readRDS(ITED500_PATH)
tdmat <- readRDS(ITEDA_PATH)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")

# filter to only primary and more than one tumour primary region
annotation <- annotation %>%
    filter(type_collapsed == "PRIMARY")
annotation <- annotation[annotation$sample != "K207-R4", ]

keep_patients <- annotation %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise(n = n()) %>% # count number of different samples of each patient
    filter(n > 1) %>%
    dplyr::select(Patient) %>%
    as_vector() # patients with +1 primary region
keep_samples <- annotation[annotation$Patient %in% keep_patients, "sample"] %>%
    as_vector()
annotation <- annotation[annotation$sample %in% keep_samples, ]

# get ITED
ited <- tdmat2ited(tdmat, keep_patients)
ited500 <- tdmat2ited(tdmat500, keep_patients)
ited1000 <- load_ited_ngenes(1000, keep_patients)
ited2500 <- load_ited_ngenes(2500, keep_patients)
ited5000 <- load_ited_ngenes(5000, keep_patients)
ited10000 <- load_ited_ngenes(10000, keep_patients)
ited12500 <- load_ited_ngenes(12500, keep_patients)
ited15000 <- load_ited_ngenes(15000, keep_patients)

cor.test(ited500[, 4], ited[, 4], method = "pearson")

crt1000 <- cor.test(ited500[, 4], ited1000[, 4], method = "pearson")
crt2500 <- cor.test(ited500[, 4], ited2500[, 4], method = "pearson")
crt5000 <- cor.test(ited500[, 4], ited5000[, 4], method = "pearson")
crt10000 <- cor.test(ited500[, 4], ited10000[, 4], method = "pearson")
crt12500 <- cor.test(ited500[, 4], ited12500[, 4], method = "pearson")
crt15000 <- cor.test(ited500[, 4], ited15000[, 4], method = "pearson")

p <- data.frame(estimate = c(
    crt1000$estimate, crt2500$estimate, crt5000$estimate,
    crt10000$estimate, crt12500$estimate, crt15000$estimate
), ngenes = c(1000, 2500, 5000, 10000, 12500, 15000)
) %>% 
    ggplot(aes(x = ngenes, y = estimate)) + 
    geom_point(size = 2, pch = 21, fill = "grey60") + 
    lims(y = c(0,1)) + 
    labs(x = "Number genes in I-TED", y = "Pearson's r with I-TED 500")

save_ggplot(p, file.path(FIG_DIR, "pearson_ited500_vs_ngenes"), w = 50, h = 50)

plot_cor <- function(i1, i2, x_axis, y_axis) {
    df <- data.frame(i1 = i1, i2 = i2)
    p <- scatter_plot(
        df = df,
        y_str = "i2", x_str = "i1",
        x_title = x_axis,
        y_title = y_axis
    ) +
        geom_smooth(method = "lm", se = FALSE, col = "red") +
        ggpubr::stat_cor(method = "spearman", size = 2)

    return(p)
}

p1000 <- plot_cor(ited500[, 4], ited1000[, 4], "I-TED 500", "I-TED 1000")
p2500 <- plot_cor(ited500[, 4], ited2500[, 4], "I-TED 500", "I-TED 2500")
p5000 <- plot_cor(ited500[, 4], ited5000[, 4], "I-TED 500", "I-TED 5000")
p10000 <- plot_cor(ited500[, 4], ited10000[, 4], "I-TED 500", "I-TED 10000")
p12500 <- plot_cor(ited500[, 4], ited12500[, 4], "I-TED 500", "I-TED 12500")
p15000 <- plot_cor(ited500[, 4], ited15000[, 4], "I-TED 500", "I-TED 15000")

lp <- list(p1000, p2500, p5000, p10000, p12500, p15000)

save_plist(lp, file.path(FIG_DIR, "ited_scatter_plots_ngenes"), w = 150, h = 150, ncol = 2)
