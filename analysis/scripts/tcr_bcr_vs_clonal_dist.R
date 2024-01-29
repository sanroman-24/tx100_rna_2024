# Association between TCR and BCR similarity and clonal distance

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(immunarch)
library(here)
library(tidyverse)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
OUT_DIR <- file.path(BASE, "analysis", "outputs")
PLOT_DIR <- file.path(BASE, "analysis", "figures")
TCR_SIM_PATH <- file.path(OUT_DIR, "tcr_repertoire_similarity.RDS")
BCR_SIM_PATH <- file.path(OUT_DIR, "bcr_repertoire_similarity.RDS")
CLONES_PATH <- file.path(META_DIR, "clones_annotation.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_clonal_dist.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
tcr_sim_mt <- readRDS(TCR_SIM_PATH)
bcr_sim_mt <- readRDS(BCR_SIM_PATH)
clones_annotation <- read_delim(CLONES_PATH, delim = "\t")
clones_annotation <- clones_annotation[!is.na(clones_annotation$parent_clone), ]
monoclonal_regions <- get_monoclonal_regions_df(clones_annotation)
clones_annotation <- clones_annotation[
    clones_annotation$sample %in% monoclonal_regions$smp,
]

clonal_dist_df <- data.frame(
    pat = c(), smp1 = c(), smp2 = c(), clonal_dist = c(),
    tcr_sim_mt = c(), bcr_sim_mt = c()
)

for (pat in unique(clones_annotation$patient)) {
    smps <- unique(clones_annotation[clones_annotation$patient == pat, ]$sample)
    if (length(smps) < 2) {
        next()
    }
    pairs <- combn(smps, 2)
    for (pair in 1:ncol(pairs)) {
        s1 <- pairs[1, pair]
        s2 <- pairs[2, pair]

        clonal_dist <- get_dist(s1, s2, clones_annotation)
        if (s1 %in% rownames(tcr_sim_mt) & s2 %in% rownames(tcr_sim_mt)) {
            tcr_sim <- tcr_sim_mt[rownames(tcr_sim_mt) == s1, colnames(tcr_sim_mt) == s2]
        } else {
            tcr_sim <- NA
        }
        if (s1 %in% rownames(bcr_sim_mt) & s2 %in% rownames(bcr_sim_mt)) {
            bcr_sim <- bcr_sim_mt[rownames(bcr_sim_mt) == s1, colnames(bcr_sim_mt) == s2]
        } else {
            bcr_sim <- NA
        }
        clonal_dist_df <- rbind(
            clonal_dist_df,
            data.frame(
                pat = pat, smp1 = s1, smp2 = s2,
                clonal_dist = clonal_dist,
                tcr_sim = tcr_sim,
                bcr_sim = bcr_sim
            )
        )
    }
}

# PLOT ASSOCIATION CLONAL DISTANCE TCR/BCR SIMILARITY -------------------------

clonal_dist_df$dist <- as.character(clonal_dist_df$clonal_dist)

# For TCR
p <- violin_plot(
    df = clonal_dist_df, x_str = "dist",
    y_str = "tcr_sim", x_title = "Clonal distance",
    y_title = "TCR similarity",
    labs_x = as.character(0:6),
    fun.y = "mean",
    ylim = c(0, 1)
)

save_ggplot(p, file.path(PLOT_DIR, "Fig6D_TCR_similarity"), w = 75, h = 45)

summary(run_lme(
    "tcr_sim", "clonal_dist", "pat",
    clonal_dist_df[!is.na(clonal_dist_df$tcr_sim), ]
))
# p = 7e-04 negative association TCR similarity and clonal dist

# For BCR
p <- violin_plot(
    df = clonal_dist_df, x_str = "dist",
    y_str = "bcr_sim", x_title = "Clonal distance",
    y_title = "BCR similarity",
    labs_x = as.character(0:6),
    fun.y = "mean",
    ylim = c(0, 1)
)

save_ggplot(p, file.path(PLOT_DIR, "SuppFigX_BCR_similarity"), w = 75, h = 45)

summary(run_lme(
    "bcr_sim", "clonal_dist", "pat",
    clonal_dist_df[!is.na(clonal_dist_df$bcr_sim), ]
))
# p = 0.21 negative association BCR similarity and clonal dist
