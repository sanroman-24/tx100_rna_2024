# Association between clonal distance and TME distance

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(ggthemes)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
FIG_DIR <- file.path(BASE, "analysis", "figures")
TME_MAT_PATH <- file.path(BASE, "analysis", "outputs", "TME_dist_primary.rds")
CLONES_PATH <- file.path(META_DIR, "clones_annotation.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_clonal_dist.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------

tme_mat <- readRDS(TME_MAT_PATH)
colnames(tme_mat) <- str_remove(colnames(tme_mat), "^[RG]_") %>%
    str_remove("_r\\d") %>%
    str_replace("K328R19", "K328-R19") %>%
    str_replace("-", "_")
colnames(tme_mat) <- ifelse(colnames(tme_mat) == "K328_T1", "K328_THR1",
    ifelse(colnames(tme_mat) == "K245_T1", "K245_THR1", colnames(tme_mat))
)
rownames(tme_mat) <- colnames(tme_mat)
clones_annotation <- read_delim(CLONES_PATH, delim = "\t")
clones_annotation <- clones_annotation[!is.na(clones_annotation$parent_clone), ]
monoclonal_regions <- get_monoclonal_regions_df(clones_annotation)
clones_annotation <- clones_annotation[clones_annotation$sample %in% monoclonal_regions$smp, ]

# CALCULATE TRANSC DISTANCE BETWEEN SAMPLES -------------------------------

clonal_dist_df <- data.frame(
    pat = c(), smp1 = c(), smp2 = c(), clonal_dist = c(),
    tme = c()
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
        if (s1 %in% rownames(tme_mat) & s2 %in% rownames(tme_mat)) {
            tme <- tme_mat[rownames(tme_mat) == s1, colnames(tme_mat) == s2]
        } else {
            tme <- NA
        }

        clonal_dist_df <- rbind(
            clonal_dist_df,
            data.frame(
                pat = pat, smp1 = s1, smp2 = s2,
                clonal_dist = clonal_dist,
                tme = tme
            )
        )
    }
}

# PLOT ASSOCIATION TRANSCRIPTIONAL AND CLONAL DISTANCE --------------------
clonal_dist_df$dist <- as.character(clonal_dist_df$clonal_dist)
p <- violin_plot(
    df = clonal_dist_df, x_str = "dist",
    y_str = "tme", x_title = "Clonal distance",
    y_title = "TME distance",
    labs_x = seq(
        min(clonal_dist_df$dist),
        max(clonal_dist_df$dist)
    ),
    fun.y = "mean",
    ylim = c(0, 1)
)

save_ggplot(p, file.path(FIG_DIR, "SupFig13_tme_clonaldist"), w = 75, h = 45)
# Run LME to control for inclusion mulitple pairs from same patient
summary(run_lme(
    "clonal_dist", "tme", "pat",
    clonal_dist_df[!is.na(clonal_dist_df$tme), ]
))

# p-value negative association TME dist clonal dist = 0.39
