# Compare intra to interpatient transcriptional heterogeneity

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(ggthemes)

# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "outputs")
FIG_DIR = file.path(BASE, "analysis", "figures")
ITED500_PATH <- file.path(BASE, "analysis", "outputs", "ITED_matrix_primary.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

is_same_pt <- function(d_mat, i,j){
    pt1 = str_extract(rownames(d_mat)[i], "^[^_-]+")
    pt2 = str_extract(colnames(d_mat)[j], "^[^_-]+")
    o <- pt1 == pt2
    names(o) <- paste0(pt1, ":", pt2)
    return(o)
}

# LOAD DATA ---------------------------------------------------------------
d_mat <- readRDS(ITED500_PATH)

td_intra <- c()
td_inter <- c()

for (i in 1:nrow(d_mat)){
    for (j in 1:ncol(d_mat)){
        bl <- is_same_pt(d_mat, i, j)
        nm <- names(bl)
        if (bl){
            td_intra <- c(td_intra, d_mat[i,j])
            names(td_intra)[length(names(td_intra))] <- nm
        } else {
            td_inter <- c(td_inter, nm = d_mat[i,j])
            names(td_inter)[length(names(td_inter))] <- nm
        }
    }
}

# same sample with same sample is NA, so remove
td_intra <- na.omit(td_intra) 

# Plot the differences in distribution

dist_df = rbind(
    data.frame(same_pat = "Intra", dist = td_intra), 
    data.frame(same_pat = "Inter", dist = td_inter)
)

p <- ggplot(dist_df, aes(x = same_pat, y = dist)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.1, alpha = .3, outlier.shape = NA) +
  ggpubr::stat_compare_means(size = 2) +
  labs(x = "Sample pair from same patient", y = "Transcriptional distance")

p <- change_axes(p)

save_ggplot(p, file.path(FIG_DIR, "SupFig2_ith_inter_ratio"), w = 35, h = 35)


# Estimate how much higher ITH compared to intertumour heterogeneity
B <- 1000

b_ratios <- c()
for (b in 1:B){
    b_same <- sample(td_intra, replace = T)
    b_diff <- sample(td_inter, replace = T)
    b_ratios <- c(b_ratios, median(b_diff) / median(b_same))
}

# 95% CI for the mean will be between 2.5% and 97.5% quantiles
quantile(b_ratios, .025)
quantile(b_ratios, 0.975)
mean(b_ratios)
hist(b_ratios)

