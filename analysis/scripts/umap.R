# Generate UMAP TRACERx Renal samples & compare ratio
# transcriptional intertumour to intratumour heterogeneity

rm(list = ls(all = TRUE))


# PACKAGES ----------------------------------------------------------------
library(umap)
library(DESeq2)
library(tidyverse)
library(data.table)
library(lemon)
library(ggthemes)
library(ggpubr)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "figures")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")


# FUNCTIONS ---------------------------------------------------------------

source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

euclidean_dist <- function(x1, y1, x2, y2) {
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
}

calculate_dist_umap <- function(s1, s2, umap_df) {
  x1 <- umap_df[umap_df$sample == s1, ]$UMAP1
  y1 <- umap_df[umap_df$sample == s1, ]$UMAP2
  x2 <- umap_df[umap_df$sample == s2, ]$UMAP1
  y2 <- umap_df[umap_df$sample == s2, ]$UMAP2
  return(euclidean_dist(x1, y1, x2, y2))
}

# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
vst <- read_rds(VST_PATH)
annotation$ITH <- as.numeric(annotation$ITH)
annotation$wgii <- as.numeric(annotation$wgii)
annotation$Regions <- as.numeric(annotation$Regions)

# K207-R4 is duplicated, with one sample on each batch, so remove one of both samples
annotation <- annotation[annotation$sample != "K207-R4", ]
vst <- vst[, colnames(vst) != "K207-R4"]

# check that it is good to go
all(annotation$sample == colnames(vst))


# FIT UMAP ----------------------------------------------------------------

set.seed(142)
umap_fit <- assay(vst) %>%
  t() %>%
  umap()

umap_df <- umap_fit$layout %>%
  as.data.frame() %>%
  dplyr::rename(
    UMAP1 = "V1",
    UMAP2 = "V2"
  )

umap_df$sample <- colnames(vst)
umap_df <- merge(umap_df, annotation, by = "sample")


# PLOT UMAP ---------------------------------------------------------------
umap_df$highlight <- ifelse(umap_df$Patient == "K243", "K243",
  ifelse(umap_df$Patient == "K153", "K153",
    ifelse(umap_df$Patient == "K390", "K390", "other")
  )
)
p <- umap_df %>%
  ggplot(aes(
    x = UMAP1,
    y = UMAP2,
    col = highlight
  )) +
  geom_point(alpha = 0.6, size = .4) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme(legend.position = "none") +
  scale_color_manual(values = c(
    "other" = "grey60",
    "K390" = "blue",
    "K243" = "red",
    "K153" = "orange"
  ))
p <- change_axes(p)

save_ggplot(p, file.path(OUT_DIR, "Fig1a_UMAP"), w = 40, h = 40)

# ESTIMATE DISTANCE IN UMAP SPACE SAMPLE SAME PATIENT VS DIFFERENT --------
dist_df <- data.frame(s1 = c(), s2 = c(), same_pat = c(), dist = c())

pairs <- combn(umap_df$sample, 2)

dist_df <- data.table::rbindlist(
  lapply(1:ncol(pairs), function(pair) {
    s1 <- pairs[1, pair]
    s2 <- pairs[2, pair]
    dist <- calculate_dist_umap(s1, s2, umap_df)
    p1 <- str_remove(s1, "-.*$")
    p2 <- str_remove(s2, "-.*$")
    same_pat <- ifelse(p1 == p2, "Yes", "No")
    data.frame(s1, s2, same_pat, dist)
  })
)


p <- ggplot(dist_df, aes(x = same_pat, y = dist)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.1, alpha = .3, outlier.shape = NA) +
  ggpubr::stat_compare_means(size = 2) +
  labs(x = "Sample pair from same patient", y = "Distance in UMAP")

p <- change_axes(p)

save_ggplot(p, file.path(OUT_DIR, "SupFig2_ith_inter_ratio_umap"), w = 35, h = 35)

# Estimate how much higher ITH compared to intertumour heterogeneity
B <- 1000
dist_same_pat <- dist_df$dist[dist_df$same_pat == "Yes"]
dist_diff_pat <- dist_df$dist[dist_df$same_pat == "No"]

# Use Bootstrap to get 95% CI
b_ratios <- c()
for (b in 1:B) {
  b_same <- sample(dist_same_pat, replace = T)
  b_diff <- sample(dist_diff_pat, replace = T)
  b_ratios <- c(b_ratios, mean(b_diff) / mean(b_same))
}

# 95% CI for the mean will be between 2.5% and 97.5% quantiles
quantile(b_ratios, .025)
quantile(b_ratios, 0.975)
mean(b_ratios)
hist(b_ratios)

# CHECK IF CLUSTERING MORE OFTEN THAN EXPECTED BY CHANCE --------
dist_df$Patient <- str_remove(dist_df$s1, "[-_][A-Z0-9]+")

closest <- dist_df %>%
  dplyr::group_by(s1) %>%
  slice_min(dist)

ct_tb <- matrix(c(table(closest$same_pat), table(dist_df$same_pat)),
  nrow = 2, byrow = T
)

# Closest sample in UMAP space is from the same patient more often
# than expected --> p-value < 2.2e-16
chisq.test(ct_tb)

# 28.8x more often than expected
o_fq <- ct_tb[1, 2] / sum(ct_tb[1, ])
e_fq <- ct_tb[2, 2] / sum(ct_tb[2, ])
o_fq / e_fq

df <- data.frame(
  same_patient = c("yes", "no", "yes", "no"),
  t = c("observed", "observed", "expected", "expected"),
  fq = c(o_fq, 1 - o_fq, e_fq, 1 - e_fq)
)

p <- ggplot(df, aes(x = t, fill = same_patient, y = fq)) +
  geom_col() +
  scale_fill_manual(values = c("yes" = "#ADD8E6", "no" = "grey80")) +
  coord_capped_flip() +
  theme(legend.position = "none")

save_ggplot(p, file.path(OUT_DIR, "SupFig5b_UMAP_clustering"), w = 60, h = 35)
