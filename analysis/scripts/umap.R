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
    "K390" = tx_palette[["darkblue"]],
    "K243" = tx_palette[["darkred"]],
    "K153" = tx_palette[["darkpurple"]]
  ))
p <- change_axes(p)

save_ggplot(p, file.path(OUT_DIR, "Fig1a_UMAP"), w = 40, h = 40)