# Association between primary-normal transcriptional distance and number of clones from MRCA


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
NORMAL_PRIM_TD_PATH <- file.path(OUT_DIR, "td_primprim_primnorm.tsv")
CLONES_PATH <- file.path(META_DIR, "clones_annotation.txt")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_clonal_dist.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH)
td_df <- read_delim(NORMAL_PRIM_TD_PATH, delim = "\t")
rownames(td_df) <- colnames(td_df)
clones_annotation <- read_delim(CLONES_PATH, delim = "\t")
clones_annotation <- clones_annotation[!is.na(clones_annotation$parent_clone), ]
clones_annotation$clone_code <- paste0(
    clones_annotation$patient,
    "-",
    clones_annotation$clone
)

# subset only to primary-normal pairs
td_df <- td_df[str_detect(td_df$sample1, "N1t1") |
    str_detect(td_df$sample2, "N1t1"), ]


td_df$normal_sample <- ifelse(str_detect(td_df$sample1, "N1t1"),
    td_df$sample1, td_df$sample2
)
td_df$sample <- ifelse(str_detect(td_df$sample1, "N1t1"),
    td_df$sample2, td_df$sample1
)

td_df <- td_df[, -c(1, 2)]
td_df <- merge(td_df, annotation, by = "sample")
td_df$sample <- str_replace(td_df$sample, "-", "_")

# Estimate the number of edges from MRCA to current clone
clones_annotation$clonal_age <-
    sapply(seq_len(nrow(clones_annotation)), function(i) {
        pt <- clones_annotation$patient[i]
        clone_id <- clones_annotation$clone_code[i]
        calculate_clonal_age(clones_annotation, pt, clone_id, 0)
    })

clones_annotation <- clones_annotation %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(clonal_age = max(clonal_age))

td_df <- merge(td_df, clones_annotation, by = "sample")

# CORRELATE WITH CLONAL DISTANCE ---------------------------

# Clonal distance
td_df <- td_df[!is.na(td_df$clonal_age), ]

summary(
    run_lme("t_d", "clonal_age", "Patient", td_df)
)
# p-value MRCA to clonal age = 0.69

td_df$clonal_age <- as.character(td_df$clonal_age)

p <- violin_plot(td_df,
    x_str = "clonal_age", y_str = "t_d",
    x_title = "Clonal distance to MRCA",
    y_title = "Transcriptional distance",
    labs_x = min(td_df$clonal_age):max(td_df$clonal_age),
    ylim = c(0, 1)
)

save_ggplot(
    p,
    file.path(FIG_DIR, "SupFig8_td_prim_norm_clonal"),
    w = 75, h = 45
)
