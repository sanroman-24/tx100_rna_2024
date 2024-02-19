# Compare ssGSEA scores of early/late clones in TRACERx Renal

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(DESeq2)
library(ggpubr)
library(nlme)
library(ggthemes)
library(gridExtra)

# PATHS -------------------------------------------------------------------


BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
TME_PATH <- file.path(BASE, "data", "processed", "consensustme_per_clone.rds")
CLONES_PATH <- file.path(META_DIR, "clones_annotation.txt")
HALLMARK_GROUPS_PATH <- file.path(META_DIR, "martinez_ruiz_2023_hallmark_gs_groups.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "get_clonal_dist.R"))
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------

gene_groups <- read_delim(HALLMARK_GROUPS_PATH, delim = "\t")
annotation <- read_delim(ANNOTATION_PATH)
tme <- readRDS(TME_PATH)
cells <- rownames(tme)
clones_annotation <- read_delim(CLONES_PATH, delim = "\t")
clones_annotation <- clones_annotation[!is.na(clones_annotation$clone), ]


# ANNOTATE TERMINAL CLONES ------------------------------------------------

clones_annotation$pat_terminal <- NA
clones_annotation$clone_code <- paste0(
    clones_annotation$patient,
    "-",
    clones_annotation$clone
)

# Find terminal clones
for (pat in unique(clones_annotation$patient)) {
    # find terminal (in leaf) clones in the patient
    term_clones <- get_terminal_clones_pt(clones_annotation, pat)
    # annotate accordingly
    clones_annotation[clones_annotation$patient == pat &
        clones_annotation$clone %in% term_clones, ]$pat_terminal <- "Terminal"
    clones_annotation[clones_annotation$patient == pat &
        !clones_annotation$clone %in% term_clones, ]$pat_terminal <- "No terminal"
}

# Calculate 'clonal age' or number of edges connecting root to leaf
clones_annotation$clonal_age <-
    sapply(seq_len(nrow(clones_annotation)), function(i) {
        pt <- clones_annotation$patient[i]
        clone_id <- clones_annotation$clone_code[i]
        calculate_clonal_age(clones_annotation, pt, clone_id, 0)
    })

# COMPARE ALL THE POPULATIONS -------------------------------------------

tme <- t(tme) %>%
    as.data.frame() %>%
    rownames_to_column("clone_code")
tme$clone_code[83] <- "K207-C4"
tme <- left_join(tme, clones_annotation, by = "clone_code") %>%
    {
        .[!duplicated(.$clone_code), ]
    }

tme <- tme[!is.na(tme$pat_terminal), ]

tme <- merge(tme, annotation, by = "sample")

# make all the paired boxplots
plist <- lapply(
    cells,
    function(cell) {
        df_to_ggpaired(
            df = tme, pat_var = "patient",
            cond = "pat_terminal", var_str = cell
        ) %>%
            plot_paired_boxplot(
                cond1 = "No terminal", cond2 = "Terminal",
                ylab = str_replace_all(cell, "_", " "),
                xlab = "", ylim = c(NA)
            )
    }
)

save_plist(plist,
    fp = file.path(FIG_DIR, "SupFig14_TME_term_noterm"),
    w = 275, h = 275, ncol = 4
)

# Compare also differences in purity
p <- df_to_ggpaired(
    df = tme, pat_var = "patient",
    cond = "pat_terminal", var_str = "purity"
) %>%
    plot_paired_boxplot(
        cond1 = "No terminal", cond2 = "Terminal",
        ylab = "Tumour Purity",
        xlab = "", ylim = c(NA)
    )

save_ggplot(p,
    file.path(FIG_DIR, "SupFig15_pur_term_noterm"),
    w = 75, h = 75
)
