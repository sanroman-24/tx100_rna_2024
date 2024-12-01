# Extend TCR and BCR differences to other genomic metrics (reviewer #1)

rm(list = ls(all = TRUE))

# PACKAGES
library(here)
library(DESeq2)
library(tidyverse)
library(nlme)
library(GenomicRanges)
library(data.table)
library(biomaRt)

# PATHS
BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
TCR_BCR_SIM_PATH <- file.path(BASE, "analysis", "outputs", "summary_tcr_bcr_sim.rds")
CNDIST_PATH <- file.path(META_DIR, "cndist_pairs_samps.tsv")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "src", "scna_ith.R"))

clean_cn_ids <- function(cndist_pairs, sample_variable) {
    cndist_pairs[[sample_variable]] <- sapply(
        1:nrow(cndist_pairs),
        function(i) {
            cnid2smp(
                cndist_pairs$pt[i],
                substr(cndist_pairs[[sample_variable]][i], 8, 40)
            )
        }
    )
    return(cndist_pairs)
}


# LOAD DATA ---------------------------------------------------------------
tbcr_similarity <- readRDS(TCR_BCR_SIM_PATH)
annotation <- read_delim(ANNOTATION_PATH)
cndist_pairs <- read_delim(CNDIST_PATH)
subclonal_scna <- get_pat_subSCNA(cndist_pairs)
annotation$sample <- clean_ids(annotation$sample)
annotation$patient <- annotation$Patient

# Add genetic ITH and CNH to the TCR and BCR similarity
genetic_ith <- annotation %>%
    mutate(patient = Patient) %>%
    group_by(patient) %>%
    dplyr::summarise(ITH = median(ITH), wgii = median(wgii)) %>%
    dplyr::select(patient, ITH, wgii)

tbcr_similarity <- left_join(tbcr_similarity, genetic_ith, by = "patient") %>%
    left_join(subclonal_scna)

tbcr_similarity$ith_scale <- log(tbcr_similarity$ITH + 1)

metrics <- c("subclonal_scna", "ith_scale", "wgii")
names(metrics) <- c("% Genome with subclonal CNA", "Genetic ITH", "wGII")
repertoires <- c("tcr_median_sim", "bcr_median_sim")
names(repertoires) <- c("TCR similarity", "BCR similarity")

e <- expand.grid(metrics, repertoires, stringsAsFactors = F)

lp <- lapply(1:nrow(e), function(i) {
    metric <- e[i, 1]
    repertoire <- e[i, 2]
    y_lab <- names(metrics)[metrics == metric]
    x_lab <- names(repertoires)[repertoires == repertoire]
    p <- scatter_cor(tbcr_similarity,
        repertoire, metric,
        y_lab, x_lab,
        add_smooth = T
    )
})

save_plist(lp, file.path(FIG_DIR, "bcr_tcr_genetic_correlates"), w = 150, h = 150, ncol = 3)

