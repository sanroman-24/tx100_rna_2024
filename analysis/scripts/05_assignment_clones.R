# Assignment of gene expression profiles, TME and purity
# to TRACERx Renal clones

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(DESeq2)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
DATA_DIR <- file.path(BASE, "data", "processed")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
VST_PATH <- file.path(DATA_DIR, "tum_vst_exp_filtered.rds")
CONSENSUSTME_PATH = file.path(DATA_DIR, "tum_consensustme.rds")
SSGSEA_PATH <- file.path(DATA_DIR, "tx_ssGSEA.tsv")
CLONES_PATH <- file.path(META_DIR, "clones_annotation.txt")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "get_clonal_dist.R"))
source(file.path(BASE, "src", "assign_to_clones.R"))

# LOAD DATA ---------------------------------------------------------------

sample_annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
vst <- readRDS(VST_PATH)
tme <- readRDS(CONSENSUSTME_PATH)
ssgsea <- read_delim(SSGSEA_PATH, delim = "\t")
clones_annotation <- read_delim(CLONES_PATH, delim = "\t")
clones_annotation <- clones_annotation[!is.na(clones_annotation$parent_clone), ]

# FIND MONOCLONAL REGIONS -------------------------------------------------
clone_regions <- get_monoclonal_regions_df(clones_annotation)
clone_regions$patient <- str_remove(clone_regions$clone, "-.*$")
write_delim(clone_regions,
    file.path(DATA_DIR, "monoclonal_regions_per_clone.tsv"),
    delim = "\t"
)


# ASSIGN EXPRESSION TO EACH CLONE -----------------------------------------

colnames(vst) <- str_replace(colnames(vst), "-", "_")
colnames(vst) <- ifelse(colnames(vst) == "K328_T1", "K328_THR1",
    ifelse(colnames(vst) == "K245_T1", "K245_THR1", colnames(vst))
)

clone_regions <- clone_regions[clone_regions$smp %in% colnames(vst), ]

vst_per_clone <- assign_expression_to_clones(vst, clone_regions)
saveRDS(
    vst_per_clone,
    file.path(DATA_DIR, "vst_per_clone.rds")
)

# ASSIGN PURITY TO EACH CLONE ---------------------------------------------

sample_annotation <- sample_annotation[sample_annotation$sample != "K207_R4", ]
sample_annotation$sample <- str_replace(sample_annotation$sample, "-", "_")
sample_annotation$sample <- ifelse(sample_annotation$sample == "K328_T1",
    "K328_THR1",
    ifelse(sample_annotation$sample == "K245_T1",
        "K245_THR1",
        sample_annotation$sample
    )
)
sample_annotation$purity <- as.numeric(sample_annotation$purity)

pur_per_clone <- assign_purity_to_clones(sample_annotation, clone_regions)

saveRDS(pur_per_clone, file.path(DATA_DIR, "purity_per_clone.rds"))

# ASSIGN TME VALUES FOR EACH CLONE ----------------------------------------
colnames(tme) <- str_replace(colnames(tme), "-", "_")
colnames(tme) <- ifelse(colnames(tme) == "K328_T1", "K328_THR1",
    ifelse(colnames(tme) == "K245_T1", "K245_THR1", colnames(tme))
)

clone_regions <- clone_regions[clone_regions$smp %in% colnames(tme), ]

consensustme_per_clone <- assign_expression_to_clones(tme, clone_regions)

saveRDS(
    consensustme_per_clone,
    file.path(DATA_DIR, "consensustme_per_clone.rds")
)


# ASSIGN SSGSEA TO EACH CLONE ---------------------------------------------

ssgsea_mat <- t(ssgsea[, -1])
colnames(ssgsea_mat) <- ssgsea$sample
colnames(ssgsea_mat) <- str_replace(colnames(ssgsea_mat), "-", "_")
colnames(ssgsea_mat) <- ifelse(colnames(ssgsea_mat) == "K328_T1",
    "K328_THR1",
    ifelse(colnames(ssgsea_mat) == "K245_T1",
        "K245_THR1",
        colnames(ssgsea_mat)
    )
)
ssgsea_mat <- ssgsea_mat[, colnames(ssgsea_mat) %in% clone_regions$smp]

ssgsea_per_clone <- assign_expression_to_clones(ssgsea_mat, clone_regions)

saveRDS(
    ssgsea_per_clone,
    file.path(DATA_DIR, "ssGSEA_per_clone.rds")
)
