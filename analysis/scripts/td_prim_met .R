# Compare I-TED between primary-primary and primary-met pairs in TRACERx Renal

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(ggthemes)
library(lemon)
library(DESeq2)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")
TOP500_GENES_PATH <- file.path(OUT_DIR, "top500_variable_genes.rds")

# FUNCTIONS ---------------------------------------------------------------

source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "td_by_sample_type.R"))
source(file.path(BASE, "src", "get_ited.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
# remove normal samples (we only want primaries & mets)
annotation <- annotation[annotation$type_collapsed %in% c("PRIMARY", "RENAL_MET", "THROMBUS", "LYMPH_NODE", "METASTASIS"), ]
annotation <- annotation[annotation$sample != "K207-R4", ]
top500_genes <- read_rds(TOP500_GENES_PATH)
vst <- assay(read_rds(VST_PATH))
vst <- vst[rownames(vst) %in% top500_genes, ]
vst <- vst[, colnames(vst) != "K207-R4"]

# RE-CALCULATE I-TED INCLUDING METASTASIS SAMPLES TOO ------------------------
d_mat <- get_dist_matrix(vst, top500_genes, f = "dcor")
write_rds(d_mat, file.path(OUT_DIR, "ITED_matrix_with_mets.rds"))
d_mat <- read_rds(file.path(OUT_DIR, "ITED_matrix_with_mets.rds"))


# GET PRIMARY-PRIMARY AND PRIMARY-METASTASIS DISTS ------------------------

met_label <- "METASTASIS"
td_pairs <- get_ited_by_tissue(d_mat, annotation, "PRIMARY", met_label)
p <- plot_paired_boxplot(td_pairs,
    cond1 = "ITH_prim_prim", cond2 = "ITH_prim_noprim",
    ylab = "I-TED", xlab = "", ylim = c(0, 1)
)

save_ggplot(p, file.path(FIG_DIR, "Fig3B_td_primprim_primmet"), w = 40, h = 40)

all_pairs_prim_met <- get_all_td_pairs_prim_noprim(
    d_mat, annotation, 
    prim_labs = "PRIMARY", 
    noprim_labs = c("PRIMARY", "RENAL_MET", "THROMBUS", "LYMPH_NODE", "METASTASIS")
    )

write_delim(all_pairs_prim_met, file.path(OUT_DIR, "td_primprim_primmet.tsv"), delim = "\t")
