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
VST_PATH <- file.path(BASE, "data", "processed", "tum_nor_vst_exp_filtered.rds")
TOP500_GENES_PATH <- file.path(OUT_DIR, "top500_variable_genes.rds")

# FUNCTIONS ---------------------------------------------------------------

source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "td_by_sample_type.R"))
source(file.path(BASE, "src", "get_ited.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
top500_genes <- read_rds(TOP500_GENES_PATH)
vst <- assay(read_rds(VST_PATH))
vst <- vst[rownames(vst) %in% top500_genes, ]
normals <- str_detect(colnames(vst), "N1t1")
colnames(vst) <- str_remove(colnames(vst), "^[A-Z]_")
normal_ids <- colnames(vst)[normals]
pats_wth_normals <- str_remove(normal_ids, "_N1t1")
vst <- vst[, str_remove(colnames(vst), "[-_].*$") %in% pats_wth_normals]
top500_genes <- read_rds(TOP500_GENES_PATH)
vst <- vst[rownames(vst) %in% top500_genes, ]


# Only primary or normal samples
prims <- annotation$sample[annotation$type_collapsed == "PRIMARY"]
vst <- vst[, c(colnames(vst[!normals]) %in% prims, normals[normals])]


# RE-CALCULATE I-TED INCLUDING PRIMARY AND NORMALS ------------------------
d_mat <- get_dist_matrix(vst, top500_genes, f = "dcor")
write_rds(d_mat, file.path(OUT_DIR, "ITED_matrix_with_normals.rds"))
d_mat <- read_rds(file.path(OUT_DIR, "ITED_matrix_with_normals.rds"))


# GET PRIMARY-PRIMARY AND PRIMARY-NORMAL DISTS ------------------------

annotation <- annotation %>% dplyr::select(Patient, sample, type_collapsed)

annotation <- rbind(annotation, data.frame(
    Patient = pats_wth_normals,
    sample = normal_ids, type_collapsed = "NORMAL"
))

normal_label = "NORMAL"
td_pairs <- get_ited_by_tissue(d_mat, annotation, "PRIMARY", normal_label)
p <- plot_paired_boxplot(td_pairs,
    cond1 = "ITH_prim_prim", cond2 = "ITH_prim_noprim",
    ylab = "I-TED", xlab = "", ylim = c(0, 1)
)

save_ggplot(p, file.path(FIG_DIR, "Fig3A_td_primprim_primnorm"), w = 40, h = 40)

all_pairs_prim_norm <- get_all_td_pairs_prim_noprim(
    d_mat, annotation,
    prim_labs = "PRIMARY",
    noprim_labs = c("NORMAL")
)

write_delim(all_pairs_prim_norm, file.path(OUT_DIR, "td_primprim_primnorm.tsv"), delim = "\t")
