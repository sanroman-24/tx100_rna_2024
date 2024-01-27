################################################################################
################ NORMALISE COUNTS USING VST METHOD #############################
################################################################################

rm(list = ls(all = TRUE))

# PACKAGES ----------------------------------------------------------------
library(DESeq2)
library(tidyverse)
library(here)


# FUNCTIONS ---------------------------------------------------------------

# Returns normalised vst counts from raw counts
#' @param counts gene expression matrix (genes x samples)
#' @param annotation dataframe with sample metadata
#' @param string_formula formula used to create DESeq2 object; doesn't influence VST normalization
vst_norm <- function(counts, annotation, string_formula) {
  dds <- DESeqDataSetFromMatrix(counts, annotation, as.formula(string_formula))
  vsd <- vst(dds, blind = TRUE)
}

# Returns genes that pass expression filter
# expression higher than `tpm_count` in more than `freq`  of patients
#' @param tpms gene expression matrix - TPM - (genes x samples)
#' @param tpm_count minimum expression of gene per sample to pass filter
#' @param freq relative frequency of samples that have to pass tpm_count filter
e_filter <- function(tpms, tpm_count = 1, freq = 0.2) {
  return(rowSums(tpms > tpm_count) > (0.2 * ncol(tpms)))
}

# PATHS -------------------------------------------------------------------
BASE <- here::here()
IN_DIR <- file.path(BASE, "data", "raw")
OUT_DIR <- file.path(BASE, "data", "processed")
TUMOUR_COUNTS_PATH <- file.path(IN_DIR, "counts_tumour.RDS")
TUMOUR_NORMAL_COUNTS_PATH <- file.path(BASE, "data", "raw", "counts_tumour_normal.RDS")
TUM_TPM_PATH <- file.path(OUT_DIR, "tumour_tpm.rds")
TUM_NOR_TPM_PATH <- file.path(OUT_DIR, "tumour_normal_tpm.rds")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")

# LOAD DATA ---------------------------------------------------------------

tum_tpm <- readRDS(TUM_TPM_PATH)
tum_nor_tpm <- readRDS(TUM_NOR_TPM_PATH)
tum_counts <- readRDS(TUMOUR_COUNTS_PATH)
tum_nor_counts <- readRDS(TUMOUR_NORMAL_COUNTS_PATH)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")


# APPLY VST NORMALIZATION -------------------------------------------------
tum_vst <- vst_norm(tum_counts, annotation, "~VHL_type")
saveRDS(tum_vst, file.path(OUT_DIR, "tum_vst.rds"))

# only in genes that pass expression filter
is_gt <- e_filter(tum_tpm)
# sum(is_gt)
# length(is_gt)
# 16716 out of 23299 genes pass expression filter
tum_vst <- vst_norm(tum_counts[is_gt, ], annotation, "~VHL_type")
saveRDS(tum_vst, file.path(OUT_DIR, "tum_vst_exp_filtered.rds"))

# add to the tumour sample annotation information of the normal samples
is_normal <- grepl("N1t1$", colnames(tum_nor_counts))
annotation <- data.frame(
  sample = c(annotation$sample, colnames(tum_nor_counts)[is_normal]),
  type = c(
    rep("tumor", sum(!is_normal)),
    rep("normal", sum(is_normal))
  )
)

tum_nor_vst <- vst_norm(tum_nor_counts, annotation, "~type")
saveRDS(tum_nor_vst, file.path(OUT_DIR, "tum_nor_vst.rds"))

# only in genes that pass expression filter
is_gt <- e_filter(tum_nor_tpm)
# sum(is_gt)
# length(is_gt)
# 16716 out of 23299 genes pass expression filter
tum_nor_vst <- vst_norm(tum_nor_counts[is_gt, ], annotation, "~type")
saveRDS(tum_nor_vst, file.path(OUT_DIR, "tum_nor_vst_exp_filtered.rds"))
