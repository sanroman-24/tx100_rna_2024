################################################################################
################ NORMALISE RAW COUNTS TO TPM ###################################
################################################################################

rm(list = ls(all = TRUE))


# PACKAGES ----------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(here)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
IN_DIR <- file.path(BASE, "data", "raw")
OUT_DIR <- file.path(BASE, "data", "processed")
TUMOUR_COUNTS_PATH <- file.path(IN_DIR, "counts_tumour.RDS")
TUMOUR_NORMAL_COUNTS_PATH <- file.path(BASE, "data", "raw", "counts_tumour_normal.RDS")
TXI_TUMOUR_PATH <- file.path(IN_DIR, "txi_tumour.rds")
TXI_TUMOUR_NORMAL_PATH <- file.path(IN_DIR, "txi_tumour_normal.rds")

# FUNCTIONS ---------------------------------------------------------------

#' @param counts gene expression matrix (genes x samples)
#' @param length_mat average gene length information (based on transcript isoform
#'  usage) obtained from RSEM output after using tximport package
compute_tpm <- function(counts, length_mat) {
  rpk <- counts / length_mat
  per_mill_scalar <- colSums(rpk) / 1e6
  tpms <- sweep(rpk, 2, per_mill_scalar, "/")
  return(tpms)
}

# LOAD DATA ---------------------------------------------------------------

tum_counts <- readRDS(TUMOUR_COUNTS_PATH)
tum_nor_counts <- readRDS(TUMOUR_NORMAL_COUNTS_PATH)
txi_tum <- readRDS(TXI_TUMOUR_PATH)
txi_tum_nor <- readRDS(TXI_TUMOUR_NORMAL_PATH)

# OBTAIN tpms -------------------------------------------------------------

tum_tpm <- computetpm(tum_counts, txi_tum$length)
saveRDS(tum_tpm, file.path(OUT_DIR, "tumour_tpm.rds"))
tum_nor_tpm <- computetpm(tum_nor_counts, txi_tum_nor$length)
saveRDS(tum_nor_tpm, file.path(OUT_DIR, "tumour_normal_tpm.rds"))
