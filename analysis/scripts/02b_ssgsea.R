### Calculate ssGSEA scores in TRACERx Renal tumour and normal samples ###

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(GSVA)
library(msigdbr)
library(DESeq2)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
TUM_TPM_PATH <- file.path(BASE, "data", "processed", "tumour_tpm.rds")
TUM_NOR_TPM_PATH <- file.path(BASE, "data", "processed", "tumour_normal_tpm.rds")
GS_PATH <- file.path(BASE, "data", "meta", "gene_signatures.rds")
OUT_DIR <- file.path(BASE, "data", "processed")

# LOAD DATA ---------------------------------------------------------------

tum_tpm <- readRDS(TUM_TPM_PATH)
tum_nor_tpm <- readRDS(TUM_NOR_TPM_PATH)
all_sign <- readRDS(GS_PATH)

# RUN SSGSEA --------------------------------------------------------------
# Tumour samples
ssGSEA <- gsva(expr = tum_tpm, gset.idx.list = all_sign, method = "ssgsea") %>%
  # transverse to get different signatures in the columns
  t() %>%
  # transform to dataframe so that results can be easily merged to sample annotation
  as.data.frame() %>%
  rownames_to_column(var = "sample")

write_delim(ssGSEA, file.path(OUT_DIR, "tx_ssGSEA.tsv"), delim = "\t")

# Normal samples

# subset only to normals to run ssGSEA here
nor_tpm <- tum_nor_tpm[, str_detect(colnames(tum_nor_tpm), "_N")]
ssGSEA_nor <- gsva(expr = nor_tpm, gset.idx.list = all_sign, method = "ssgsea") %>%
  # transverse to get different signatures in the columns
  t() %>%
  # transform to dataframe so that results can be easily merged to sample annotation
  as.data.frame() %>%
  rownames_to_column(var = "sample")

write_delim(ssGSEA_nor, file.path(OUT_DIR, "tx_normal_ssGSEA.tsv"), delim = "\t")
