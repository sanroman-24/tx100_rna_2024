### Run consensusTME in TRACERx Renal tumour samples ###

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(ConsensusTME)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
TUM_TPM_PATH <- file.path(BASE, "data", "processed", "tumour_tpm.rds")
OUT_DIR <- file.path(BASE, "data", "processed")


# LOAD DATA ---------------------------------------------------------------
tum_tpm <- readRDS(TUM_TPM_PATH)

# RUN CONSENSUSTME --------------------------------------------------------
tum_tme <- ConsensusTME::consensusTMEAnalysis(
    bulkExp = tum_tpm,
    cancer = "KIRC",
    statMethod = "ssgsea"
)

saveRDS(tum_tme, file.path(OUT_DIR, "tum_consensustme.rds"))
