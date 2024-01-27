# Calculate top 500 genes with highest variation in TRACERx Renal
# (Used to calculate I-TED as in Martinez-Ruiz et al, Nature 2023)

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(here)


# PATHS -------------------------------------------------------------------
BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "outputs")
VST_PATH = file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")

# LOAD DATA ---------------------------------------------------------------
vst = readRDS(VST_PATH)

# FIND TOP 500 WITH HIGHEST VARIATION -------------------------------------

# Filter to genes with at least 5 counts in at least 20% of the cohort
# The filter here is more strict to avoid high variability in genes with
# problematic detection only in a subset of samples
vst = assay(vst)
keep_genes = (rowMeans(vst > 5) > .2) 
vst = vst[keep_genes, ] # 16704 genes passed filtering
top500_genes = apply(vst, 1, sd) %>% sort(decreasing = TRUE) %>% head(n = 500) %>% {names(.)}
saveRDS(top500_genes, file.path(OUT_DIR, "top500_variable_genes.rds"))

