# Compare I-TED scores obtained using top500 most variable genes or more

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(energy)
library(doMC)
library(argparser)
library(here)

# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "outputs")
PLOT_DIR = file.path(BASE, "analysis", "figures")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")
ITED500_PATH <- file.path(BASE, "analysis", "outputs", "ITED_matrix_primary.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "get_ited.R"))
argparser <- arg_parser("")
argparser <- add_argument(argparser, "--ngenes", help = "number of top variable genes where to run I-TED", type = "numeric")
args <- parse_args(argparser)
n_genes <- args$ngenes

# LOAD DATA ---------------------------------------------------------------
annotation = read_delim(ANNOTATION_PATH, delim = "\t")
vst = readRDS(VST_PATH)
annotation$ITH <- as.numeric(annotation$ITH)
annotation$wgii <- as.numeric(annotation$wgii)
annotation$Regions <- as.numeric(annotation$Regions)

annotation <- annotation[annotation$sample != "K207-R4", ]
vst <- vst[,colnames(vst) != "K207-R4"]

# check that it is good to go
all(annotation$sample == colnames(vst))

print("all loaded")

# Filter to genes with at least 5 counts in at least 20% of the cohort
# The filter here is more strict to avoid high variability in genes with
# problematic detection only in a subset of samples
# Take number of genes for which I-TED is wished to be estimated
vst <- assay(vst)
keep_genes <- (rowMeans(vst > 5) > .2)
vst <- vst[keep_genes, ] # 16704 genes passed filtering
top_genes <- apply(vst, 1, sd) %>%
    sort(decreasing = TRUE) %>%
    head(n = n_genes) %>%
    {
        names(.)
    }

# CALCULATE TRANSCRIPTIONAL DISTANCE --------------------------------------
d_mat = get_dist_matrix_pl(assay(vst), rownames(vst), f = "dcor", cores = 6)
saveRDS(d_mat, file.path(OUT_DIR, "ITED_matrix_allgenes.rds"))

d_mat = readRDS(file.path(OUT_DIR, "ITED_matrix_allgenes.rds"))

# Now run for PRIMARY only
annotation <- annotation %>% 
  filter(type_collapsed == "PRIMARY")

# get samples with more than one tumour primary regions
keep_patients <- annotation %>% 
  dplyr::group_by(Patient) %>% dplyr::summarise(n = n()) %>% # count number of different samples of each patient
  filter(n > 1) %>% dplyr::select(Patient) %>% as_vector() # patients with +1 primary region

keep_samples <- annotation[annotation$Patient %in% keep_patients, "sample"] %>%
  as_vector()

annotation <- annotation[annotation$sample %in% keep_samples, ]
vst <- vst[ ,colnames(vst) %in% keep_samples]

print("Running I-TED for only primary tumour samples...")
d_mat = get_dist_matrix_pl(vst, top_genes, f = "dcor", cores = 6)

write_rds(d_mat, file.path(OUT_DIR, glue::glue("ITED_matrix_primary_ngenes_{n_genes}.rds")))