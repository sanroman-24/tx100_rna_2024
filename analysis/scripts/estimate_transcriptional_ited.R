# Estimate transcriptional intratumour heterogeneity in primary ccRCC

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(DESeq2)
library(energy)
library(pheatmap)
library(ggthemes)


# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "outputs")
PLOT_DIR = file.path(BASE, "analysis", "figures")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")
TOP500_GENES_PATH = file.path(OUT_DIR, "top500_variable_genes.rds")
  
# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))


# LOAD DATA ---------------------------------------------------------------
annotation = read_delim(ANNOTATION_PATH, delim = "\t")
vst = read_rds(VST_PATH)
annotation$ITH <- as.numeric(annotation$ITH)
annotation$wgii <- as.numeric(annotation$wgii)
annotation$Regions <- as.numeric(annotation$Regions)
top500_genes = read_rds(TOP500_GENES_PATH)

annotation <- annotation[annotation$sample != "K207-R4", ]
vst <- vst[,colnames(vst) != "K207-R4"]

# check that it is good to go
all(annotation$sample == colnames(vst))

# CALCULATE TRANSCRIPTIONAL DISTANCE --------------------------------------
d_mat = get_dist_matrix(assay(vst), top500_genes, f = "dcor")
write_rds(d_mat, file.path(OUT_DIR, "ITED_matrix.rds"))

# FILTER TO PATIENTS WITH MORE THAN 1 PRIMARY REGION SAMPLED

annotation <- annotation %>% 
  filter(type_collapsed == "PRIMARY")

# get samples with more than one tumour primary regions
keep_patients <- annotation %>% 
  dplyr::group_by(Patient) %>% dplyr::summarise(n = n()) %>% # count number of different samples of each patient
  filter(n > 1) %>% dplyr::select(Patient) %>% as_vector() # patients with +1 primary region

keep_samples <- annotation[annotation$Patient %in% keep_patients, "sample"] %>%
  as_vector()

annotation_ <- annotation[annotation$sample %in% keep_samples, ]
vst <- vst[ ,colnames(vst) %in% keep_samples]
vst = assay(vst)

d_mat = get_dist_matrix(vst, top500_genes, f = "dcor")
write_rds(d_mat, file.path(OUT_DIR, "ITED_matrix_primary.rds"))

# p = pheatmap(d_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)
# png(file.path(OUTPUT_PLOTS_DIR, "transcriptional_distance_matrix.png"), height = 100, width = 100, units = "mm", res = 300)
# p
# dev.off()


# CALCULATE I-TED ---------------------------------------------------------
summary_ited = data.frame(patient = keep_patients)

summary_ited$min_ited = summarise_ited(d_mat, summary_ited$patient, "min")
summary_ited$max_ited = summarise_ited(d_mat, summary_ited$patient, "max")
summary_ited$median_ited = summarise_ited(d_mat, summary_ited$patient, "median")

ited_pairs = unlist(summarise_ited(d_mat, unname(keep_patients), "c"))
ited_pairs_df = data.frame(pair_id = names(ited_pairs), t_d = unname(ited_pairs))
ited_pairs_df$patient = str_sub(ited_pairs_df$pair_id, 1, 4)


# PLOT I-TED --------------------------------------------------------------
ited_pairs_df$patient = factor(ited_pairs_df$patient)
ited_pairs_df = ited_pairs_df %>% mutate(patient_ord = fct_reorder(patient, t_d))

p = ggplot(ited_pairs_df, aes(x = patient_ord, y = t_d)) + 
  geom_point(alpha = 0.3, col = "lightblue") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

p = p + 
  geom_point(data = summary_ited, aes(x = patient, y = max_ited), shape = 15, col = "orchid", alpha = 0.5) + 
  geom_point(data = summary_ited, aes(x = patient, y = min_ited), shape = 15, col = "orchid", alpha = 0.5) + 
  geom_point(data = summary_ited, aes(x = patient, y = median_ited), shape = 23, col = "black", fill = "blueviolet") + 
  labs(x = "", y = "Transcriptional distance") + 
  theme(axis.text.x = element_text(size = 8))

p = change_axes(p)

save_ggplot(p, file.path(PLOT_DIR, "Fig2E_Transcriptional_ITED_primary"), w = 180, h = 70)


