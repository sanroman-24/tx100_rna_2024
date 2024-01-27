# Get TCR and BCR similarity in TRACERx Renal

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(immunarch)
library(tidyverse)
library(here)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
PLOT_DIR <- file.path(BASE, "analysis", "figures")
BCR_CLONES_PATH <- file.path(BASE, "data", "raw", "tumour_BCR_clones.RDS")
TCR_CLONES_PATH <- file.path(BASE, "data", "raw", "tumour_TCR_clones.RDS")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))

# LOAD DATA ---------------------------------------------------------------
tcr_clonotypes <- readRDS(TCR_CLONES_PATH)
bcr_clonotypes <- readRDS(BCR_CLONES_PATH)
annotation <- read_delim(ANNOTATION_PATH)
annotation$sample <- str_replace(annotation$sample, "-", "_")


# CALCULATE BCR/TCR OVERLAPS ----------------------------------------------

# Takes a while to run; uncomment below if you want to rerun

# bcr_rep_similarity <- immunarch::repOverlap(
#     bcr_clonotypes$data,
#     .method = "morisita"
# )
# colnames(bcr_rep_similarity) <- str_remove(colnames(bcr_rep_similarity), "^[RG]_") %>%
#     str_remove(".clonotypes.IGK_IGL") %>%
#     str_remove("_r\\d") %>%
#     str_replace("K328R19", "K328-R19") %>%
#     str_replace("-", "_")

# colnames(bcr_rep_similarity) <- ifelse(colnames(bcr_rep_similarity) == "K328_T1",
#     "K328_THR1",
#     ifelse(colnames(bcr_rep_similarity) == "K245_T1",
#         "K245_THR1", colnames(bcr_rep_similarity)
#     )
# )
# rownames(bcr_rep_similarity) <- colnames(bcr_rep_similarity)
# saveRDS(bcr_rep_similarity, file.path(OUT_DIR, "bcr_repertoire_similarity.RDS"))
bcr_sim = readRDS(file.path(OUT_DIR, "bcr_repertoire_similarity.RDS"))

# tcr_rep_similarity <- immunarch::repOverlap(
#     tcr_clonotypes$data,
#     .method = "cosine"
# )

# colnames(tcr_rep_similarity) <- str_remove(colnames(tcr_rep_similarity), "^[RG]_") %>%
#     str_remove(".clonotypes.TCRA_TCRB") %>%
#     str_remove("_r\\d") %>%
#     str_replace("K328R19", "K328-R19") %>%
#     str_replace("-", "_")

# colnames(tcr_rep_similarity) <- ifelse(colnames(tcr_rep_similarity) == "K328_T1",
#     "K328_THR1",
#     ifelse(colnames(tcr_rep_similarity) == "K245_T1",
#         "K245_THR1", colnames(tcr_rep_similarity)
#     )
# )
# rownames(tcr_rep_similarity) <- colnames(tcr_rep_similarity)
# saveRDS(tcr_rep_similarity, file.path(OUT_DIR, "tcr_repertoire_similarity.RDS"))
tcr_sim = readRDS(file.path(OUT_DIR, "tcr_repertoire_similarity.RDS"))


# GET SUMMARY SIMILARITY PER PATIENT -------------------------------------------
# consider only primary regions
annotation <- annotation %>% 
  filter(type_collapsed == "PRIMARY")

# get samples with more than one tumour primary regions
keep_patients <- annotation %>% 
  dplyr::group_by(Patient) %>% 
  summarise(n = n()) %>%
  filter(n > 1) %>% dplyr::select(Patient) %>% as_vector() 


summary_sim = data.frame(patient = keep_patients)

summary_sim$tcr_min_sim = summarise_sim(tcr_sim, summary_sim$patient, "min")
summary_sim$tcr_max_sim = summarise_sim(tcr_sim, summary_sim$patient, "max")
summary_sim$tcr_median_sim = summarise_sim(tcr_sim, summary_sim$patient, "median")

summary_sim$bcr_min_sim = summarise_sim(bcr_sim, summary_sim$patient, "min")
summary_sim$bcr_max_sim = summarise_sim(bcr_sim, summary_sim$patient, "max")
summary_sim$bcr_median_sim = summarise_sim(bcr_sim, summary_sim$patient, "median")

tcr_sim_pairs = unlist(summarise_sim(tcr_sim, unname(keep_patients), "c"))
bcr_sim_pairs = unlist(summarise_sim(bcr_sim, unname(keep_patients), "c"))
sim_pairs_df = data.frame(pair_id = names(tcr_sim_pairs), 
                          tcr_sim = unname(tcr_sim_pairs),
                          bcr_sim = unname(bcr_sim_pairs))
sim_pairs_df$patient = str_sub(sim_pairs_df$pair_id, 1, 4)

sim_pairs_df$patient = factor(sim_pairs_df$patient)
sim_pairs_df = sim_pairs_df %>% mutate(patient_ord = fct_reorder(patient, tcr_sim))

saveRDS(summary_sim, file.path(OUT_DIR, "summary_tcr_bcr_sim.rds"))

# PLOT TCR SIMILARITY -----------------------------------------------------

p = ggplot(sim_pairs_df, aes(x = patient_ord, y = tcr_sim)) +
  geom_point(alpha = 0.3, col = "lightblue") 

p = p + 
  geom_point(data = summary_sim, aes(x = patient, y = tcr_max_sim), 
  shape = 15, col = "orchid", alpha = 0.5) + 
  geom_point(data = summary_sim, aes(x = patient, y = tcr_min_sim), 
  shape = 15, col = "orchid", alpha = 0.5) + 
  geom_point(data = summary_sim, aes(x = patient, y = tcr_median_sim), 
  shape = 23, col = "black", fill = "blueviolet") + 
  labs(x = "", y = "TCR Similarity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  lemon::coord_capped_cart(bottom = 'both', left = 'both') 

save_ggplot(p, file.path(PLOT_DIR, "Fig6A_TCR_sim_summary"), w = 180, h = 70)

# PLOT BCR SIMILARITY -----------------------------------------------------
p = ggplot(sim_pairs_df, aes(x = patient_ord, y = bcr_sim)) +
  geom_point(alpha = 0.3, col = "lightblue") 

p = p + 
  geom_point(data = summary_sim, aes(x = patient, y = bcr_max_sim), 
  shape = 15, col = "orchid", alpha = 0.5) + 
  geom_point(data = summary_sim, aes(x = patient, y = bcr_min_sim), 
  shape = 15, col = "orchid", alpha = 0.5) + 
  geom_point(data = summary_sim, aes(x = patient, y = bcr_median_sim), 
  shape = 23, col = "black", fill = "blueviolet") + 
  labs(x = "", y = "BCR Similarity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + 
  lemon::coord_capped_cart(bottom = 'both', left = 'both') 

save_ggplot(p, file.path(PLOT_DIR, "Fig6A_BCR_sim_summary"), w = 180, h = 70)