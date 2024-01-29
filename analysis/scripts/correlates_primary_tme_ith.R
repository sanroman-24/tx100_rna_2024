# Associate TME ITH in multivariable regression with:
# - % Subclonal SCNA
# - Genetic ITH
# - Purity ITH
# - Subclonal epigenetic mutations
# - Subclonal genome doubling event
# - Clinical stage
# - No. regions

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(ggthemes)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
CNDIST_PATH <- file.path(META_DIR, "cndist_pairs_samps.tsv")
TD_MAT_PATH <- file.path(BASE, "analysis", "outputs", "TME_dist_primary.rds")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))

get_pat_subSCNA <- function(cndist_pairs) {
  cndist_pairs %>%
    mutate(patient = str_extract(s1, "[A-Z]*_K\\d{3}")) %>%
    dplyr::group_by(patient) %>%
    dplyr::summarise(subclonal_scna = median(fga_subclonal_wgd_corrected)) %>%
    mutate(patient = str_remove(patient, "^[A-Z]_"))
}

calculate_purity_ith <- function(annotation) {
  pats <- unique(annotation$Patient)
  pat_pur_ith <- c()
  for (pat in pats) {
    samps <- annotation[annotation$Patient == pat, ]$sample
    all_pairs <- combn(samps, 2)
    p_dists <- c()
    for (pair in 1:ncol(all_pairs)) {
      s1 <- all_pairs[1, pair]
      s2 <- all_pairs[2, pair]
      p_s1 <- annotation$purity[annotation$sample == s1]
      p_s2 <- annotation$purity[annotation$sample == s2]
      p_dists <- c(p_dists, abs(p_s1 - p_s2))
    }
    pat_pur_ith <- c(pat_pur_ith, median(p_dists, na.rm = TRUE))
  }
  df_out <- data.frame(patient = pats, pur_ith = pat_pur_ith)
  return(df_out)
}

# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH)
annotation$ITH <- as.numeric(annotation$ITH)
annotation$wgii <- as.numeric(annotation$wgii)
annotation$Regions <- as.numeric(annotation$Regions)
annotation$purity <- as.numeric(annotation$purity)
d_mat <- readRDS(TD_MAT_PATH)
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]
cndist_pairs <- read_delim(CNDIST_PATH, delim = "\t")
pat_subSCNA <- get_pat_subSCNA(cndist_pairs)


# SURVIVAL ANALYSIS BASED ON I-TED ----------------------------------------
pats <- unique(str_extract(rownames(d_mat), "^[A-Z0-9]+"))

tme_ith <- sapply(pats, function(patient) {
  summarise_ited_per_patient(d_mat, patient, "median")
})
tme_ith_df <- data.frame(patient = pats, ited = tme_ith)

# FILTER TO PATIENTS WITH MORE THAN 1 PRIMARY REGION SAMPLED --------------

# consider only primary regions
annotation <- annotation %>%
  filter(type_collapsed == "PRIMARY")

# get samples with more than one tumour primary regions
keep_patients <- annotation %>%
  dplyr::group_by(Patient) %>%
  dplyr::summarise(n = n()) %>% # count number of different samples of each patient
  filter(n > 1) %>%
  dplyr::select(Patient) %>%
  as_vector() # patients with +1 primary region

# get samples IDs
keep_samples <- annotation[annotation$Patient %in% keep_patients, "sample"] %>%
  as_vector()

annotation <- annotation[annotation$sample %in% keep_samples, ]
annotation$patient <- annotation$Patient

subclonal_alts <- annotation %>%
  dplyr::group_by(patient) %>%
  dplyr::summarise(
    SETD2_alt = ifelse(all(SETD2 == "1"), "clonal",
      ifelse(sum(SETD2 == "1") > 0, "subclonal", "no")
    ),
    alt_9p = ifelse(all(loss_9p == "1"), "clonal",
      ifelse(sum(loss_9p == "1") > 0, "subclonal", "no")
    ),
    alt_14q = ifelse(all(loss_14q == "1"), "clonal",
      ifelse(sum(loss_14q == "1") > 0, "subclonal", "no")
    ),
    PBRM1_alt = ifelse(all(PBRM1 == "1"), "clonal",
      ifelse(sum(PBRM1 == "1") > 0, "subclonal", "no")
    ),
    subclonal_wgd = sum(Genome.doublings == "1") > 0 &
      sum(Genome.doublings == "0") > 0
  )

subclonal_alts$epi_alt <- ifelse(
  subclonal_alts$SETD2_alt != "no" &
    subclonal_alts$PBRM1_alt != "no", "yes", "no"
)
subclonal_alts$subclonal_9p <- ifelse(
  subclonal_alts$alt_9p == "subclonal", 1, 0
)
subclonal_alts$subclonal_14q <- ifelse(
  subclonal_alts$alt_14q == "subclonal", 1, 0
)
annotation$ITH_scaled <- log(annotation$ITH + 1)

purity_ith <- calculate_purity_ith(annotation)

# MULTIVARIABLE LINEAR REGRESSION -----------------------------------------
tme_ith_df <- merge(tme_ith_df, purity_ith, by = "patient") %>%
  merge(pat_subSCNA, by = "patient") %>%
  merge(subclonal_alts, by = "patient") %>%
  left_join(annotation[!duplicated(annotation$patient), ], by = "patient")
tme_ith_df$Subject <- tme_ith_df$patient
tme_ith_df <- left_join(tme_ith_df, clinical_data, by = "Subject")
tme_ith_df$Stage <- ifelse(
  tme_ith_df$`Overall Stage` %in% c("I", "II"), "I-II", "III-IV"
)

form <- "ited ~ pur_ith + subclonal_scna + Stage + Regions + ITH_scaled + SETD2_alt + PBRM1_alt + subclonal_9p + subclonal_14q + subclonal_wgd"
af <- get_perc_variation(form, df = tme_ith_df)
p <- plot_perc_variation(af)
p <- p + theme(legend.position = "none")

save_ggplot(p, file.path(FIG_DIR, "Fig5B_correlates_tme_ith"), w = 70, h = 45)
