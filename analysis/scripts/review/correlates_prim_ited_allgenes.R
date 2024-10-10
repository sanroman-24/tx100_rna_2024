# Using I-TED estimated with different total number of genes...
# Associate I-TED in multivariable regression with:
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
OUT_DIR <- file.path(BASE, "analysis", "outputs", "review")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
CNDIST_PATH <- file.path(META_DIR, "cndist_pairs_samps.tsv")

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
    pat_pur_ith <- c(
      pat_pur_ith,
      sd(
        annotation$purity[annotation$Patient == pat]
      )
    )
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
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]
cndist_pairs <- read_delim(CNDIST_PATH, delim = "\t")
pat_subSCNA <- get_pat_subSCNA(cndist_pairs)

wrap_correlates_ngenes <- function(ngenes) {
  TD_MAT_PATH <- file.path(
    BASE, "analysis", "outputs",
    "review", "ITED_matrix_ngenes",
    glue::glue("ITED_matrix_primary_ngenes_{ngenes}.rds")
  )
  d_mat <- readRDS(TD_MAT_PATH)
  pats <- unique(str_extract(rownames(d_mat), "^[A-Z0-9]+"))

  iteds <- sapply(pats, function(patient) {
    summarise_ited_per_patient(d_mat, patient, "median")
  })
  iteds_df <- data.frame(patient = pats, ited = iteds)

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

  drivers <- which(colnames(annotation) == "VHL"):(which(colnames(annotation) == "loss_14q"))
  drivers <- colnames(annotation)[drivers]
  subclonal_alts <- annotation %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise_at(
      drivers,
      is_sub
    )

  subclonal_alts$epi_alt <- sapply(1:nrow(subclonal_alts), function(i) {
    any(
      c(
        subclonal_alts[["SETD2"]][i] == 1,
        subclonal_alts[["PBRM1"]][i] == 1,
        subclonal_alts[["BAP1"]][i] == 1,
        subclonal_alts[["ARID1A"]][i] == 1,
        subclonal_alts[["KDM5C"]][i] == 1
      )
    )
  })

  subclonal_alts <- add_subclonal_cn(subclonal_alts)
  subclonal_alts$patient <- subclonal_alts$Patient

  annotation$ITH_scaled <- log(annotation$ITH + 1)

  purity_ith <- calculate_purity_ith(annotation)


  # MULTIVARIABLE LINEAR REGRESSION -----------------------------------------
  iteds_df <- merge(iteds_df, purity_ith, by = "patient") %>%
    merge(pat_subSCNA, by = "patient") %>%
    merge(subclonal_alts, by = "patient") %>%
    left_join(annotation[!duplicated(annotation$patient), ], by = "patient")
  iteds_df$Subject <- iteds_df$patient
  iteds_df <- left_join(iteds_df, clinical_data, by = "Subject")
  iteds_df$Stage <- ifelse(
    iteds_df$`Overall Stage` %in% c("I", "II"), "I-II", "III-IV"
  )
  iteds_df$size <- as.numeric(iteds_df[["Size of primary tumour (mm)"]])

  form <- "ited ~ pur_ith + subclonal_scna + Stage + Regions + size + ITH_scaled + epi_alt + subclonal_loss_9p + subclonal_loss_14q"
  af <- get_perc_variation(form, df = iteds_df)
  rownames(af) <- c(
    "Purity ITH", "% Genome Subclonal CNA", "Stage", "Regions",
    "Size", "Genetic ITH", "Subclonal epigenetic alteration",
    "Subclonal 9p loss",
    "Subclonal 14q loss", "Residuals"
  )
  p <- plot_perc_variation(af)
  p <- p + ggtitle(glue::glue("I-TED {ngenes}"))
  p <- p + theme(legend.position = "none")

  return(p)
}


p1000 <- wrap_correlates_ngenes(1000)
p2500 <- wrap_correlates_ngenes(2500)
p5000 <- wrap_correlates_ngenes(5000)
p10000 <- wrap_correlates_ngenes(10000)
p12500 <- wrap_correlates_ngenes(12500)
p15000 <- wrap_correlates_ngenes(15000)

lp <- list(p1000, p2500, p5000, p10000, p12500, p15000)

save_plist(lp, file.path(FIG_DIR, "ITED_correlates_ngenes"), w = 140, h = 200, ncol = 2)
