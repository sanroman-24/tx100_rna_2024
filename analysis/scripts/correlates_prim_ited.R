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
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
CNDIST_PATH <- file.path(META_DIR, "cndist_pairs_samps.tsv")
TD_MAT_PATH <- file.path(BASE, "analysis", "outputs", "ITED_matrix_primary.rds")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))
source(file.path(BASE, "src", "scna_ith.R"))

is_sub <- function(v) {
  return(ifelse(sum(v > 0) > 0 & sum(v == 0) > 0, 1, 0))
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

add_subclonal_cn <- function(subclonal_alts) {
  cn_drivers <- colnames(subclonal_alts)
  msk <- grepl("gain", cn_drivers) | grepl("loss", cn_drivers)
  cn_drivers <- cn_drivers[msk]
  for (cn_driver in cn_drivers) {
    v <- paste0("subclonal_", cn_driver)
    subclonal_alts[[v]] <- subclonal_alts[[cn_driver]]
  }
  return(subclonal_alts)
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

sapply(drivers, function(driver) mean(subclonal_alts[[driver]]))
sapply(drivers, function(driver) sum(subclonal_alts[[driver]]))

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


# selecting variables for analysis

# - Variation in purity that can confound I-TED
# - % genome with subclonal SCNAs
# - Overall genetic variation (genetic ITH)
# - Clinical variables: stage (more aggressive tumors)
# - Number of regions: rule out confounder
# - Size: requested by reviewer #1
# - Subclonal epigenetic driver mutation: KDM5C, ARID1A, SETD2, PBRM1, BAP1

form <- "ited ~ pur_ith + subclonal_scna + Stage + Regions + size + ITH_scaled + epi_alt + subclonal_loss_9p + subclonal_loss_14q"
af <- get_perc_variation(form, df = iteds_df)
rownames(af) <- c(
    "Purity ITH", "% Genome with subclonal CNA", "Stage", "Regions",
    "Size", "Genetic ITH", "Subclonal epigenetic alteration",
    "Subclonal 9p loss",
    "Subclonal 14q loss", "Residuals"
  )
p <- plot_perc_variation(af)
p <- p + theme(legend.position = "none")

save_ggplot(p, file.path(FIG_DIR, "Fig1C_correlates_ited"), w = 70, h = 45)



# Fit random I-TED models with other subclonal CN to see relevance 9p loss better
# 14 driver SCNAs previously identified in TRACERx Renal
# See Figure S2 here: https://www.cell.com/cell/fulltext/S0092-8674(18)30389-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418303891%3Fshowall%3Dtrue
subclonal_cn <- colnames(iteds_df)
msk <- grepl("subclonal_loss", subclonal_cn) | grepl("subclonal_gain", subclonal_cn)
subclonal_cn <- subclonal_cn[msk]


df <- data.table::rbindlist(lapply(subclonal_cn, function(event) {
  form <- glue::glue("ited ~ pur_ith + subclonal_scna + Stage + Regions + size + ITH_scaled + epi_alt + {event}")
  af <- get_perc_variation(form, df = iteds_df)
  af <- af[rownames(af) == event, ]
  pval <- af[["Pr(>F)"]]
  perc_var <- af[["perc_var"]]
  return(data.frame(event = event, pval = pval, perc_var = perc_var))
}))

df$asterisks <- sapply(df$pval, function(p) pval_to_asterisks(p))

df$sig <- df$pval < .05

p <- df %>%
  arrange(perc_var) %>%
  mutate(event = str_to_title(str_replace_all(event, "_", " "))) %>%
  mutate(event = factor(event, levels = event)) %>%
  ggplot(aes(x = perc_var, y = event)) +
  geom_col(aes(fill = sig), alpha = 0.7, col = "black") +
  geom_text(aes(label = asterisks, x = perc_var + 0.5),
    vjust = 0.6, hjust = 0, size = 4
  ) +
  lims(x = c(0, 15)) +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = tx_palette[["lightblue"]])) +
  labs(x = "% of variance in I-TED explained", y = "") + 
  theme(legend.position = "none")

save_ggplot(p, file.path(FIG_DIR, "ited_9p_vs_others"), w = 70, h = 50)
