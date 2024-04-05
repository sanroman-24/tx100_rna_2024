# Survival analysis in TRACERx Renal based on TME-ITH estimation

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
PLOT_DIR <- file.path(BASE, "analysis", "figures")
TD_MAT_PATH <- file.path(BASE, "analysis", "outputs", "TME_dist_primary.rds")
CLINICAL_ANNOTATION_PATH <- file.path(BASE, "data", "meta", "TRACERx_s1_1_clinical.txt")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))

calculate_ited <- function(d_mat, patient) {
  samps <- colnames(d_mat)[str_detect(colnames(d_mat), patient)]
  pat_pairs <- combn(samps, 2)
  t_distances <- c()
  for (pair in 1:ncol(pat_pairs)) {
    s1 <- pat_pairs[1, pair]
    s2 <- pat_pairs[2, pair]
    t_distance <- d_mat[rownames(d_mat) == s1, colnames(d_mat) == s2]
    t_distances <- c(t_distances, t_distance)
  }
  return(median(t_distances))
}

# LOAD DATA ---------------------------------------------------------------
d_mat <- readRDS(TD_MAT_PATH)
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]


# SURVIVAL ANALYSIS BASED ON TME-ITH ----------------------------------------
pats <- unique(str_extract(rownames(d_mat), "^[A-Z0-9]+"))

tme_ith <- sapply(pats, function(patient) {
  summarise_ited_per_patient(d_mat, patient, "median")
})
tme_ith_df <- data.frame(patient = pats, ited = tme_ith)

tme_ith_df$Subject <- tme_ith_df$patient
tme_ith_df <- merge(tme_ith_df, clinical_data, by = "Subject") %>%
  mutate(pfs = ifelse(`PFS (months)` == "-", 0, 1)) %>%
  mutate(pfs_time = as.numeric(ifelse(`PFS (months)` == "-",
    `Total follow up (months)`, `PFS (months)`
  ))) %>%
  mutate(OS = ifelse(Outcome == "Death", 1, 0)) %>%
  mutate(OS_time = as.numeric(`Total follow up (months)`))

# PFS KM
tme_ith_df$high_ited <- tme_ith_df$ited > median(tme_ith_df$ited)
km_fit <- survfit(Surv(pfs_time, pfs) ~ high_ited,
  data = tme_ith_df
)

km_plot <- plot_km(tme_ith_df, km_fit, plt = c("black", "orange"))

save_baseplot(km_plot, file.path(PLOT_DIR, "Fig4C_tme_ith_PFS"), h = 90, w = 70)

# OS KM
km_fit <- survfit(Surv(OS_time, OS) ~ high_ited,
  data = tme_ith_df
)

km_plot <- plot_km(tme_ith_df, km_fit, plt = c("black", "orange"))

# save_baseplot(km_plot, file.path(PLOT_DIR, "SuppFigXX_tme_ith_OS"), h = 90, w = 70)

# Cox regression
coxph(Surv(OS_time, OS) ~ ited, data = tme_ith_df)
