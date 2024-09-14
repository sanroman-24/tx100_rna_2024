# Survival analysis in TRACERx Renal based on I-TED estimation

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(gridExtra)
library(grid)

# PATHS -------------------------------------------------------------------

BASE = here::here()
PLOT_DIR = file.path(BASE, "analysis", "figures")
TD_MAT_PATH = file.path(BASE, "analysis", "outputs", "ITED_matrix_primary.rds")
CLINICAL_ANNOTATION_PATH <- file.path(BASE, "data", "meta", "TRACERx_s1_1_clinical.txt")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))

# LOAD DATA ---------------------------------------------------------------
d_mat = readRDS(TD_MAT_PATH)
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data<- clinical_data[-1, ]


# SURVIVAL ANALYSIS BASED ON I-TED ----------------------------------------
pats = unique(str_extract(rownames(d_mat), "^[A-Z0-9]+"))

iteds = sapply(pats, function(patient){summarise_ited_per_patient(d_mat, patient, "median")})
iteds_df = data.frame(patient = pats, ited = iteds)

iteds_df$Subject = iteds_df$patient
iteds_df = merge(iteds_df, clinical_data, by = "Subject") %>% 
  mutate(pfs = ifelse(`PFS (months)` == "-", 0, 1)) %>% 
  mutate(pfs_time = as.numeric(ifelse(`PFS (months)` == "-", `Total follow up (months)`, `PFS (months)`))) %>% 
  mutate(OS = ifelse(Outcome == "Death", 1, 0)) %>% 
  mutate(OS_time =  as.numeric(`Total follow up (months)`)) %>% 
  mutate(Stage = ifelse(`Overall Stage` %in% c("I", "II"), "I-II", "III-IV"))

# PFS KM
iteds_df$ITED = ifelse(iteds_df$ited > median(iteds_df$ited), "High I-TED", "Low I-TED") 
km_fit <- survfit(Surv(pfs_time, pfs) ~ ITED, 
                  data=iteds_df)

km_plot1 <- plot_km(iteds_df, km_fit, plt = c(tx_palette[["darkred"]], "grey50"))

save_baseplot(km_plot1, file.path(PLOT_DIR, "Fig1F_transcriptional_ited_PFS"), h = 90, w = 70)

# OS KM
km_fit <- survfit(Surv(OS_time, OS) ~ ITED, 
                  data=iteds_df)

km_plot2 <- plot_km(iteds_df, km_fit, plt = c(tx_palette[["darkred"]], "grey50"))

save_baseplot(km_plot2, file.path(PLOT_DIR, "SuppFig7_transcriptional_ited_OS"), h = 90, w = 70)

