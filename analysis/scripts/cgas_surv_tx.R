# Survival analysis in TRACERx Renal by:
# 1) Expression putative repressors cGAS-STING
# 2) WGII

rm(list = ls(all = TRUE))

# LIBRARIES ------------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggthemes)
library(survival)
library(survminer)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]

vst <- readRDS(VST_PATH)

# WRANGLE DATA ---------------------------------------------------------------
# filter to genes of interest
vst <- assay(vst[rownames(vst) %in% c("ENPP1", "SLC19A1"), ]) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
rownames(vst) <- str_sub(rownames(vst), 1, 12)
# get the maximum wgii, ENPP1 and SLC19A1 per patient
# and merge with annotation for subsequent clinical analysis
vst <- rownames_to_column(vst, "sample") %>%
    mutate(Patient = str_sub(sample, 1, 4))
vst <- vst %>%
    group_by(Patient) %>%
    summarise(ENPP1 = max(ENPP1), SLC19A1 = max(SLC19A1))
annotation <- merge(vst, annotation[!duplicated(annotation$Patient), ])
max_wgii <- annotation %>%
    group_by(Patient) %>%
    summarise(max_wgii = max(wgii, na.rm = TRUE)) %>%
    filter(max_wgii > 0)
annotation <- merge(annotation, max_wgii, by = "Patient")
max_wgii_thresh <- median(annotation$max_wgii)
annotation <- annotation %>%
    mutate(high_wgii = ifelse(max_wgii > max_wgii_thresh,
        "high_wgii", "low_wgii"
    ))

SLC19A1_thresh <- quantile(annotation$SLC19A1, 0.75)
annotation <- annotation %>%
    mutate(high_SLC19A1 = ifelse(SLC19A1 > SLC19A1_thresh, "high_SLC19A1", "low_SLC19A1"))

ENPP1_thresh <- quantile(annotation$ENPP1, 0.75)
annotation <- annotation %>%
    mutate(high_ENPP1 = ifelse(ENPP1 > ENPP1_thresh, "high_ENPP1", "low_ENPP1"))

# find when any of the genes has high expression
annotation <- annotation %>%
    mutate(high_SLCENPP = case_when(
        (high_SLC19A1 == "high_SLC19A1" | high_ENPP1 == "high_ENPP1") ~ "high_SLCENPP",
        (high_SLC19A1 == "low_SLC19A1" & high_ENPP1 == "low_ENPP1") ~ "low_SLCENPP"
    ))

# find when aneuploidy and/or any of the genes high
annotation <- annotation %>%
    mutate(SLCENPP_wgii = case_when(
        (high_SLCENPP == "high_SLCENPP" & high_wgii == "low_wgii") ~ "SLCENPP_high_wgii_low",
        (high_SLCENPP == "high_SLCENPP" & high_wgii == "high_wgii") ~ "SLCENPP_high_wgii_high",
        (high_SLCENPP == "low_SLCENPP" & high_wgii == "low_wgii") ~ "SLCENPP_low_wgii_low",
        (high_SLCENPP == "low_SLCENPP" & high_wgii == "high_wgii") ~ "SLCENPP_low_wgii_high",
    ))

# Add the clinical information from each subject
annotation <- merge(clinical_data, annotation, by.x = "Subject", by.y = "Patient")
annotation <- annotation %>%
    mutate(PFS = ifelse(`PFS (months)` == "-", 0, 1)) %>%
    mutate(`PFS (months)` = as.numeric(ifelse(`PFS (months)` == "-", `Total follow up (months)`, `PFS (months)`))) %>%
    mutate(OS = ifelse(Outcome == "Death", 1, 0)) %>%
    mutate(`Total follow up (months)` = as.numeric(`Total follow up (months)`)) %>%
    dplyr::rename(
        "progression_free_survival" = "PFS", "follow_up_PFS_month" = "PFS (months)",
        "Overall_Survival" = "OS", "total follow up (month)" = "Total follow up (months)"
    )

# KAPLAN MEIER PLOTS ---------------------------------------------------------------
plt <- c("tomato3", "grey25", "tomato1", "grey50")

# PFS
km <- survfit(Surv(follow_up_PFS_month, progression_free_survival) ~ SLCENPP_wgii,
    data = annotation
)

km_plot <- plot_km(annotation, km, plt = plt)
save_baseplot(km_plot, file.path(FIG_DIR, "Fig4E_cgas_wgii_pfs"), h = 90, w = 70)

# OS
km <- survfit(Surv(`total follow up (month)`, Overall_Survival) ~ SLCENPP_wgii,
    data = annotation
)
km_plot <- plot_km(annotation, km, plt = plt)
save_baseplot(km_plot, file.path(FIG_DIR, "Fig4E_cgas_wgii_os"), h = 90, w = 70)
