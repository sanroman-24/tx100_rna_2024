# Get TCR and BCR similarity in TRACERx Renal

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(immunarch)
library(tidyverse)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
OUT_DIR <- file.path(BASE, "analysis", "outputs")
PLOT_DIR <- file.path(BASE, "analysis", "figures")
SUM_TCR_BCR_PATH <- file.path(OUT_DIR, "summary_tcr_bcr_sim.rds")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_ited.R"))

# LOAD DATA ---------------------------------------------------------------
summary_sim <- readRDS(SUM_TCR_BCR_PATH)
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]
clinical_data$patient <- clinical_data$Subject
summary_sim <- merge(summary_sim, clinical_data, by = "patient") %>%
    mutate(pfs = ifelse(`PFS (months)` == "-", 0, 1)) %>%
    mutate(pfs_time = as.numeric(ifelse(`PFS (months)` == "-",
        `Total follow up (months)`, `PFS (months)`
    ))) %>%
    mutate(os = ifelse(Outcome == "Death", 1, 0)) %>%
    mutate(os_time = as.numeric(`Total follow up (months)`))

# KAPLAN MEIER PLOTS ---------------------------------------------------------------

plt <- c("black", "orange")

# TCR
summary_sim$high_tcr_sim <-
    summary_sim$tcr_median_sim > median(summary_sim$tcr_median_sim, na.rm = T)

km <- survfit(Surv(pfs_time, pfs) ~ high_tcr_sim,
    data = summary_sim
)

km_plot <- plot_km(summary_sim, km, plt = plt)
save_baseplot(km_plot, file.path(PLOT_DIR, "SuppFig21_TCR_sim_pfs"), h = 90, w = 70)



# BCR
summary_sim$high_bcr_sim <-
    summary_sim$bcr_median_sim > median(summary_sim$bcr_median_sim, na.rm = T)

km <- survfit(Surv(pfs_time, pfs) ~ high_bcr_sim,
    data = summary_sim
)

km_plot <- plot_km(summary_sim, km, plt = plt)
save_baseplot(km_plot, file.path(PLOT_DIR, "SuppFig21_BCR_sim_pfs"), h = 90, w = 70)

# BCR & TCR
summary_sim$tcr_bcr_group <- paste0(
    "tcr_high_", summary_sim$high_tcr_sim,
    "bcr_high_", summary_sim$high_bcr_sim
)
summary_sim <- summary_sim %>%
    filter(!is.na(high_tcr_sim) & !is.na(high_bcr_sim))

km <- survfit(Surv(pfs_time, pfs) ~ tcr_bcr_group,
    data = summary_sim
)
plt <- "lancet"
km_plot <- plot_km(summary_sim, km, plt = plt)
save_baseplot(km_plot, file.path(PLOT_DIR, "Fig5B_BCR_TCR_sim_pfs"),
    h = 90, w = 70
)
