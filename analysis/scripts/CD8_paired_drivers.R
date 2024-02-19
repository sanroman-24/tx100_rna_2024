# Analysis of the effect of subclonal drivers in TRACERx Renal cohort


rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(here)
library(harmonicmeanp)
library(psych)
library(ggthemes)


# PATHS -------------------------------------------------------------------


BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
TME_PATH <- file.path(BASE, "data", "processed", "tum_consensustme.rds")
HALLMARK_GROUPS_PATH <- HALLMARK_GROUPS_PATH <- file.path(META_DIR, "martinez_ruiz_2023_hallmark_gs_groups.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
tme <- as.data.frame(t(readRDS(TME_PATH)))
tme <- rownames_to_column(tme, "sample")
cells <- colnames(tme)[-1]
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")

# merge datasets to have the sample annotation combined with the ssGSEA scores
tme <- left_join(tme, annotation, by = "sample")
tme <- tme[tme$type_collapsed == "PRIMARY", ]

# CALCULATE DIFFERENCES IN CD8 BASED ON MUT STATUS ---------------------

drivers <- c("PBRM1", "SETD2", "loss_9p", "loss_14q")

plist <- lapply(drivers, function(driver) {
    tme[[driver]] <- ifelse(tme[[driver]] == 2, 1, tme[[driver]])
    sub_tme_df <- find_subclonal_alteration(tme, driver)
    sub_tme_df <- dplyr::select(sub_tme_df, c("sample", "Patient", driver, "T_cells_CD8"))
    sub_tme_df[[driver]] <- ifelse(sub_tme_df[[driver]] == 0, "WT", "Mutant")
    p <- df_to_ggpaired(
        df = sub_tme_df, pat_var = "Patient",
        cond = driver, var_str = "T_cells_CD8"
    ) %>%
        plot_paired_boxplot(
            cond1 = "WT", cond2 = "Mutant",
            ylab = "CD8 T cells",
            xlab = "", ylim = c(NA)
        )
    p = p + ggtitle(driver)
    # save_ggplot(p,
    #     file.path(FIG_DIR, paste0("paired_", driver, "_T_CD8")),
    #     w = 80, h = 80
    # )
})

save_plist(plist,
    fp = file.path(FIG_DIR, "SupFig16_CD8_paired_driver"),
    w = 150, h = 150, ncol = 2
)

