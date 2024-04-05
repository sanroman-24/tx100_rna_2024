# Analysis of the abundance of cell types in TME per TRACERx evo subtype

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(here)
library(tidyverse)
library(fastDummies)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
PLOT_DIR <- file.path(BASE, "analysis", "figures")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
TME_PATH <- file.path(BASE, "data", "processed", "tum_consensustme.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))


# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
tme <- read_rds(TME_PATH)
cells <- rownames(tme)

annotation <- annotation[annotation$sample != "K207-R4", ]
tme <- tme[, colnames(tme) != "K207-R4"]
tme <- t(tme) %>%
    as.data.frame() %>%
    rownames_to_column("sample")
tme <- merge(tme, annotation, by = "sample")

tme <- dummy_cols(tme, select_columns = c("EvoType"))
# filter to only primary tumour regions
tme <- tme[str_detect(tme$type, "(PRIMARY)|Primary"), ]

# COMPARE IMMUNE INFILTRATION ---------------------------------------------

# COMPUTE WITH LINEAR MIXED EFFECTS MODEL COMPARISON 1 VS ALL
out <- data.frame(row.names = cells)
evo_dummies <- which(str_detect(colnames(tme), "EvoType_"))
# remove EvoType_undefined, which is not of interest
evo_dummies <- evo_dummies[-8]

for (j in 1:length(evo_dummies)) {
    dummy_evo <- colnames(tme)[evo_dummies[j]]
    for (i in seq_along(cells)) {
        cell <- cells[i]
        lme_fit <- run_lme(cell, dummy_evo, "Patient", tme)
        t_value <- coef(summary(lme_fit))[2, 4]
        p_value <- coef(summary(lme_fit))[2, 5]
        out[i, j] <- -log10(p_value) * sign(t_value)
    }
}

# set column names of the dataframe storing the results from the comparisons
colnames(out) <- colnames(tme)[evo_dummies]

# reshape into a format suitable for geom_tile() to work with for plotting
out <- out %>%
    rownames_to_column(var = "immune_cell") %>%
    reshape2::melt(value.name = "signif")

out <- mutate(out, signif = case_when(
    signif > 3 ~ 3,
    signif < -3 ~ -3,
    TRUE ~ signif
))

# plot considering all the values
p = plot_tile(out, "immune_cell", "variable", "signif", lgd = "yes")
# Change to a more sensible order
p = p + scale_y_discrete(
        expand = c(0, 0),
        limits = c(
            "EvoType_Multiple_Truncal_Drivers", "EvoType_VHL_wt",
            "EvoType_BAP1_driven", "EvoType_PBRM1_sCNA",
            "EvoType_PBRM1_mTOR", "EvoType_PBRM1_SETD2",
            "EvoType_Monodriver_VHL", "EvoType_non_driver_subtype"
        )
    )

# save plot
save_ggplot(p, file.path(PLOT_DIR, "Fig4E_immune_evotypes"), w = 130, h = 60)


# no legend
p = plot_tile(out, "immune_cell", "variable", "signif", lgd = "no")
# Change to a more sensible order
p = p + scale_y_discrete(
        expand = c(0, 0),
        limits = c(
            "EvoType_Multiple_Truncal_Drivers", "EvoType_VHL_wt",
            "EvoType_BAP1_driven", "EvoType_PBRM1_sCNA",
            "EvoType_PBRM1_mTOR", "EvoType_PBRM1_SETD2",
            "EvoType_Monodriver_VHL", "EvoType_non_driver_subtype"
        )
    )

# save plot
save_ggplot(p, file.path(PLOT_DIR, "Fig4E_immune_evotypes_no_lgd"), w = 130, h = 60)