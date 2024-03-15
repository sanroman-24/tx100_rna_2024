# Compare ssGSEA scores of early/late clones in TRACERx Renal

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(DESeq2)
library(ggpubr)
library(nlme)
library(harmonicmeanp)
library(ggthemes)

# PATHS -------------------------------------------------------------------


BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
SSGSEA_CLONE_PATH <- file.path(BASE, "data", "processed", "ssGSEA_per_clone.rds")
CLONES_PATH <- file.path(META_DIR, "clones_annotation.txt")
HALLMARK_GROUPS_PATH <- file.path(META_DIR, "martinez_ruiz_2023_hallmark_gs_groups.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "get_clonal_dist.R"))
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------

gene_groups <- read_delim(HALLMARK_GROUPS_PATH, delim = "\t")
annotation <- read_delim(ANNOTATION_PATH)
ssgsea <- readRDS(SSGSEA_CLONE_PATH)
clones_annotation <- read_delim(CLONES_PATH, delim = "\t")
clones_annotation <- clones_annotation[!is.na(clones_annotation$clone), ]


# ANNOTATE TERMINAL CLONES ------------------------------------------------

clones_annotation$pat_terminal <- NA
clones_annotation$clone_code <- paste0(
    clones_annotation$patient,
    "-",
    clones_annotation$clone
)

# Find terminal clones
for (pat in unique(clones_annotation$patient)) {
    # find terminal (in leaf) clones in the patient
    term_clones <- get_terminal_clones_pt(clones_annotation, pat)
    # annotate accordingly
    clones_annotation[clones_annotation$patient == pat &
        clones_annotation$clone %in% term_clones, ]$pat_terminal <- "Terminal"
    clones_annotation[clones_annotation$patient == pat &
        !clones_annotation$clone %in% term_clones, ]$pat_terminal <- "No terminal"
}

# Calculate 'clonal age' or number of edges connecting root to leaf
clones_annotation$clonal_age <-
    sapply(seq_len(nrow(clones_annotation)), function(i) {
        pt <- clones_annotation$patient[i]
        clone_id <- clones_annotation$clone_code[i]
        calculate_clonal_age(clones_annotation, pt, clone_id, 0)
    })

# COMPARE ALL THE SSGSEA --------------------------------------------------

ssgsea <- t(ssgsea) %>%
    as.data.frame() %>%
    rownames_to_column("clone_code")
ssgsea$clone_code[83] <- "K207-C4"
ssgsea <- left_join(ssgsea, clones_annotation, by = "clone_code") %>%
    {
        .[!duplicated(.$clone_code), ]
    }

ssgsea <- ssgsea[!is.na(ssgsea$pat_terminal), ]

# make all the paired boxplots
lp = lapply(
    colnames(ssgsea)[str_detect(colnames(ssgsea), "motzer")],
    function(sign) {
        df_to_ggpaired(
            df = ssgsea, pat_var = "patient",
            cond = "pat_terminal", var_str = sign
        ) %>%
            plot_paired_boxplot(
                cond1 = "No terminal", cond2 = "Terminal",
                ylab = sign, xlab = "", ylim = c(NA)
            )
    }
)

# lp = lapply(
#   colnames(ssgsea)[str_detect(colnames(ssgsea), "HALLMARK")],
#   function(sign) {
#     df_to_ggpaired(
#       df = ssgsea, pat_var = "patient",
#       cond = "pat_terminal", var_str = sign
#     ) %>%
#       plot_paired_boxplot(
#         cond1 = "No terminal", cond2 = "Terminal",
#         ylab = sign, xlab = "", ylim = c(NA)
#       )
#   }
# )

# Highlight Motzer T effector and Myeloid inflammation for Fig4D
p = df_to_ggpaired(
            df = ssgsea, pat_var = "patient",
            cond = "pat_terminal", var_str = "motzer_t_effector"
        ) %>%
            plot_paired_boxplot(
                cond1 = "No terminal", cond2 = "Terminal",
                ylab = sign, xlab = "", ylim = c(NA)
            )

save_ggplot(p, file.path(FIG_DIR, "Fig5D_T_eff"), w = 45, h = 45)

p = df_to_ggpaired(
            df = ssgsea, pat_var = "patient",
            cond = "pat_terminal", var_str = "motzer_myeloid_inflammation"
        ) %>%
            plot_paired_boxplot(
                cond1 = "No terminal", cond2 = "Terminal",
                ylab = sign, xlab = "", ylim = c(NA)
            )

save_ggplot(p, file.path(FIG_DIR, "Fig5D_myeloid"), w = 45, h = 45)


# Perform all the comparisons between MOTZER signature scores using LME

hm_sigs <- colnames(ssgsea)[str_detect(colnames(ssgsea), "motzer")]
ssgsea_difs_df <- data.frame(sig = c(), t_value = c(), p_value = c())
for (sig in hm_sigs) {
  form <- as.formula(paste0(sig, " ~ clonal_age"))
  res <- coef(summary(lme(form, random = ~ 1 | patient, data = ssgsea)))
  ssgsea_difs_df <- rbind(ssgsea_difs_df, data.frame(sig = sig, t_value = res[2, 4], p_value = res[2, 5]))
}

ssgsea_difs_df$p_adj <- p.adjust(ssgsea_difs_df$p_value)

# Perform all the comparisons between HALLAMARK signature scores using LME

hm_sigs <- colnames(ssgsea)[str_detect(colnames(ssgsea), "HALLMARK")]
ssgsea_difs_df <- data.frame(sig = c(), t_value = c(), p_value = c())
for (sig in hm_sigs) {
    form <- as.formula(paste0(sig, " ~ clonal_age"))
    res <- coef(summary(lme(form, random = ~ 1 | patient, data = ssgsea)))
    ssgsea_difs_df <- rbind(ssgsea_difs_df, data.frame(sig = sig, t_value = res[2, 4], p_value = res[2, 5]))
}

ssgsea_difs_df$p_adj <- p.adjust(ssgsea_difs_df$p_value)

# collapse by gene groups
ssgsea_difs_df$Hallmark <- str_remove(str_to_lower(ssgsea_difs_df$sig), "hallmark_")
ssgsea_difs_df <- merge(ssgsea_difs_df, gene_groups, by = "Hallmark")
write_delim(ssgsea_difs_df, file.path(OUT_DIR, "ST4_ssgsea_difs_early_late.tsv"), delim = "\t")

ssgsea_difs_df <- ssgsea_difs_df %>%
    filter(Functional_group != "immune") %>% # we focus in next sections on TME
    group_by(Functional_group) %>%
    summarise(p_value = hmp.stat(p_value), t_value = mean(t_value)) %>%
    mutate(signif = (p_value <= .05))

p <- plot_lolli(
    df = ssgsea_difs_df, x_str = "t_value",
    y_str = "Functional_group", fill_str = "signif",
    fill_manual = c("TRUE" = "coral", "FALSE" = "grey80"),
    xlim = c(-2, 2)
)

save_ggplot(p, file.path(FIG_DIR, "Fig3F_ssgsea_early_late_lolli"), w = 50, h = 45)
