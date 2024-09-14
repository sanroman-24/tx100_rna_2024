# Analysis of the effect of subclonal drivers in TRACERx Renal cohort


rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(stringr)
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
SSGSEA_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")
HALLMARK_GROUPS_PATH <- file.path(META_DIR, "martinez_ruiz_2023_hallmark_gs_groups.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))

pair_comparison <- function(df, pth_str, gene) {
    dir.create(file.path(FIG_DIR, "driver_comparison"))
    df[[gene]] <- ifelse(df[[gene]] == 0, "wt", "mt")
    df_to_ggpaired(
        df = df, pat_var = "Patient",
        cond = gene, var_str = pth_str
    ) %>%
        plot_paired_boxplot(
            cond1 = "wt", cond2 = "mt",
            ylab = pth_str, xlab = "", ylim = c(NA)
        )

    # exp %>%
    #     ggpubr::ggpaired(cond1 = "mt", cond2 = "wt") +
    #     ggpubr::stat_compare_means(paired = TRUE, size = 2.5) +
    #     labs(y = pth_str, x = "")
    ggsave(file.path(FIG_DIR, "driver_comparison", paste0(gene, "paired_comparison_", pth_str, ".pdf")), height = 45, width = 45, units = "mm", dpi = 300)

    exp <- df %>%
        group_by(Patient, !!sym(gene)) %>%
        summarise(avg_exp = mean(!!sym(pth_str)))
    exp <- exp %>% pivot_wider(names_from = !!sym(gene), values_from = avg_exp)
    wtest <- wilcox.test(exp$mt, exp$wt, paired = TRUE)
    return(data.frame(
        gene = gene, pathway = pth_str,
        avg_mt_minus_wt = mean(exp$mt - exp$wt),
        p_val = wtest$p.value
    ))
}


# LOAD DATA ---------------------------------------------------------------
gene_groups <- read_delim(HALLMARK_GROUPS_PATH, delim = "\t")
ssgsea <- read_delim(SSGSEA_PATH, delim = "\t")
pathways <- colnames(ssgsea)[-1]
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")

# merge datasets to have the sample annotation combined with the ssGSEA scores
ssgsea <- left_join(ssgsea, annotation, by = "sample")
ssgsea <- ssgsea[ssgsea$type_collapsed == "PRIMARY", ]

# CALCULATE DIFFERENCES IN SSGSEA BASED ON MUT STATUS ---------------------

drivers <- c("loss_9p", "loss_14q", "SETD2", "PBRM1")

ssgsea_difs_df <- data.frame(gene = c(), pathway = c(), avg_dif = c(), p_val = c())

for (driver in drivers) {
    ssgsea[[driver]] <- ifelse(ssgsea[[driver]] == 2, 1, ssgsea[[driver]])
    sub_ssgsea_df <- find_subclonal_alteration(ssgsea, driver)
    for (pathway in pathways) {
        ssgsea_difs_df <- rbind(ssgsea_difs_df, pair_comparison(sub_ssgsea_df, pathway, driver))
    }
}

ssgsea_difs_df <- ssgsea_difs_df[ssgsea_difs_df$gene %in% c("loss_9p", "loss_14q"), ]

write_delim(ssgsea_difs_df, file.path(OUT_DIR, "ssgsea_difs_9p_14q.tsv"), delim = "\t")

# EXPLORE RESULTS WITH VISUALIZATION -LOG10 FDR * SIGN DIFF ---------------
ssgsea_difs_df$Hallmark <- str_remove(str_to_lower(ssgsea_difs_df$pathway), "hallmark_")
ssgsea_difs_df <- merge(ssgsea_difs_df, gene_groups, by = "Hallmark")

ssgsea_difs_df <- ssgsea_difs_df %>%
    arrange(Functional_group) %>%
    # filter(Functional_group != "immune") %>%
    mutate(id = 1:nrow(.)) %>%
    mutate(sig = p_val < 0.05) %>%
    mutate(sign = -log10(p_val) * sign(avg_mt_minus_wt)) %>%
    mutate(sign = case_when(
        sign < -3 ~ -3,
        sign > 3 ~ 3,
        TRUE ~ sign
    )) %>%
    mutate(pathway = str_replace_all(pathway, "_", " ")) %>%
    mutate(pathway = str_remove_all(pathway, "HALLMARK")) %>%
    mutate(pathway = fct_reorder(pathway, id))

p <- plot_tile(
    df = ssgsea_difs_df, x_str = "pathway", y_str = "gene",
    fill_str = "sign", lgd = "no", sig = "sig"
)

save_ggplot(p, file.path(FIG_DIR, "Fig3B_ssgsea_paired_driver"), w = 120, h = 50)

p <- plot_tile(
    df = ssgsea_difs_df, x_str = "pathway", y_str = "gene",
    fill_str = "sign", lgd = "yes", sig = "sig"
)

save_ggplot(p, file.path(FIG_DIR, "Fig3B_ssgsea_paired_driver_legend"), w = 180, h = 70)

# same plot, but reducing the number of pathways with the data from Carlos Martinez-Ruiz et al, Nature 2023
ssgsea_difs_df <- ssgsea_difs_df %>%
    # filter(Functional_group != "immune") %>% # focus on immune on different section of the paper
    group_by(Functional_group, gene) %>%
    summarise(p_val = hmp.stat(p_val), avg_mt_minus_wt = mean(avg_mt_minus_wt)) %>%
    mutate(sig = p_val < 0.05) %>%
    mutate(sign = -log10(p_val) * sign(avg_mt_minus_wt)) %>%
    mutate(sign = case_when(
        sign < -3 ~ -3,
        sign > 3 ~ 3,
        TRUE ~ sign
    ))

ssgsea_difs_df$Functional_group <- ssgsea_difs_df$Functional_group %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    str_replace("Dna", "DNA")

ssgsea_difs_df$gene <- ssgsea_difs_df$gene %>%
    str_replace("_", " ") %>%
    str_to_title()


p <- plot_tile(ssgsea_difs_df,
    x_str = "Functional_group", y_str = "gene",
    fill_str = "sign", lgd = "yes", sig = "sig"
)

save_ggplot(p, file.path(FIG_DIR, "ssgsea_paired_groups_lgd"), w = 50, h = 50)

p <- plot_tile(ssgsea_difs_df,
    x_str = "Functional_group", y_str = "gene",
    fill_str = "sign", lgd = "no", sig = "sig"
)

save_ggplot(p, file.path(FIG_DIR, "ssgsea_paired_groups"), w = 90, h = 50)
