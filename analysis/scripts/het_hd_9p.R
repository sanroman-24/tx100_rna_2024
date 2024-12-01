# Compare transcriptional programmes in 9p HD vs het loss

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(here)
library(ggthemes)
library(RColorBrewer)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
HALLMARK_GROUPS_PATH <- file.path(META_DIR, "martinez_ruiz_2023_hallmark_gs_groups.txt")
CLINICAL_ANNOTATION_PATH <- file.path(BASE, "data", "meta", "TRACERx_s1_1_clinical.txt")
SSGSEA_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")
COUNTS_PATH <- file.path(BASE, "data", "processed", "tumour_tpm.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))

# MAIN ---------------------------------------------------------------

cdkn2_genes <- c("CDKN2A", "CDKN2B")
ifn_genes <- c(
    "IFNB1", "IFNW1", "IFNA21", "IFNA4",
    "IFNA16", "IFNA7", "IFNA14", "IFNA5",
    "IFNA6", "IFNA13", "IFNA2", "IFNA8",
    "IFNA1", "IFNE"
)

tpm <- readRDS(COUNTS_PATH)
annotation <- read_delim(ANNOTATION_PATH)

## Evaluate global expression CDKN2A/B and IFN I 9p21 genes
exp <- apply(tpm, 1, median)
df <- data.frame(gene = names(exp), mean_tpm = log10(exp + 1)) %>%
    mutate(is_cdkn2 = gene %in% cdkn2_genes, is_ifn = gene %in% ifn_genes) %>%
    mutate(highlight = ifelse(is_cdkn2, "CDKN2A/B", ifelse(is_ifn, "9p IFN cluster", "Other"))) %>%
    mutate(gene = fct_reorder(gene, mean_tpm))

p <- ggplot(df, aes(x = gene, y = mean_tpm, col = highlight, alpha = highlight)) +
    geom_point() +
    scale_color_manual(values = c("Other" = "grey80", "CDKN2A/B" = "red", "9p IFN cluster" = "blue")) +
    scale_alpha_manual(values = c("Other" = 0.01, "CDKN2A/B" = 1, "9p IFN cluster" = 1)) +
    labs(x = "Gene Rank", y = "Average TRACERx Renal log10 (TPM + 1)", color = "Gene") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    guides(alpha = "none")

ggsave(file.path(FIG_DIR, "expression_level_9p_genes.pdf"),
    w = 80, h = 80, units = "mm"
)

# COMPARE DIFFERENCES IN EXPRESSION BY LOH AND HD
df <- tpm[rownames(tpm) %in% c(cdkn2_genes, ifn_genes), ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    merge(annotation)

hd_9p <- read_delim(file.path(META_DIR, "9p_status.csv"))
hd_9p$sample <- hd_9p$ID
hd_9p$sample <- clean_ids(hd_9p$sample)
hd_9p$Patient <- str_extract(hd_9p$sample, "K\\d+")
df$sample <- clean_ids(df$sample)
df <- merge(hd_9p, df, by = "sample")
df <- df[!duplicated(df$sample),]

# CDKN2 GENES

p <- df %>%
    dplyr::select(cdkn2_genes, status_9p, sample) %>%
    pivot_longer(cdkn2_genes, names_to = "cdkn_gene", values_to = "tpm") %>%
    ggplot(aes(x = status_9p, y = log10(tpm + 1))) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_point(alpha = 0.3, size = .2, position = position_dodge2(width = .2)) + 
    ggpubr::stat_compare_means(comparisons = list(
        c("HD", "LOH"),
        c("HD", "WT"),
        c("LOH", "WT")
    ), size = 1.9) +
    facet_wrap(~cdkn_gene, scales = "free_y") +
    theme(
        strip.text = element_text(size = 6, face = "bold")
    ) +
    labs(x = "9p21 status", y = "log10(TPM+1)")

save_ggplot(p, file.path(FIG_DIR, "CDKN2_loss_9p"), w = 60, h = 40)

# IFN GENES
p <- df %>%
    dplyr::select(ifn_genes, status_9p, sample) %>%
    pivot_longer(ifn_genes, names_to = "ifn_gene", values_to = "tpm") %>%
    ggplot(aes(x = status_9p, y = log10(tpm + 1))) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_point(alpha = 0.3, size = .2, position = position_dodge2(width = .2)) + 
    ggpubr::stat_compare_means(comparisons = list(
        c("HD", "LOH"),
        c("HD", "WT"),
        c("LOH", "WT")
    ), size = 2.5) +
    facet_wrap(~ifn_gene, scales = "free_y") +
    theme(
        strip.text = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)
    ) +
    labs(x = "9p21 status", y = "log10(TPM+1)")

save_ggplot(p, file.path(FIG_DIR, "IFN_loss_9p"), w = 150, h = 150)

## COMPARE "PHENOTYPIC EFFECT" ---------------------------------------

# COMPARE TRANSCRIPTIONAL IMPACT IN PATHWAYS
gene_groups <- read_delim(HALLMARK_GROUPS_PATH)
ssgsea <- read_delim(SSGSEA_PATH)
gene_groups$pathway <- paste0("HALLMARK_", str_to_upper(gene_groups$Hallmark))
pf_paths <- gene_groups$pathway[gene_groups$Functional_group == "proliferation"]
ifn_paths <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
scores <- data.frame(
    sample = ssgsea[, 1],
    proliferation = rowMeans(ssgsea[, pf_paths]),
    ifn = rowMeans(ssgsea[, ifn_paths]), 
    ifn_alpha = ssgsea[[ifn_paths[1]]],
    ifn_gamma = ssgsea[[ifn_paths[2]]]
)

scores$sample <- clean_ids(scores$sample)
df <- merge(scores, df)

# PROLIFERATION

p <- ggplot(df, aes(x = status_9p, y = proliferation)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(comparisons = list(
        c("HD", "LOH"),
        c("HD", "WT"),
        c("LOH", "WT")
    ), size = 2.5) +
    labs(x = "9p21 status", y = "ssGSEA proliferation")

ggsave(
    file.path(FIG_DIR, "ssgsea_prolif_hd_loh_wt_9p.pdf"),
    w = 50, h = 50, units = "mm"
)

# IFN PATHWAY
p <- ggplot(df, aes(x = status_9p, y = ifn)) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_point(alpha = 0.3, size = .2, position = position_dodge2(width = .2)) +
    ggpubr::stat_compare_means(comparisons = list(
        c("HD", "LOH"),
        c("HD", "WT"),
        c("LOH", "WT")
    ), size = 2.5) +
    labs(x = "9p21 status", y = "ssGSEA IFN")

ggsave(
    file.path(FIG_DIR, "ssgsea_ifn_hd_loh_wt_9p.pdf"),
    w = 50, h = 50, units = "mm"
)

# IFN alpha
p <- ggplot(df, aes(x = status_9p, y = ifn_alpha)) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_point(alpha = 0.3, size = .2, position = position_dodge2(width = .2)) +
    ggpubr::stat_compare_means(comparisons = list(
        c("HD", "LOH"),
        c("HD", "WT"),
        c("LOH", "WT")
    ), size = 2.5) +
    labs(x = "9p21 status", y = "ssGSEA IFN alpha")

ggsave(
    file.path(FIG_DIR, "ssgsea_ifn_alpha_hd_loh_wt_9p.pdf"),
    w = 50, h = 50, units = "mm"
)

# IFN gamma
p <- ggplot(df, aes(x = status_9p, y = ifn_gamma)) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    geom_point(gamma = 0.3, size = .2, position = position_dodge2(width = .2)) +
    ggpubr::stat_compare_means(comparisons = list(
        c("HD", "LOH"),
        c("HD", "WT"),
        c("LOH", "WT")
    ), size = 2.5) +
    labs(x = "9p21 status", y = "ssGSEA IFN gamma")

ggsave(
    file.path(FIG_DIR, "ssgsea_ifn_gamma_hd_loh_wt_9p.pdf"),
    w = 50, h = 50, units = "mm"
)