# Compare ssGSEA scores between matched primary and normal samples in TRACERx Renal


rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(nlme)
library(harmonicmeanp)
library(ggthemes)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
SSGSEA_NORMALS_PATH <- file.path(BASE, "data", "processed", "tx_normal_ssGSEA.tsv")
SSGSEA_TUMOURS_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
HALLMARK_GROUPS_PATH <- file.path(BASE, "data", "meta", "martinez_ruiz_2023_hallmark_gs_groups.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))


# LOAD DATA ---------------------------------------------------------------
gene_groups <- read_delim(HALLMARK_GROUPS_PATH, delim = "\t")
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
ssgsea_tumour <- read_delim(SSGSEA_TUMOURS_PATH)
ssgsea_tumour$patient <- str_remove(ssgsea_tumour$sample, "[-_].*$")
ssgsea_normal <- read_delim(SSGSEA_NORMALS_PATH)
ssgsea_normal$patient <- str_extract(ssgsea_normal$sample, "[A-Z]_([^_]+)_.*$", 1)

# subset ssgsea_tumour to only patients with matched primary normal sample
prim_smps <- annotation[annotation$type_collapsed == "PRIMARY", ]$sample
ssgsea_tumour <- ssgsea_tumour[ssgsea_tumour$sample %in% prim_smps, ]
ssgsea_tumour <- ssgsea_tumour[ssgsea_tumour$patient %in% ssgsea_normal$patient, ]

ssgsea_tumour$sample_type <- "primary_tumour"
ssgsea_normal$sample_type <- "normal"
ssgsea_df <- rbind(ssgsea_normal, ssgsea_tumour)


# LME COMPARING SAMPLE TYPES ----------------------------------------------

hm_sigs <- colnames(ssgsea_df)[str_detect(colnames(ssgsea_df), "HALLMARK")]
ssgsea_difs_df <- data.frame(sig = c(), t_value = c(), p_value = c())
for (sig in hm_sigs) {
    res <- coef(summary(run_lme(sig, "sample_type", "patient", ssgsea_df)))
    ssgsea_difs_df <- rbind(ssgsea_difs_df, data.frame(sig = sig, t_value = res[2, 4], p_value = res[2, 5]))
}

ssgsea_difs_df$p_adj <- p.adjust(ssgsea_difs_df$p_value)

# collapse by gene groups
ssgsea_difs_df$Hallmark <- str_remove(str_to_lower(ssgsea_difs_df$sig), "hallmark_")
ssgsea_difs_df <- merge(ssgsea_difs_df, gene_groups, by = "Hallmark")

write_delim(ssgsea_difs_df, file.path(OUT_DIR, "ST1_ssgsea_difs_paired_primary_normal.tsv"), delim = "\t")
write_delim(ssgsea_difs_df, file.path(OUT_DIR, "ST2_ssgsea_difs_paired_primary_normal.csv"), delim = ",")

p <- ssgsea_difs_df %>%
    filter(Functional_group != "immune") %>% # we focus in next sections on TME
    group_by(Functional_group) %>%
    summarise(p_value = hmp.stat(p_value), t_value = mean(t_value)) %>%
    mutate(signif = (p_value <= .05)) %>%
    ggplot(aes(x = t_value, y = Functional_group, fill = signif)) +
    geom_segment(aes(x = 0, y = Functional_group, xend = t_value, yend = Functional_group), color = "grey50") +
    geom_vline(aes(xintercept = 0), color = "grey35", size = 1) +
    geom_point(size = 3, pch = 21) +
    scale_fill_manual(values = c("TRUE" = "coral", "FALSE" = "grey80")) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    lemon::coord_capped_cart(bottom = "left", left = "both") +
    theme(
        axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "none"
    ) +
    theme(panel.border = element_blank(), axis.line = element_line())


save_ggplot(p, file.path(FIG_DIR, "ssgsea_paired_primary_normal_lolli"), h = 45, w = 50)
