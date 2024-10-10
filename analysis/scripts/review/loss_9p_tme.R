# Analysis of the underpinnings of 9p loss TME (#8 reviewer #1)

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(here)
library(ggthemes)
library(patchwork)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs", "review")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
SSGSEA_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")
TME_PATH <- file.path(BASE, "data", "processed", "tum_consensustme.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
ssgsea <- read_delim(SSGSEA_PATH)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
annotation <- annotation[annotation$type_collapsed == "PRIMARY",]

# 1) COMPARE PURITY IN 9P LOSS VS WILD-TYPE -------------------------------
# subset to only patients with subclonal 9p loss
annotation <- find_subclonal_alteration(annotation, "loss_9p")


# LME -> lower purity in 9p loss; p-value = 0.05
summary(run_lme(
    "purity", "loss_9p", "Patient", annotation[!is.na(annotation$purity),]
))

p1 <- annotation %>% 
    drop_na(purity) %>%
    mutate(loss_9p = ifelse(loss_9p == "1", "9p21 loss", "9p WT")) %>%
    df_to_ggpaired(
        pat_var = "Patient",
        cond = "loss_9p", var_str = "purity"
    ) %>%
        plot_paired_boxplot(
            cond1 = "9p WT", cond2 = "9p21 loss",
            ylab = "Purity",
            xlab = "", ylim = c(NA)
        )

# save_ggplot(p, file.path(FIG_DIR, "purity_paired_9p"), w = 50, h = 50)

# 2) COMPARE MACROPHAGE ABUNDANCE
tme <- as.data.frame(t(readRDS(TME_PATH)))
tme <- rownames_to_column(tme, "sample")
cells <- colnames(tme)[-1]
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")

# merge datasets to have the sample annotation combined with the ssGSEA scores
tme <- left_join(tme, annotation, by = "sample")
tme <- tme[tme$type_collapsed == "PRIMARY", ]
sub_tme_df <- find_subclonal_alteration(tme, "loss_9p")
sub_tme_df <- dplyr::select(sub_tme_df, c("sample", "Patient", "loss_9p", "wgii", "Macrophages"))
sub_tme_df$loss_9p <- ifelse(sub_tme_df$loss_9p == 0, "WT", "9p21 loss")
summary(run_lme(
    "Macrophages", "loss_9p", "Patient", sub_tme_df
))
wgii_thr <- median(tme$wgii, na.rm = T) # 0.23
summary(nlme::lme(Macrophages ~ loss_9p + (wgii > 0.23), random = ~ 1|Patient, data = sub_tme_df[!is.na(sub_tme_df$wgii),]))

p2 <- df_to_ggpaired(
        df = sub_tme_df, pat_var = "Patient",
        cond = "loss_9p", var_str = "Macrophages"
    ) %>%
        plot_paired_boxplot(
            cond1 = "WT", cond2 = "9p21 loss",
            ylab = "Macrophages",
            xlab = "", ylim = c(NA)
        )

# save_ggplot(p,
#     fp = file.path(FIG_DIR, "Macrophages_paired_9p_loss"),
#     w = 50, h = 50
# )

sub_tme_df <- find_subclonal_alteration(tme, "loss_9p")
sub_tme_df <- dplyr::select(sub_tme_df, c("sample", "Patient", "loss_9p", "wgii", "T_cells_CD8"))
sub_tme_df$loss_9p <- ifelse(sub_tme_df$loss_9p == 0, "WT", "9p21 loss")
summary(run_lme(
    "T_cells_CD8", "loss_9p", "Patient", sub_tme_df
))

summary(nlme::lme(T_cells_CD8 ~ loss_9p + (wgii > 0.23), random = ~ 1|Patient, data = sub_tme_df[!is.na(sub_tme_df$wgii),]))


p3 <- df_to_ggpaired(
        df = sub_tme_df, pat_var = "Patient",
        cond = "loss_9p", var_str = "T_cells_CD8"
    ) %>%
        plot_paired_boxplot(
            cond1 = "WT", cond2 = "9p21 loss",
            ylab = "CD8 T Cells",
            xlab = "", ylim = c(NA)
        )

p <- p1 + p2 + p3 + plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_DIR, "loss_9p_tme.png"), width = 200, height = 70, units = "mm", dpi = 300)

# 3) COMPARE SSGSEA SCORES IN MYELOID INFLAMMATION VS T CELL EFFECTOR
ssgsea <- ssgsea[ssgsea$sample %in% annotation$sample,]
ssgsea <- merge(ssgsea, annotation, "sample")

summary(run_lme(
    "motzer_t_effector", "loss_9p", "Patient", ssgsea
))

p <- ssgsea %>% 
    mutate(loss_9p = ifelse(loss_9p == "1", "Loss 9p", "9p WT")) %>%
    df_to_ggpaired(
        pat_var = "Patient",
        cond = "loss_9p", var_str = "motzer_t_effector"
    ) %>%
        plot_paired_boxplot(
            cond1 = "9p WT", cond2 = "Loss 9p",
            ylab = "T cell effector",
            xlab = "", ylim = c(NA)
        )

save_ggplot(p, file.path(FIG_DIR, "motzer_teff_paired_9p"), w = 50, h = 50)

summary(run_lme(
    "motzer_myeloid_inflammation", "loss_9p", "Patient", ssgsea
))

p <- ssgsea %>% 
    mutate(loss_9p = ifelse(loss_9p == "1", "Loss 9p", "9p WT")) %>%
    df_to_ggpaired(
        pat_var = "Patient",
        cond = "loss_9p", var_str = "motzer_myeloid_inflammation"
    ) %>%
        plot_paired_boxplot(
            cond1 = "9p WT", cond2 = "Loss 9p",
            ylab = "ssGSEA Myeloid inflammation",
            xlab = "", ylim = c(NA)
        )

save_ggplot(p, file.path(FIG_DIR, "motzer_myelinf_paired_9p"), w = 50, h = 50)


dif <- read_delim(file.path(BASE, "analysis", "outputs", "ssgsea_difs_9p_14q.tsv"), delim = "\t")

dif %>% filter(gene == "loss_9p") %>% 
    filter(str_detect(pathway, "HALLMARK_")) %>%
    arrange(p_val)


dif$high <- ifelse(
    str_detect(dif$pathway, "HALLMARK_INTERFERON"),
    "yes", "no"
)
dif$path <- ifelse(dif$high == "yes", dif$pathway, "")

p <- dif %>% filter(gene == "loss_9p") %>% 
    filter(str_detect(pathway, "HALLMARK_")) %>%
    ggplot(aes(x = avg_mt_minus_wt, y = -log10(p_val))) + 
    geom_point(aes(fill = high), pch = 21) + 
    ggrepel::geom_text_repel(aes(label = path), size = 1) + 
    geom_hline(aes(yintercept = -log10(0.05)), lty = "dashed") + 
    geom_vline(aes(xintercept = 0), lty = "dashed") + 
    scale_fill_manual(values = c("yes" = "red", "no" = "grey60")) + 
    theme(legend.position = "none") + 
    labs(x = "Difference in ssGSEA 9p loss", y = "-log10 FDR") + 
    xlim(c(-0.04, 0.04))

save_ggplot(p, file.path(FIG_DIR, "IFN_paired_9p"), w = 50, h = 50)
