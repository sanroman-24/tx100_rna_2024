# Compare ssGSEA scores between matched primary and metastasis samples in TRACERx Renal
# for CIN signatures, based on reviewer #1 feedback


rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(nlme)
library(harmonicmeanp)
library(ggthemes)
library(patchwork)

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

# subset ssgsea_tumour to only patients with matched primary normal sample
prim_smps <- annotation[annotation$type_collapsed == "PRIMARY", ]$sample
ssgsea_primary <- ssgsea_tumour[ssgsea_tumour$sample %in% prim_smps, ]
met_smps <- annotation[annotation$type_collapsed == "METASTASIS", ]$sample
ssgsea_met <- ssgsea_tumour[ssgsea_tumour$sample %in% met_smps, ]

keep_pats <- dplyr::intersect(ssgsea_primary$patient, ssgsea_met$patient)
ssgsea_primary <- ssgsea_primary[ssgsea_primary$patient %in% keep_pats, ]
ssgsea_met <- ssgsea_met[ssgsea_met$patient %in% keep_pats, ]

ssgsea_primary$sample_type <- "primary_tumour"
ssgsea_met$sample_type <- "met"
ssgsea_df <- rbind(ssgsea_met, ssgsea_primary)

# LME CIN signatures primary vs metastasis (reviewer #1)

summary(run_lme("cin70", "sample_type", "patient", ssgsea_df))
summary(run_lme("HET70", "sample_type", "patient", ssgsea_df))
summary(run_lme("cin_bakhoum", "sample_type", "patient", ssgsea_df))

cin_df <- ssgsea_df %>%
  dplyr::select(cin70, HET70, cin_bakhoum, sample_type, patient, sample) %>%
  mutate(sample_type = ifelse(sample_type == "met", "Metastasis", "Primary tumor"))

dodge <- position_dodge2(width = 0.2)
p1 <- ggplot(cin_df, aes(x = sample_type, y = cin70)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(position = dodge, alpha = 0.2) +
  geom_text(aes(
    x = 1.5, y = max(cin70) + 0.1 * max(cin70),
    label = "LME p-value = 0.006"
  ), size = 2, family = "Arial") +
  labs(x = "", y = "CIN70") + 
  theme(axis.text.x = element_text(size = 4))

p2 <- ggplot(cin_df, aes(x = sample_type, y = cin_bakhoum)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(position = dodge, alpha = 0.2) +
  geom_text(aes(
    x = 1.5, y = max(cin_bakhoum) + 0.1 * max(cin_bakhoum),
    label = "LME p-value = 0.74"
  ), size = 2, family = "Arial") +
  labs(x = "", y = "CIN Bakhoum") + 
  theme(axis.text.x = element_text(size = 4))

p3 <- ggplot(cin_df, aes(x = sample_type, y = HET70)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_point(position = dodge, alpha = 0.2) +
  geom_text(aes(
    x = 1.5, y = max(HET70) + 0.1 * max(HET70),
    label = "LME p-value = 0.43"
  ), size = 2, family = "Arial") +
  labs(x = "", y = "HET70") + 
  theme(axis.text.x = element_text(size = 4))

p <- p1 + p2 + p3 + plot_annotation(tag_levels = 'A')

ggsave(file.path(FIG_DIR, "cin_metastasis_primary.png"), width = 100, height = 50, units = "mm", dpi = 300)
