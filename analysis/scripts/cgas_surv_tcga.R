# Survival analysis in TCGA-KIRC by:
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
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "tcga_annotation.tsv")
TCGA_WGII_PATH <- file.path(META_DIR, "tcga_kirc_ith_wgii.tsv")
RNA_PATH <- file.path(BASE, "data", "raw", "TCGA_KIRC_counts.rds")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

# LOAD DATA ---------------------------------------------------------------

rnaseq_data <- read_rds(RNA_PATH)
annotation <- read_delim(CLINICAL_ANNOTATION_PATH)
wgii_data <- read_delim(TCGA_WGII_PATH, delim = "\t") %>%
    dplyr::select(Sample, clonal_WGII, subclonal_WGII) %>%
    mutate(total_WGII = clonal_WGII + subclonal_WGII)

# merge tcga clinical and wgii annotation.
# The latter does not have all the TCGA cases, so left join.
annotation <- merge(annotation, wgii_data,
    by.x = "bcr_patient_barcode", by.y = "Sample", all.x = TRUE
)

# Consider only primary solid tumour samples
rnaseq_data <- rnaseq_data[ , rnaseq_data$definition == "Primary solid Tumor"]
rnaseq_data <- rnaseq_data[, !duplicated(rnaseq_data$patient)]

# NORMALISE COUNTS TO VST -------------------------------------------------
counts <- rnaseq_data@assays@data$`HTSeq - Counts`
rownames(counts) <- rnaseq_data@rowRanges$external_gene_name
colnames(counts) <- rnaseq_data@colData$patient
rownames(rnaseq_data@colData) = str_sub(rownames(rnaseq_data@colData), 1, 12)

dds_tcga <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = rnaseq_data@colData, 
                                   design = ~ patient)

vst <- vst(dds_tcga, blind = TRUE)

# WRANGLE DATA ---------------------------------------------------------------
# filter to genes of interest
vst <- assay(vst[rownames(vst) %in% c("ENPP1", "SLC19A1"), ]) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
rownames(vst) <- str_sub(rownames(vst), 1, 12)
vst <- rownames_to_column(vst, "bcr_patient_barcode") 
annotation = merge(vst, annotation, by = "bcr_patient_barcode")

wgii_thresh <- median(annotation$total_WGII, na.rm = T)
annotation <- annotation %>%
    mutate(high_wgii = ifelse(total_WGII > wgii_thresh,
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
annotation$PFS.time = as.numeric(annotation$PFS.time)
annotation$OS.time = as.numeric(annotation$OS.time)

# KAPLAN MEIER PLOTS ---------------------------------------------------------------
plt <- c("tomato3", "grey25", "tomato1", "grey50")

# PFS
km <- survfit(Surv(PFS.time, PFS) ~ SLCENPP_wgii,
    data = annotation
)

km_plot <- plot_km(annotation, km, plt = plt)
save_baseplot(km_plot, file.path(FIG_DIR, "Fig4F_cgas_wgii_pfs"), h = 90, w = 70)

# OS
km <- survfit(Surv(OS.time,OS) ~ SLCENPP_wgii,
    data = annotation
)
km_plot <- plot_km(annotation, km, plt = plt)
save_baseplot(km_plot, file.path(FIG_DIR, "Fig4F_cgas_wgii_os"), h = 90, w = 70)
