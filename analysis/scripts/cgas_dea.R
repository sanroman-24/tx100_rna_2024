# Differential expression analysis cGAS-STING genes by wGII

rm(list = ls(all = TRUE))

# LIBRARIES ------------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggthemes)
library(here)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
COUNTS_PATH <- file.path(BASE, "data", "raw", "counts_tumour.RDS")
BAKHOUM_GENES_DIR <- file.path(META_DIR, "bakhoum_genes")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))


# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
counts <- readRDS(COUNTS_PATH)

# GET CGAS-STING GENES -----------------------------------------------------
STING_genes <- c("TMEM173")
NFKB1_genes <- c("NFKB1", "IKBKG", "TRAF6", "RELA")
NFKB2_genes <- c("NFKB2", "RELB", "MAP3K14")
STING_act_genes <- c("CCL5", "CXCL9", "CXCL10", "CXCL11")
STING_supp_genes <- c("MB21D1", "ENPP1", "TREX1", "IFI16", "IRF3", "ATM", "TBK1", "PARP1", "NT5E", "SLC19A1")

# get them into a list format
cGAS_STING_genes <- list(
    STING = STING_genes,
    NFKB1 = NFKB1_genes,
    NFKB2 = NFKB2_genes,
    STING_act = STING_act_genes,
    STING_supp = STING_supp_genes
)

# Add the genes used by Bakhoum's group in one of their CIN papers
bak_signatures_files <- list.files(BAKHOUM_GENES_DIR, pattern = "CIN_.*Bakhoum.txt", full = T)
bak_signatures <- gsub("\\.txt", "", basename(bak_signatures_files))

bak_genes <- lapply(bak_signatures_files, function(filePath) as.vector(read.table(filePath, header = TRUE)[, 1]))
names(bak_genes) <- bak_signatures
bak_genes <- bak_genes[names(bak_genes) != "CIN_InflammationGenes_Bakhoum"]
cGAS_STING_genes <- c(cGAS_STING_genes, bak_genes)


# subset to primaries and match annotation and counts
annotation <- annotation %>% filter(type == "PRIMARY")
counts <- counts[, colnames(counts) %in% annotation$sample]
# check that the samples in the counts matrix are the same to those in the sample annotation
all(colnames(counts) == annotation$sample)

# DIFFERENTIAL EXPRESSION -----------------------------------------------------
annotation$wgii <- as.numeric(annotation$wgii)
annotation$purity <- as.numeric(annotation$purity)
dds <- DESeqDataSetFromMatrix(
    countData = counts, colData = annotation,
    design = ~ purity + wgii
)
dds <- DESeq(dds)
res <- results(dds)

# subset results to only cGAS-STING genes, creating new variable that specifies the "cGAS" group of each gene
res <- res[rownames(res) %in% unlist(cGAS_STING_genes), ] %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_name") %>%
    mutate(gene_cat = case_when(
        gene_name %in% STING_genes ~ "STING_genes",
        gene_name %in% NFKB1_genes ~ "NFKB1_genes",
        gene_name %in% NFKB2_genes ~ "NFKB2_genes",
        gene_name %in% STING_act_genes ~ "STING_act_genes",
        gene_name %in% STING_supp_genes ~ "STING_supp_genes",
        gene_name %in% bak_genes$CIN_InterferonRegulators_Bakhoum ~ "Bakhoum_InterferonRegulators",
        gene_name %in% bak_genes$CIN_NFKB_Targets_Bakhoum ~ "Bakhoum_NFKB_Targets",
        gene_name %in% bak_genes$CIN_NNKFB_Regulators_Bakhoum ~ "Bakhoum_NFKB_Regulators",
        gene_name %in% bak_genes$CIN_NonCanonical_NFKB_Bakhoum ~ "Bakhoum_NonCanonical_NFKB",
        gene_name %in% bak_genes$CIN_Signature_Bakhoum ~ "Bakhoum_CINSignature"
    )) %>% filter(gene_cat != "Bakhoum_CINSignature")

res$log_p <- -log10(res$padj)
write_delim(res, file.path(OUT_DIR, "ST5_cGAS_dea.tsv"), delim = "\t")

res <- res %>% mutate(
        gene_name =
            ifelse(gene_name %in% c(paste0("CXCL", 9:11), "SLC19A1", "ENPP1", "NFKB1"),
                gene_name, ""
            )
    ) %>%
    mutate(alpha = ifelse(gene_name == "", "soft", "dark"))

# plot the results into a Volcano Plot
p <- plot_volcano(
    df = res, x_str = "log2FoldChange",
    y_str = "log_p", fill_str = "gene_cat", lab_str = "gene_name",
    alpha_str = "alpha", lgd = "yes", xlim = c(-5, 5)
)

save_ggplot(p, file.path(FIG_DIR, "Fig4D_dea_wgii_legend"), w = 45, h = 65)


p <- plot_volcano(
    df = res, x_str = "log2FoldChange",
    y_str = "log_p", fill_str = "gene_cat", lab_str = "gene_name",
    alpha_str = "alpha", lgd = "no", xlim = c(-5, 5)
)

save_ggplot(p, file.path(FIG_DIR, "Fig4D_dea_wgii"), w = 45, h = 65)
