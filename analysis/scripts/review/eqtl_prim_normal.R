## eQTL of DIFFERENCES in expression between two samples
rm(list = ls(all = TRUE))

# PACKAGES
library(here)
library(DESeq2)
library(tidyverse)
library(nlme)
library(GenomicRanges)
library(data.table)
library(MuMIn)
library(biomaRt)
library(RColorBrewer)
library(doMC)
library(foreach)


# PATHS
BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs", "review")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
CN_DIR <- file.path(BASE, "data", "raw", "cndir")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_nor_vst_exp_filtered.rds")
CHROMOSOME_ARMS_PATH <- file.path(BASE, "data", "meta", "chromosome_arms.tsv")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "analysis", "scripts", "review", "eQTL_functions.R"))
source(file.path(BASE, "analysis", "scripts", "review", "eQTL_diff_functions.R"))

# LOAD DATA
annotation <- read_delim(ANNOTATION_PATH)
annotation$sample <- clean_ids(annotation$sample)
annotation <- annotation %>%
    dplyr::select(
        Patient, sample, BAP1,
        SETD2, PBRM1, ARID1A, KDM5C, BAP1, loss_9p, loss_14q,
        Genome.doublings,
        purity, type_collapsed
    )
vst <- assay(readRDS(VST_PATH))
colnames(vst) <- clean_ids(colnames(vst))
# only patients with tumor-normal pairs
normal_samples <- colnames(vst)[grepl("_N1t1", colnames(vst))]
patients <- get_pt(normal_samples)
vst <- vst[, get_pt(colnames(vst)) %in% patients]
chromosome_arms <- read_delim(CHROMOSOME_ARMS_PATH)

# add some dummy information to annotation for normal samples
annotation <- rbind(
    annotation,
    data.frame(
        Patient = get_pt(normal_samples), sample = normal_samples,
        BAP1 = 0, SETD2 = 0, PBRM1 = 0, ARID1A = 0, KDM5C = 0,
        Genome.doublings = 0, loss_9p = 0,
        loss_14q = 0, purity = 0, type_collapsed = "NORMAL"
    )
)


gene_regions <- get_gene_coordinates(rownames(vst))
gene_regions.gr <- GRanges(
    seqnames = gene_regions$chromosome_name,
    IRanges(start = gene_regions$start_position, end = gene_regions$end_position)
)
chr.gr <- GRanges(seqnames = "9", IRanges(
    start = chromosome_arms$start[chromosome_arms$chromosome == "9" & chromosome_arms$arm == "p"],
    end = chromosome_arms$end[chromosome_arms$chromosome == "9" & chromosome_arms$arm == "p"]
))

genes_9p <- findOverlaps(gene_regions.gr, chr.gr)@from %>%
    {
        gene_regions$hgnc_symbol[.]
    } %>%
    unique()
# Get matrix of CN per gene
# genes.gr <- GRanges(
#     seqnames = gene_regions$chromosome_name,
#     IRanges(start = gene_regions$start, end = gene_regions$end)
# )

# genes.gr$id <- gene_regions$hgnc_symbol

# cnmat <- get_cnmat(smps, CN_DIR, genes.gr)
# saveRDS(cnmat, file.path(OUT_DIR, "cnmat.rds"))
cnmat <- readRDS(file.path(OUT_DIR, "cnmat.rds"))
# Assume CN = 2 for normal samples
cnmat_normal <- matrix(2, nrow = nrow(cnmat), ncol = length(normal_samples))
colnames(cnmat_normal) <- normal_samples

cnmat <- cbind(cnmat, cnmat_normal)

# Match to only CN and exp data
genes <- dplyr::intersect(rownames(vst), rownames(cnmat))
samples <- dplyr::intersect(colnames(vst), colnames(cnmat))

vst <- vst[rownames(vst) %in% genes, colnames(vst) %in% samples]
cnmat <- cnmat[rownames(cnmat) %in% genes, colnames(cnmat) %in% samples]

vst <- vst[match(rownames(cnmat), rownames(vst)), ]
vst <- vst[, match(colnames(cnmat), colnames(vst))]
all(rownames(vst) == rownames(cnmat))
all(colnames(vst) == colnames(cnmat))

# Get exp diff matrix
# exp_difs <- get_exp_dif(vst)
# saveRDS(exp_difs, file.path(OUT_DIR, "exp_difs_prim_normal.rds"))
exp_difs <- readRDS(file.path(OUT_DIR, "exp_difs_prim_normal.rds"))
pairs <- colnames(exp_difs)

# Get difference CN
# cn_dif_mat <- get_cn_dif(cnmat, pairs)

# saveRDS(cn_dif_mat, file.path(OUT_DIR, "cn_dif_mat_prim_normal.rds"))
cn_dif_mat <- as.matrix(readRDS(file.path(OUT_DIR, "cn_dif_mat_prim_normal.rds")))
cn_dif_mat <- round(cn_dif_mat)

# Get purity differences vector
pur_dif <- get_pur_dif(annotation, "purity", pairs)

# Get patient of origin
patients <- get_pat_vector(pairs)

# Get differences in driver status
# SETD2
setd2_dif <- get_driver_dif(annotation, "SETD2", pairs)

# BAP1
bap1_dif <- get_driver_dif(annotation, "BAP1", pairs)

# PBRM1
pbrm1_dif <- get_driver_dif(annotation, "PBRM1", pairs)

# ARID1A
arid1_dif <- get_driver_dif(annotation, "ARID1A", pairs)

# KDM5C
kdm5c_dif <- get_driver_dif(annotation, "KDM5C", pairs)

epi_dif <- kdm5c_dif + arid1_dif + pbrm1_dif + bap1_dif + setd2_dif

# 9p loss
chr9p_dif <- get_driver_dif(annotation, "loss_9p", pairs)

# 14q loss
chr14q_dif <- get_driver_dif(annotation, "loss_14q", pairs)

# WGD difference -> to better account for differences in CN
wgd_dif <- get_driver_dif(annotation, "Genome.doublings", pairs)

# Get primary-normal status
annotation$is_normal <- annotation$type_collapsed == "NORMAL"
is_normal <- get_normal_dif(annotation, "is_normal", pairs)

# define random variable to find if contribution sign
random_binary <- sample(0:1, size = length(setd2_dif), replace = T, prob = c(0.5, 0.5))

# define random continuous to find if contribution sign
random_cont <- sample(seq(0, 1, by = 0.05), size = length(setd2_dif), replace = T)

# FIT LME FOR EACH OF THE INDIVIDUAL GENES IN THE MODEL
fixed <- c(
    "cn", "pur", "epi_dif",
    "loss_9p", "random_binary",
    "random_continuous", "is_normal"
)

registerDoMC(5)
lfits <- foreach(i = 1:nrow(exp_difs)) %dopar% {
    # print(i)
    # if gene in 9p, don't add CN to avoid collinearity
    is_9p <- rownames(exp_difs)[i] %in% genes_9p
    tryCatch(
        fit <- run_lme_eqtl(
            exp = scale(exp_difs[i, ]), cn = cn_dif_mat[i, ],
            pur = pur_dif, epi_dif = epi_dif,
            pbrm1 = pbrm1_dif, loss_9p = chr9p_dif,
            loss_14q = chr14q_dif, is_normal = is_normal,
            wgd_dif = wgd_dif,
            random_binary = random_binary,
            pt = patients, random_continuous = random_cont,
            fixed = fixed,
            random = c("pt", "wgd_dif")
        ), error = function(e){print(i);return(NA)}
    )
}
# Get fraction of genes where normal vs tumor explains differences in expression
df_stats <- pull_all_eqtl_stats(lfits, "is_normal", genes = rownames(exp_difs), na.rm = T)


write_delim(df_stats, file.path(OUT_DIR, "normal_eqtl.tsv"), delim = "\t")

df_stats <- read_delim(file.path(OUT_DIR, "normal_eqtl.tsv"), delim = "\t")

mean(df_stats$fdr < .05) # 57.6%
sum(df_stats$fdr < .05) # 8122
mean(df_stats$p < .05) # 61.6%
sum(df_stats$p < .05) # 8676


# Bar plot normal genes
sum(df_stats$fdr <= .05 & df_stats$t > 0)
mean(df_stats$fdr <= .05 & df_stats$t > 0)
sum(df_stats$fdr <= .05 & df_stats$t < 0)
mean(df_stats$fdr <= .05 & df_stats$t < 0)

normal_genes <- table(df_stats$fdr < 0.05 & df_stats$variable == "is_normal")

p <- data.frame(
    total_genes = sum(normal_genes),
    nsignificant = normal_genes[2],
    proportion_significant = normal_genes[2] / sum(normal_genes)
) %>%
    ggplot(aes(x = 1)) +
    geom_col(aes(y = 1), fill = "grey95") +
    geom_col(aes(y = proportion_significant), fill = tx_palette[["lightred"]]) +
    geom_text(aes(y = proportion_significant, label = nsignificant), vjust = -0.5, size = 1.5) +
    ylim(c(0, 1)) +
    labs(y = "Proportion genes significantly dysregulated in primary", x = "")

save_ggplot(p, file.path(FIG_DIR, "Fig_2B"), w = 15, h = 70)