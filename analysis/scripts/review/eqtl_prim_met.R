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
library(foreach)
library(doMC)


# PATHS
BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs", "review")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
CN_DIR <- file.path(BASE, "data", "raw", "cndir")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")
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
smps <- annotation$sample
vst <- assay(readRDS(VST_PATH))
colnames(vst) <- clean_ids(colnames(vst))
chromosome_arms <- read_delim(CHROMOSOME_ARMS_PATH)

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
    {gene_regions$hgnc_symbol[.]} %>% unique()
# Get matrix of CN per gene
genes.gr <- GRanges(
    seqnames = gene_regions$chromosome_name,
    IRanges(start = gene_regions$start, end = gene_regions$end)
)

genes.gr$id <- gene_regions$hgnc_symbol

# cnmat <- get_cnmat(smps, CN_DIR, genes.gr)
# saveRDS(cnmat, file.path(OUT_DIR, "cnmat_prim_mets.rds"))
cnmat <- readRDS(file.path(OUT_DIR, "cnmat_prim_mets.rds"))

# Match to only CN and exp data
vst <- vst[match(rownames(cnmat), rownames(vst)), ]
vst <- vst[, colnames(vst) %in% colnames(cnmat)]
all(rownames(vst) == rownames(cnmat))
all(colnames(vst) == colnames(cnmat))

# Get exp diff matrix
# exp_difs <- get_exp_dif(vst)
# saveRDS(exp_difs, file.path(OUT_DIR, "exp_difs_primary_and_mets.rds"))
exp_difs <- readRDS(file.path(OUT_DIR, "exp_difs_primary_and_mets.rds"))
pairs <- colnames(exp_difs)

# Get difference CN
# cn_dif_mat <- get_cn_dif(cnmat, pairs)

# saveRDS(cn_dif_mat, file.path(OUT_DIR, "cn_dif_mat_primary_and_mets.rds"))
cn_dif_mat <- as.matrix(readRDS(file.path(OUT_DIR, "cn_dif_mat_primary_and_mets.rds")))
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

# define random variable to find if contribution sign
random_binary <- sample(0:1, size = length(setd2_dif), replace = T, prob = c(0.5, 0.5))

# define random continuous to find if contribution sign
random_cont <- sample(seq(0, 1, by = 0.05), size = length(setd2_dif), replace = T)

# define random variable to define prim-met/prim-prim/met-met pairs
annotation$is_metastasis <- annotation$type_collapsed == "METASTASIS"
is_metastasis <- get_met_dif(annotation, "is_metastasis", pairs)

# FIT LME FOR EACH OF THE INDIVIDUAL GENES IN THE MODEL
fixed <- c(
    "cn", "pur", "epi_dif",
    "loss_9p", "is_metastasis",
    "random_binary",
    "random_continuous"
)

registerDoMC(5)
lfits <- foreach(i = 1:nrow(exp_difs)) %dopar% {
    print(i)
    # if gene in 9p, don't add CN to avoid collinearity
    is_9p = rownames(exp_difs)[i] %in% genes_9p
    fit <- run_lme_eqtl(
        exp = scale(exp_difs[i, ]), cn = cn_dif_mat[i, ],
        pur = pur_dif, epi_dif = epi_dif, loss_9p = chr9p_dif,
        loss_14q = chr14q_dif, is_metastasis = is_metastasis,
        random_binary = random_binary, wgd_dif = wgd_dif,
        pt = patients, random_continuous = random_cont,
        fixed = fixed,
        random = c("pt", "wgd_dif")
    )
}


# Get fraction of genes where prim-met explains expression
df_stats <- pull_all_eqtl_stats(lfits, "is_metastasis", genes = rownames(exp_difs), na.rm = T)

write_delim(df_stats, file.path(OUT_DIR, "metastasis_eqtl.tsv"), delim = "\t")

df_stats <- read_delim(file.path(OUT_DIR, "metastasis_eqtl.tsv"), delim = "\t")

mean(df_stats$fdr < .05) # 34.6%
sum(df_stats$fdr < .05) # 4886
mean(df_stats$p < .05) # 44.44%
sum(df_stats$p < .05) # 6276

# Bar plot metastasis genes
sum(df_stats$fdr <= .05 & df_stats$t > 0)
mean(df_stats$fdr <= .05 & df_stats$t > 0)
sum(df_stats$fdr <= .05 & df_stats$t < 0)
mean(df_stats$fdr <= .05 & df_stats$t < 0)
met_genes <- table(df_stats$fdr < 0.05 & df_stats$variable == "is_metastasis")

p <- data.frame( 
    total_genes = sum(met_genes), 
    nsignificant = met_genes[2],
    proportion_significant = met_genes[2] / sum(met_genes)
) %>% 
    ggplot(aes(x = 1)) + 
    geom_col(aes(y = 1), fill = "grey95") + 
    geom_col(aes(y = proportion_significant), fill = tx_palette[["lightred"]]) + 
    geom_text(aes(y = proportion_significant, label = nsignificant), vjust = -0.5, size = 1.5) + 
    ylim(c(0,1)) + 
    labs(y = "Proportion genes significantly dysregulated in metastasis", x = "")

save_ggplot(p, file.path(FIG_DIR, "Fig_2D"), w = 15, h = 70)
