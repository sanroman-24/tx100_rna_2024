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
vst <- assay(readRDS(VST_PATH))
colnames(vst) <- clean_ids(colnames(vst))
smps <- annotation$sample[annotation$type_collapsed == "PRIMARY"]
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
# genes.gr <- GRanges(
#     seqnames = gene_regions$chromosome_name,
#     IRanges(start = gene_regions$start, end = gene_regions$end)
# )

# genes.gr$id <- gene_regions$hgnc_symbol

# cnmat <- get_cnmat(smps, CN_DIR, genes.gr)
# saveRDS(cnmat, file.path(OUT_DIR, "cnmat.rds"))
cnmat <- readRDS(file.path(OUT_DIR, "cnmat.rds"))

# Match to only CN and exp data
vst <- vst[match(rownames(cnmat), rownames(vst)), ]
vst <- vst[, colnames(vst) %in% colnames(cnmat)]
all(rownames(vst) == rownames(cnmat))
all(colnames(vst) == colnames(cnmat))

# Get exp diff matrix
# exp_difs <- get_exp_dif(vst)
# saveRDS(exp_difs, file.path(OUT_DIR, "exp_difs_primary.rds"))
exp_difs <- readRDS(file.path(OUT_DIR, "exp_difs_primary.rds"))
pairs <- colnames(exp_difs)

# Get difference CN
# cn_dif_mat <- get_cn_dif(cnmat, pairs)

# saveRDS(cn_dif_mat, file.path(OUT_DIR, "cn_dif_mat_primary.rds"))
cn_dif_mat <- as.matrix(readRDS(file.path(OUT_DIR, "cn_dif_mat_primary.rds")))
cn_dif_mat <- round(cn_dif_mat)


# Get purity differences vector
pur_dif <- get_pur_dif(annotation, "purity", pairs)

# Get patient of origin
patients <- get_pat_vector(pairs)

# Get differences in driver status
# Epigenetic drivers
# KDM5C, ARID1A, SETD2, PBRM1, BAP1

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

# FIT LME FOR EACH OF THE INDIVIDUAL GENES IN THE MODEL

fixed <- c(
    "cn", "pur", "epi_dif",
    "loss_9p", "random_binary",
    "random_continuous"
)

registerDoMC(5)
lfits <- foreach(i = 1:nrow(exp_difs)) %dopar% {
    # if gene in 9p, don't add CN to avoid collinearity
    is_9p = rownames(exp_difs)[i] %in% genes_9p
    if (is_9p){fixed = c("pur", "epi_dif", "loss_9p", "random_binary", "random_continuous")}
    fit <- run_lme_eqtl(
        exp = scale(exp_difs[i, ]), cn = cn_dif_mat[i, ],
        pur = pur_dif, epi_dif = epi_dif, loss_9p = chr9p_dif,
        loss_14q = chr14q_dif, random_binary = random_binary,
        pt = patients, random_continuous = random_cont,
        wgd_dif = wgd_dif, fixed = fixed,
        random = c("pt", "wgd_dif")
    )
}

# find proportion explained by the model
p_explained <- get_proportion_explained(lfits)
quantile(p_explained, na.rm = T)
mean(p_explained, na.rm = T) # 49.24%
t.test(p_explained) # [95%CI = 49.0% - 49.5%]
median(p_explained, na.rm = T) # 50%

p <- data.frame(
    gene = rownames(exp_difs),
    p_explained = p_explained,
    total = 1
) %>%
    arrange(p_explained) %>%
    mutate(gene = factor(gene, levels = gene)) %>%
    ggplot(aes(x = gene)) +
    geom_col(aes(y = total), col = "grey95") +
    geom_col(aes(y = p_explained), col = tx_palette[["lightgreen"]]) +
    labs(x = "Genes", y = "Proportion explained by LME") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

save_ggplot(p, file.path(FIG_DIR, "proportion_prim_td_fitted_model"), w = 100, h = 50)


# Get fraction of genes where CN significantly explains expression
df_stats <- pull_all_eqtl_stats(lfits, fixed, genes = rownames(exp_difs), na.rm = T)
rename <- c(
    "cn" = "CN", "epi_dif" = "Epigenetic driver",
    "loss_14q" = "Loss 14q", "loss_9p" = "Loss 9p", "pur" = "Purity",
    "random_binary" = "Random binary",
    "random_continuous" = "Random continuous"
)

df_stats$variable <- rename[df_stats$variable]
write_delim(df_stats, file.path(OUT_DIR, "eqtl_primary.tsv"), delim = "\t")

df_stats <- read_delim(file.path(OUT_DIR, "eqtl_primary.tsv"), delim = "\t")
df_stats$fdr_model <- p.adjust(df_stats$model_p)

mean(df_stats$fdr_model[!duplicated(df_stats$gene)] < .05)
sum(df_stats$fdr_model[!duplicated(df_stats$gene)] < .05)
mean(df_stats$model_p[!duplicated(df_stats$gene)] < .05)
sum(df_stats$model_p[!duplicated(df_stats$gene)] < .05)


# Plot number of genes significantly explained by each variable
tb <- table(df_stats$fdr <= 0.05, df_stats$variable)

# tb
#           CN Epigenetic driver Loss 9p Purity Random binary Random continuous
#   FALSE  9185              9890    7851   6984         14120             14120
#   TRUE   4791              4230    6269   7136             0                 0

n_total <- tb[1, ] + tb[2, ]
n_explained <- tb[2, ]

p <- data.frame(
    n_total = n_total, n_explained = n_explained,
    variable = names(n_explained)
) %>%
    arrange(n_explained) %>%
    mutate(variable = factor(variable, levels = variable)) %>%
    ggplot(aes(x = variable)) +
    geom_col(aes(y = n_total), fill = "grey95") +
    geom_col(aes(y = n_explained), fill = tx_palette[["lightgreen"]]) +
    labs(x = "Variable in eQTL model", y = "Number significant genes") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

save_ggplot(p, file.path(FIG_DIR, "significant_genes_primary_eqtl"), w = 60, h = 60)

# Plot distributions in violin plot
df_stats$variable <- factor(df_stats$variable,
    levels = c("Purity", "CN", "Loss 9p", "Epigenetic driver", "Random binary", "Random continuous")
)

p <- df_stats %>%
    filter(!is.na(variable)) %>%
    ggplot(aes(x = variable, y = t)) +
    geom_violin(alpha = 0.2, fill = tx_palette[["darkblue"]]) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_fill_manual(values = colors) +
    labs(x = "", y = "Regression coefficient (z-score)") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

save_ggplot(p, file.path(FIG_DIR, "regression_coefficients_eqtl_prim"), w = 45, h = 60)

# Check for trans effects of 9p loss

chrs_9p_genes <- na.omit(
    df_stats$gene[df_stats$variable == "Loss 9p" & df_stats$fdr < .05]
) %>%
    {
        gene_regions$chromosome_name[gene_regions$hgnc_symbol %in% .]
    }

total_genes <- na.omit(
    df_stats$gene[df_stats$variable == "Loss 9p"]
) %>%
    {
        gene_regions$chromosome_name[gene_regions$hgnc_symbol %in% .]
    }

table(chrs_9p_genes) / table(total_genes)

table(chrs_9p_genes != "9")

# 340 with significant difference in expression in chromosome 9 (60%)
# But up to 5950 genes with difference in expression from other chromosomes

chrs_9p_genes <- table(chrs_9p_genes)
total_genes <- table(total_genes)

p <- data.frame(
    chromosome = names(chrs_9p_genes),
    total_genes = as.vector(total_genes),
    nsignificant = as.vector(chrs_9p_genes),
    proportion_significant = as.vector(chrs_9p_genes) / as.vector(total_genes)
) %>%
    ggplot(aes(x = chromosome)) +
    geom_col(aes(y = 1), fill = "grey95") +
    geom_col(aes(y = proportion_significant), fill = tx_palette[["lightgreen"]]) +
    geom_text(aes(y = proportion_significant, label = nsignificant), vjust = -0.5, size = 2) +
    scale_x_discrete(limits = as.character(1:22)) +
    labs(y = "% genes associated with 9p loss", x = "Chromosome")

save_ggplot(p, file.path(FIG_DIR, "SF_chromosomes_significant_9p"), w = 100, h = 50)
