## eQTL of HERV expression with CN at the locus
## correcting by patient of origin, VHL status and purity

rm(list = ls(all = TRUE))

# PACKAGES
library(here)
library(DESeq2)
library(tidyverse)
library(nlme)
library(GenomicRanges)
library(data.table)

# PATHS
BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs", "review")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
CN_DIR <- file.path(BASE, "data", "raw", "cndir")
META_DIR <- file.path(BASE, "data", "meta")
HERV_ANNO_PATH <- file.path(META_DIR, "meta_for_hervs.RDS")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
HERV_PATH <- file.path(BASE, "data", "processed", "herv_counts_corrected.RDS")
HERV_REGIONS_PATH <- file.path(BASE, "data", "meta", "ccrcc_hervs_annotation.tsv")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "analysis", "scripts", "review", "eQTL_functions.R"))

# LOAD DATA
annotation <- read_delim(ANNOTATION_PATH)
hervs_meta <- readRDS(HERV_ANNO_PATH)
hervs <- readRDS(HERV_PATH)
hervs_regions <- read_delim(HERV_REGIONS_PATH)

# Transform to VST to ensure normalised comparison later
smps <- colnames(hervs)
smps <- smps[!grepl("N1t1", smps)]
hervs_meta <- hervs_meta[hervs_meta$sample %in% smps, ]
hervs <- hervs[, colnames(hervs) %in% smps]
dds <- DESeqDataSetFromMatrix(hervs, hervs_meta, ~driverMut)
vst <- varianceStabilizingTransformation(dds, blind = TRUE)
vst <- as.matrix(assay(vst))

# Get matrix of CN per gene
hervs.gr <- GRanges(
    seqnames = str_remove(hervs_regions$seqnames, "chr"),
    IRanges(start = hervs_regions$start, end = hervs_regions$end)
)
hervs.gr$id <- hervs_regions$ID

cnmat <- get_cnmat(smps, CN_DIR, hervs.gr)

# Check expression matrix is OK
all(rownames(vst) == rownames(cnmat))
all(colnames(vst) == colnames(cnmat))

# Get purity vector
pur <- get_pur_vector(smps, annotation)

# Get patient of origin
pts <- get_pat_vector(smps)

is_vhl <- str_detect(hervs_meta$driverMut, "VHL")

# filter low expression
# (i.e. only hervs expressed in at least 70% of tumour samples)
hervs_filtered <- filter_low(hervs, thr = .7)
exp <- vst[rownames(vst) %in% hervs_filtered, ]
cn <- cnmat[rownames(cnmat) %in% hervs_filtered, ]
cn[cn < 1] <- 1
cn[cn > 6] <- 6

lfits <- lapply(1:nrow(exp), function(i) {
    print(i)
    run_lme_eqtl(
        exp = exp[i, ], cn = cn[i, ], pur = pur,
        pt = pts, is_vhl = is_vhl,
        fixed = c("cn", "pur", "is_vhl"), random = "pt"
    )
})

# Get fraction of HERVs where CN significantly explains expression
df <- pull_eqtl_stats(lfits, "cn", rownames(exp))
df$fdr <- p.adjust(df$p, method = "fdr")
nrow(df[df$fdr < .05 & df$t < 0, ])
nrow(df[df$fdr < .05 & df$t > 0, ])
nrow(df[df$p < 0.05 & df$t > 0, ])
nrow(df[df$p < 0.05 & df$t < 0, ])
thr_sig <- min(df$t[df$p < 0.05 & df$t > 0])
thr_sig_down <- max(df$t[df$p < 0.05 & df$t < 0])

p <- ggplot(df, aes(x = t)) +
    geom_histogram(fill = tx_palette[["lightblue"]], alpha = .5) +
    # geom_vline(xintercept = 0, linetype = "dashed", col = "red", size = 1) +
    geom_vline(xintercept = thr_sig, linetype = "dashed", col = tx_palette[["lightred"]], size = 1) +
    geom_vline(xintercept = thr_sig_down, linetype = "dashed", col = tx_palette[["lightpurple"]], size = 1) + 
    labs(x = "LME t-statistic", y = "Number HERVs")

save_ggplot(p, file.path(FIG_DIR, "herv_eqtl_cn"), w = 70, h = 70)

mean(df$t > 0)
nrow(df[p.adjust(df$p, method = "fdr") < 0.05 & df$t > 0, ])
nrow(df[df$p < 0.05 & df$t > 0, ]) / nrow(df)
nrow(df[df$p < 0.05 & df$t < 0, ]) / nrow(df)
