# Analysis of the individual association of HERVs with PFS

rm(list = ls(all = TRUE))

# PACKAGES
library(survival)
library(tidyverse)
library(DESeq2)
library(survminer)
library(data.table)

# PATHS ------------------------------------------------------------------------
BASE <- here::here()
DATA <- file.path(BASE, "data", "processed")
META <- file.path(BASE, "data", "meta")
RAW <- file.path(BASE, "data", "raw")
OUT_DIR <- file.path(BASE, "analysis", "outputs")
ANNOTATION_PATH <- file.path(META, "tx_annotation.tsv")
CLINICAL_ANNOTATION_PATH <- file.path(META, "TRACERx_s1_1_clinical.txt")

# FUNCTIONS
run_cox_pfs <- function(df, herv) {
    cox <- coxph(Surv(pfs_time, pfs) ~ e + stage_simple, data = df)
    coef <- summary(cox)$coefficients
    hr <- coef[1,2]
    p_value <- coef[1,5]
    return(data.frame(
        herv = herv, HR = hr, p_value = p_value
    ))
}

# LOAD DATA -------------------------------------------------------------------
hervCounts <- readRDS(file.path(DATA, "herv_counts_corrected.RDS"))
metaDF <- readRDS(file.path(META, "meta_for_hervs.RDS"))
annotation <- read_delim(ANNOTATION_PATH)
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]
clinical_data <- clinical_data %>%
    mutate(pfs = ifelse(`PFS (months)` == "-", 0, 1)) %>%
    mutate(pfs_time = as.numeric(ifelse(`PFS (months)` == "-",
        `Total follow up (months)`, `PFS (months)`
    ))) %>%
    mutate(
        stage_simple =
            ifelse(`Overall Stage` %in% c("I", "II"), "I-II", "III-IV")
    )


# VST normalization to obtain expression matrix
dds <- DESeqDataSetFromMatrix(hervCounts, metaDF, ~driverMut)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsd_mat <- as.matrix(assay(vsd))
vsd_mat <- vsd_mat[, !grepl("N1t1", colnames(vsd_mat))]

# Association of individual HERVs
surv_hervs <- rbindlist(lapply(1:nrow(vsd_mat), function(i) {
    herv_counts <- vsd_mat[i, ]
    df <- data.frame(sample = names(herv_counts), e = herv_counts)
    df$Subject <- str_remove(df$sample, "_.*$")
    df <- df %>%
        group_by(Subject) %>%
        dplyr::summarise(e = mean(e))
    df <- merge(df, clinical_data)
    herv <- rownames(vsd_mat)[i]
    return(run_cox_pfs(df, herv))
}))

write_delim(surv_hervs, file.path(OUT_DIR, "ST8_surv_hervs_tx100.tsv"), delim = "\t")

surv_hervs$padj <- p.adjust(surv_hervs$p_value, "fdr")
table(surv_hervs$p_value < 0.05 & surv_hervs$HR > 1)
table(surv_hervs$p_value < 0.05 & surv_hervs$HR < 1)

quantile(surv_hervs$HR, na.rm = TRUE)


ggplot(surv_hervs, aes(x = HR)) + geom_histogram() + geom_vline(xintercept = 1)
