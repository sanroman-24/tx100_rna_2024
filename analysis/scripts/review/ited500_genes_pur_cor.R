## Correlation wtih purity top variable 500 genes used for I-TED calculation

rm(list = ls(all = T))

## PACKAGES
library(here)
library(tidyverse)

## PATHS
BASE <- here::here()
TOP500_PATH <- file.path(BASE, "analysis", "outputs", "top500_variable_genes.rds")
OUT_DIR <- file.path(BASE, "analysis", "outputs", "review")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "src", "get_ited.R"))

# LOAD DATA
annotation <- read_delim(ANNOTATION_PATH)
annotation$sample <- clean_ids(annotation$sample)
vst <- readRDS(VST_PATH)
top500_genes <- readRDS(TOP500_PATH)

# CORRELATION OF INDIVIDUAL GENES WITH PURITY
vst <- vst[, vst$type_collapsed == "PRIMARY"]

purity_correlations <- sapply(top500_genes, function(gene) {
    expression <- as.vector(assay(vst[rownames(vst) == gene]))
    purity <- as.numeric(vst$purity)
    crt <- cor.test(expression, purity, method = "pearson")
    return(c("r" = crt$estimate, "p_value" = crt$p.value))
})

purity_correlations <- as.data.frame(t(purity_correlations))

sum(p.adjust(purity_correlations$p_value) < .05)
sum(purity_correlations$p_value < .05)
mean(purity_correlations$r.cor)
mean(abs(purity_correlations$r.cor))

p <- ggplot(purity_correlations, aes(x = r.cor)) +
    geom_histogram(fill = "lightblue", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", col = "red") +
    geom_vline(xintercept = mean(purity_correlations$r.cor), linetype = "dashed", col = "orchid1") +
    ggtitle("Top 500 variable genes") +
    labs(x = "Pearson correlation with tumor purity", y = "Number of I-TED500 Genes")


purity_correlations_all <- sapply(rownames(vst), function(gene) {
    expression <- as.vector(assay(vst[rownames(vst) == gene]))
    purity <- as.numeric(vst$purity)
    crt <- cor.test(expression, purity, method = "pearson")
    return(c("r" = crt$estimate, "p_value" = crt$p.value))
})

purity_correlations_all <- as.data.frame(t(purity_correlations_all))

p_all <- ggplot(purity_correlations_all, aes(x = r.cor)) +
    geom_histogram(fill = "lightblue", alpha = 0.5) +
    ggtitle("All genes") +
    geom_vline(xintercept = 0, linetype = "dashed", col = "red") +
    geom_vline(xintercept = mean(purity_correlations_all$r.cor), linetype = "dashed", col = "orchid1") +
    labs(x = "Pearson correlation with tumor purity", y = "Number of Genes")

lp <- list(p, p_all)
save_plist(lp, file.path(FIG_DIR, "purity_correlation_ited500_genes"), ncol = 2, w = 100, h = 50)


p1 <- .4
p2 <- .6

simulate_gene_expression <- function(model, purity) {
    coef <- coef(model)
    predicted_value <- coef[1] + coef[2] * purity
    simulated_value <- predicted_value
    return(simulated_value)
}

simulated_expression <- t(sapply(top500_genes, function(gene) {
    print(gene)
    expression <- as.vector(assay(vst[rownames(vst) == gene]))
    purity <- as.numeric(vst$purity)
    mod <- lm(expression ~ purity)
    e1 <- predict(mod, data.frame(purity = p1))
    e2 <- predict(mod, data.frame(purity = p2))
    return(c(e1, e2))
}))

get_dist(simulated_expression[, 1], simulated_expression[, 2], f = "dcor")

p <- data.frame(
    gene = top500_genes,
    e1 = simulated_expression[, 1],
    e2 = simulated_expression[, 2]
) %>%
    ggplot(aes(x = e1, y = e2)) + 
    geom_point(pch = 21, fill = "grey60") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
    geom_smooth(method = "lm") + 
    ggpubr::stat_cor(size = 2) + 
    labs(x = "Expression in R1", y = "Expression in R2")

save_ggplot(p, file.path(FIG_DIR, "purity_ited_mock_r1r2"), w = 50, h = 50)

# I-TED 500 vs purity
iteds_df <- read_delim(
    file.path(BASE, "analysis", "outputs", "ited_correlates.tsv"),
    delim = "\t"
)

p <- ggplot(iteds_df, aes(x = pur_ith)) +
    geom_histogram(fill = "lightblue", alpha = 0.5) +
    labs(y = "Number of patients", x = "Purity ITH") + 
    geom_vline(xintercept = median(iteds_df$pur_ith), linetype = "dashed", col = "red")

save_ggplot(p, file.path(FIG_DIR, "pur_ith_histogram"), w = 50, h = 50)

p <- ggplot(iteds_df, aes(x = pur_ith, y = ited)) +
    geom_point(alpha = 0.5) +
    labs(y = "I-TED", x = "Purity ITH") 

save_ggplot(p, file.path(FIG_DIR, "pur_ith_ited"), w = 50, h = 50)
