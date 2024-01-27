# Use LME to associate PCs from TRACERx Renal with clinical, genomic and TME/transcriptional variables

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(here)
library(nlme)
library(ggpubr)
library(lemon)



# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
TX_PCA_PATH <- file.path(BASE, "analysis", "outputs", "tx_transcriptional_pca.rds")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
SSGSEA_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")

# FUNCTIONS ---------------------------------------------------------------

source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "utils.R")) # contains run_lme function

# LOAD DATA ---------------------------------------------------------------

tx_pca <- readRDS(TX_PCA_PATH)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]
clinical_data$Patient <- clinical_data$Subject
ssgsea <- read_delim(SSGSEA_PATH, delim = "\t")



# PREPARE DATA TO RUN LME -------------------------------------------------


pcs <- tx_pca$x[, c(1, 2, 3, 4, 5)] # in a previous analysis 43% of the total variance was explained by the 5 first PC components

all_annotation <- merge(annotation, ssgsea, by = "sample") %>% left_join(clinical_data, by = "Patient")
all_annotation$purity <- as.numeric(all_annotation$purity)
pcs_annotated <- pcs %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  merge(all_annotation, by = "sample")
pcs_annotated$Age <- as.numeric(pcs_annotated$Age)
pcs_annotated$stage <- pcs_annotated$`Overall Stage`
pcs_annotated$wgii <- as.numeric(pcs_annotated$wgii)
pcs_annotated$SETD2 <- ifelse(pcs_annotated$SETD2 > 0, 1, 0)
pcs_annotated$VHL <- ifelse(pcs_annotated$VHL > 0, 1, 0)
pcs_annotated$PBRM1 <- ifelse(pcs_annotated$PBRM1 > 0, 1, 0)
pcs_annotated$KDM5C <- ifelse(pcs_annotated$KDM5C > 0, 1, 0)
pcs_annotated$MTOR <- ifelse(pcs_annotated$MTOR > 0, 1, 0)
pcs_annotated$WGD <- pcs_annotated$Genome.doublings
pcs_annotated$is_primary <- ifelse(pcs_annotated$type_collapsed == "PRIMARY", 1, 0)
pcs_annotated <- pcs_annotated[!is.na(pcs_annotated$purity), ]


vars <- c(
  "Age", "Sex", "stage", "is_primary", "purity",
  "SETD2", "PBRM1", "VHL", "BAP1", "MTOR", "KDM5C",
  "wgii", "WGD", "loss_9p", "loss_14q"
)

res_df <- data.frame(feature = c(), pc = c(), pval = c(), test = c())

for (var in vars) {
  for (pc in paste0("PC", 1:5)) {
    # correction by tumour purity, unless we are testing purity itself
    lme_res <- summary(run_lme(pc, c(var, "purity"), "Patient", pcs_annotated))
    if (var == "purity") {
      lme_res <- summary(run_lme(pc, "purity", "Patient", pcs_annotated))
    }

    pval <- lme_res$tTable[2, 5]
    test <- lme_res$tTable[2, 4]
    res_df <- rbind(res_df, data.frame(feature = var, pc = pc, pval = pval, test = test))
  }
}

res_df$padj <- p.adjust(res_df$pval)

p <- ggplot(res_df, aes(x = pc, y = feature, fill = ifelse(padj < .05, -log10(padj) * sign(test), NA))) +
  geom_tile(size = .3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "firebrick", midpoint = 0, na.value = "white") +
  scale_y_discrete(limits = rev(c("Age", "Sex", "is_primary", "stage", "purity", "VHL", "BAP1", "SETD2", "PBRM1", "KDM5C", "MTOR", "wgii", "WGD", "loss_9p", "loss_14q"))) +
  labs(x = "", y = "", fill = "-log10 FDR * sign(test)") +
  theme(legend.position = "none")

p <- change_axes(p)

save_ggplot(p, file.path(OUT_DIR, "Fig2C_pca_lme_no_legend"), w = 50, h = 40)

p <- ggplot(res_df, aes(x = pc, y = feature, fill = ifelse(padj < .05, -log10(padj) * sign(test), NA))) +
  geom_tile(size = .3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "firebrick", midpoint = 0, na.value = "white") +
  scale_y_discrete(limits = rev(c("Age", "Sex", "is_primary", "stage", "purity", "VHL", "BAP1", "SETD2", "PBRM1", "KDM5C", "MTOR", "wgii", "WGD", "loss_9p", "loss_14q"))) +
  labs(x = "", y = "", fill = "-log10 FDR * sign(test)")

p <- change_axes(p)

save_ggplot(p, file.path(OUT_DIR, "Fig2C_pca_lme"), w = 50, h = 40)
