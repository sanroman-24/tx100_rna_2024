#-----------------------------------------------------------------------------#
#                              HERV ANALYSIS                                  #      
#-----------------------------------------------------------------------------#

set.seed(920101)

# GENERAL LIBRARIES -----------------------------------------------------------
lib <- c("tidyverse", "Polychrome", "RColorBrewer", "patchwork", "janitor",
         "ComplexHeatmap", "ggpubr", "ggrepel")
invisible(lapply(lib, library, character.only = TRUE))

# COLOR PALETTE ---------------------------------------------------------------
basicPal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
              "#F781BF", "#8DD3C7", "#224F63", "#FFFF33", "#999999")

# PATHS ------------------------------------------------------------------------
BASE <- here::here()
DATA <- file.path(BASE, "data", "processed")
META <- file.path(BASE, "data", "meta")
RAW <- file.path(BASE, "data", "raw")

# LOAD DATA -------------------------------------------------------------------
hervCounts <- readRDS(file.path(DATA, "herv_counts_corrected.RDS"))
metaDF <- readRDS(file.path(META, "meta_for_hervs.RDS"))

#------------------------------------#
#   UMAP ANALYSIS (SUPP FIGURE 23)   #
#------------------------------------#

library(DESeq2)
library (umap)

# VST normalization to obtain expression matrix
dds <- DESeqDataSetFromMatrix(hervCounts, metaDF, ~driverMut)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsd_mat <- as.matrix(assay(vsd))

# Fit UMAP
vsd_umap <- umap(t(vsd_mat))
umapDF <- as.data.frame(vsd_umap$layout) %>% 
  rownames_to_column(var = "sample")
colnames(umapDF) <- c("sample", "UMAP_1", "UMAP_2")
umapDF <- left_join(umapDF, metaDF, by = "sample")

# Mark the patient IDs similar to Figure 2A
points_to_plot <- umapDF %>% 
  filter(Patient %in% c("K243", "K390", "K153"))

ggplot() +
  geom_point(umapDF, mapping = aes(x = UMAP_1, y = UMAP_2), colour = "grey80") +
  theme_classic() +
  geom_point(points_to_plot, mapping = aes(x = UMAP_1, y = UMAP_2, colour = Patient), size = 2) +
  scale_color_manual(values = c("darkgreen", "red", "blue"))

#------------------------------------------#
#   EXPRESSION ANALYSIS (FIGURE 7A & 7B)   #
#------------------------------------------#

# Get median HERV expression per sample (Fig 7A)
vsd_median <- colMedians(vsd_mat) %>% 
  as.data.frame() %>% 
  rownames_to_column()
colnames(vsd_median) <- c("sample", "medianExp")
vsd_median$sample <- colnames(vsd_mat)

vsd_median_plot <- vsd_median %>% 
  left_join(metaDF, by = "sample") %>% 
  dplyr::filter(!driverType %in% c("Multiple", "SETD2"))

driverLevels <- c("Adjacent normal", "VHL (WT)", "VHL",  "VHL + PBRM1", "VHL + PBRM1 + SETD2",  
                  "VHL + SETD2", "VHL + BAP1", "VHL(M)", "VHL(M) + PBRM1", "VHL(M) + PBRM1 + SETD2", 
                  "VHL(M) + SETD2", "VHL(M) + BAP1")
vsd_median_plot$driverType <- factor(vsd_median_plot$driverType, levels = driverLevels)

ggplot(vsd_median_plot, aes(x = driverType, y = medianExp, colour = VHL_methyl)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) + 
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 1) +
  theme_minimal() +
  labs(x = "", y = "Median expression") +
  coord_flip() +
  scale_color_manual(values = c("steelblue", "darkred"))

# Differential expression analysis (Fig 7B)
# Make new DESeq2 object with new metadata annotation to compare normal & VHL (WT) to tumour samples with driver mutations
meta_for_DE <- metaDF %>%
  filter(!driverType %in% c("Multiple", "SETD2")) %>% 
  mutate(DEgrp = ifelse(driverType %in% c("Adjacent normal", "VHL (WT)"), "NorWT", "VHL"))
hervCounts_for_DE <- hervCounts[, meta_for_DE$sample]

dds <- DESeqDataSetFromMatrix(hervCounts_for_DE, meta_for_DE, ~DEgrp)
dds <- DESeq(dds)

# Get results
resDE <- results(dds, alpha = 0.05)
resDE_lfc_normal <- lfcShrink(dds, res = resDE, type = "normal", coef = "DEgrp_VHL_vs_NorWT")

# Plot results
plotDF <- as.data.frame(resDE_lfc_normal) %>% 
  mutate(log_padj = -log(padj, base = 10)) %>% 
  mutate(colGrp = ifelse(log2FoldChange > 2 & log_padj > 15, "A", 
                         ifelse(log2FoldChange < -2 & log_padj > 15, "B", "C"))) %>% 
  drop_na()

labelDF <- plotDF %>% filter(colGrp == "A") %>% mutate(ID = rownames(.))

ggplot(plotDF, aes(x = log2FoldChange, y = log_padj, col = colGrp)) +
  geom_point(size = 3) +
  theme_minimal() +
  geom_hline(yintercept = 15, linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  scale_color_manual(values = alpha(c("red", "grey40"), alpha = 0.5)) +
  labs(x = bquote(~Log[2]~ "fold change"), y = bquote("-" ~Log[10]~ "p-adj")) +
  theme(legend.position = "none") +
  geom_text_repel(data = labelDF, mapping = aes(x = log2FoldChange, y = log_padj, label = ID), color = "red")

#-----------------------------------#
#   SURVIVAL ANALYSIS (FIGURE 7C)   #
#-----------------------------------#

# Check correlation between median HERV expression with clinical outcomes
vsd_km <- vsd_median %>% 
  left_join(metaDF, by = "sample") %>% 
  dplyr::filter(!driverType %in% c("Adjacent normal", "VHL (WT)")) %>% 
  group_by(Patient) %>% 
  mutate(medianExpOverall = median(medianExp)) %>%       # Get median expression per patient
  ungroup() %>% 
  distinct(Patient, .keep_all = TRUE) %>% 
  mutate(pfs = ifelse(PFS_months == "-", 0, 1)) %>% 
  mutate(pfs_time = as.numeric(ifelse(PFS_months == "-", Total_follow_up_months, PFS_months))) %>% 
  mutate(OS = ifelse(Outcome == "Death", 1, 0)) %>% 
  mutate(OS_time = as.numeric(Total_follow_up_months)) %>% 
  mutate(highHERV = medianExpOverall > median(medianExpOverall)) %>%     # Classify into high / low HERV expression
  mutate(stratify = case_when(highHERV == TRUE & Stage %in% c("I", "II") ~ "High HERV + Low Stage",
                              highHERV == TRUE & Stage %in% c("III", "IV") ~ "High HERV + High Stage",
                              highHERV == FALSE & Stage %in% c("I", "II") ~ "Low HERV + Low Stage",
                              highHERV == FALSE & Stage %in% c("III", "IV") ~ "Low HERV + High Stage"))

# Univariate analysis
library(survival)
library(survminer)
library(forestmodel)

km_fit <- survfit(Surv(pfs_time, pfs) ~ stratify, data = vsd_km)
km_plot <- ggsurvplot(km_fit, data = vsd_km, pval = TRUE, conf.int = FALSE, risk.table = TRUE, risk.table.col = "strata",
                      palette = basicPal, legend = "none", tables.y.text = FALSE)
km_plot$plot <- km_plot$plot + labs(x = "", y = "Progression Free Survival")
km_plot$table <- km_plot$table + labs(x = "Months", y = "")
km_plot

# Multivariate analysis using cox proportional hazard model
vsd_mva <- vsd_km %>% 
  mutate(Stage = case_when(Stage %in% c("I", "II") ~ "Low",
                           Stage %in% c("III", "IV") ~ "High")) %>% 
  mutate(Stage = factor(Stage, levels = c("Low", "High"))) %>% 
  transmute(pfs_time, pfs, medianExpOverall, purity, Stage)

mvaRes <- coxph(Surv(pfs_time, pfs) ~ ., data = vsd_mva)
forest_model(mvaRes)

#-----------------------------------------------------------#
#   CHECK RELATIONSHIP WITH TCR/BCR DIVERSITY (FIGURE 7D)   #
#-----------------------------------------------------------#

library(immunarch)

# Check correlation with of median HERV expression per patient with TCR diversity calculated using Gini index
tcr <- readRDS(file.path(DATA, "TCR_beta_for_hervs.RDS"))
gini <- as.data.frame(repDiversity(tcr$data, "gini"))
gini$Sample <- rownames(gini)
gini <- left_join(gini, tcr$meta, by = "Sample")
colnames(gini) <- c("Gini", "sample", "Total_clone_count", "Patient", "Type", "WGII", "Purity", "WGII_group")
gini <- gini[,1:3]

vsd_immRep <- vsd_median %>% 
  left_join(metaDF, by = "sample") %>% 
  dplyr::filter(!driverType %in% c("Adjacent normal", "VHL (WT)")) %>% 
  left_join(gini, by = "sample")

ggplot(vsd_immRep, aes(x = medianExp, y = Gini)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_minimal() +
  labs(x = "Median HERV expression", y = "TCR diversity (Gini coefficient)")

# Check correlation with of median HERV expression per patient with BCR diversity calculated using Gini index
bcr_normal <- readRDS(file.path(RAW, "normal_BCR_clones.RDS"))
bcr_tumor <- readRDS(file.path(RAW, "tumour_BCR_clones.RDS"))

gini_normal <- as.data.frame(repDiversity(bcr_normal$data, "gini"))
gini_normal <- gini_normal %>% 
  mutate(ID = rownames(.)) %>% 
  mutate(sample = str_split_i(ID, "\\.", 1)) %>% 
  mutate(sample = str_replace(sample, "[GRT]_", "")) %>% 
  dplyr::rename(Gini = V1) %>% 
  dplyr::select(sample, Gini)

gini_tumor <- as.data.frame(repDiversity(bcr_tumor$data, "gini"))
gini_tumor <- gini_tumor %>% 
  mutate(ID = rownames(.)) %>% 
  mutate(sample = str_split_i(ID, "\\.", 1)) %>% 
  mutate(sample = str_replace(sample, "[GRT]_", "")) %>% 
  mutate(sample = str_replace(sample, "_r[123]", "")) %>% 
  mutate(sample = str_replace(sample, "-", "_")) %>% 
  dplyr::rename(Gini = V1) %>% 
  dplyr::select(sample, Gini)

gini_bcr <- rbind(gini_normal, gini_tumor) %>% 
  filter(sample %in% vsd_immRep$sample) %>% 
  dplyr::rename(GiniBCR = Gini)

vsd_immRep <- vsd_immRep %>% 
  left_join(gini_bcr, by = "sample")

ggplot(vsd_immRep, aes(x = medianExp, y = GiniBCR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_minimal() +
  labs(x = "Median HERV expression", y = "BCR diversity (Gini coefficient)")

#--------------------------------------------------#
#   CHECK RELATIONSHIP WITH TME (SUPP FIGURE 24)   #
#--------------------------------------------------#

# Remove normal & WT VHL to keep background constant
vsd_cor <- vsd_median %>% 
  left_join(metaDF, by = "sample") %>% 
  dplyr::filter(!driverType %in% c("Adjacent normal", "VHL (WT)")) %>% 
  dplyr::select(sample, medianExp, purity)

# Read in TME results & keep only samples in vsd_cor
tme <- readRDS("data/processed/tum_consensustme.rds")
colnames(tme) <- str_replace_all(colnames(tme), "-", "_")
keep <- dplyr::intersect(vsd_cor$sample, colnames(tme))
tme <- tme[, keep]
idx <- match(colnames(tme), vsd_cor$sample)
vsd_cor <- vsd_cor[idx, ]
all(vsd_cor$sample == colnames(tme))
vsd_cor <- t(vsd_cor) %>% as.data.frame() %>% janitor::row_to_names(1)
all(colnames(vsd_cor) == colnames(tme))

# Plot correlation for known ccRCC-associated HERVs
erve4 <- c("trans439b87486cb5eaf5", "transb9be6f94f4ef3f52")
herv4700 <- c("trans90fe87ed60d8122f")

vsd_known_hervs <- vsd_mat[c(erve4, herv4700), colnames(tme)] %>% as.data.frame()
idx <- match(colnames(tme), colnames(vsd_known_hervs))
vsd_known_hervs <- vsd_known_hervs[, idx]
all(colnames(vsd_known_hervs) == colnames(tme))

# Bind all into 1 data frame
tme <- rbind(tme, vsd_cor, vsd_known_hervs)
tme_for_corr <- t(tme) %>% as.data.frame() %>% mutate_all(function(x) as.numeric(x))

# Plot correlogram
library(corrr)
library(corrplot)
library(Hmisc)

corRes <- rcorr(as.matrix(tme_for_corr), type = "spearman")

rows <- c("medianExp", "trans439b87486cb5eaf5", "transb9be6f94f4ef3f52", "trans90fe87ed60d8122f")
cols <- c("B_cells", "Cytotoxic_cells", "Dendritic_cells", "Endothelial", "Eosinophils", "Fibroblasts", "Macrophages",
          "Macrophages_M1", "Macrophages_M2", "Mast_cells", "Monocytes", "NK_cells", "Neutrophils", "Plasma_cells", 
          "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "T_regulatory_cells", "Immune_Score", "purity")

corRes$r <- corRes$r[rows, cols]
corRes$n <- corRes$n[rows, cols]
corRes$P <- corRes$P[rows, cols]

corrplot(corRes$r, p.mat = corRes$P, method = "color", is.corr = TRUE, sig.level = c(0.001, 0.01, 0.05), pch.cex = 1, 
         insig = "label_sig", col = COL2("PiYG", 40), bg = "black")



