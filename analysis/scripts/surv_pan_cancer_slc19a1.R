# Survival analysis ENPP1 and SLC19A1 pan-cancer

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(survival)
library(survminer)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
## DOWNLOAD CORRESPONDING FILES HERE: https://figshare.com/s/f9982c2ab465317e6cc9
DATA_DIR <- file.path(BASE, "data", "raw", "ext")
META_DIR <- file.path(BASE, "data", "meta")
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
CLINICAL_PATH <- file.path(META_DIR, "allClinicalTCGA.txt")


source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
# SURVIVAL ANALYSIS --------------------------------------------------------

# load clinical data
clin_anno <- read_delim(file = CLINICAL_PATH, delim = "\t")
clin_anno$PFI <- as.numeric(clin_anno$PFI)
clin_anno$OS <- as.numeric(clin_anno$OS)
clin_anno$PFI.time <- as.numeric(clin_anno$PFI.time)
clin_anno$OS.time <- as.numeric(clin_anno$OS.time)

# prepare to parse through all the RNA-Seq files and perform ttype survival analysis
fps <- list.files(DATA_DIR, "RNASeq.RDS")
hrPFSSurvivalLowHigh <- data.frame()
hrOSSurvivalLowHigh <- data.frame()

hrPFSSurvivalLowHighSLCENPP <- data.frame()
hrOSSurvivalLowHighSLCENPP <- data.frame()

hrPFSSurvivalSLC19A1 <- data.frame()
hrOSSurvivalSLC19A1 <- data.frame()

i <- 0
for (fp in fps){
  if (str_detect(fp, "\\w{4}_RNASeq.RDS")){
    i <- i + 1
    tumourType <- str_sub(fp, 1, 4)
    fp <- file.path(DATA_DIR, fp)
    dataset <- read_rds(fp)
    # In TCGA dataset, consider only primary solid tumour samples
    dataset <- dataset[ , dataset$definition == "Primary solid Tumor"]
    dataset <- dataset[, !duplicated(dataset$patient)]
    
    # compute variance stabilizing transformation
    counts <- dataset@assays@data$`HTSeq - Counts`
    rownames(counts) <- dataset@rowRanges$external_gene_name
    colnames(counts) <- dataset@colData$patient

    rownames(dataset@colData) = dataset@colData$patient
    
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = dataset@colData, 
                                  design = ~ patient)
    
    vst <- vst(dds, blind = TRUE)
    
    # prepare to merge to clinical data
    # change data type to data frame and then put genes in columns and patient IDs in rows
    SLC19A1_vst <- assay(vst[rownames(vst) %in% c("SLC19A1", "ENPP1"), ]) %>% 
      as.data.frame() %>% t() %>% as.data.frame()
    # only consider first 12 characters of patient ID as they are the way they are stored in the clinical annotation
    rownames(SLC19A1_vst) <- str_sub(rownames(SLC19A1_vst), 1, 12)
    # move rownames to a new variable to allow merging
    SLC19A1_vst <- rownames_to_column(SLC19A1_vst, "bcr_patient_barcode")
    # merge
    clinical_SLC19A1 <- merge(SLC19A1_vst, clin_anno, by = "bcr_patient_barcode") 
    
    # define high and low SLC19A1 groups
    SLC19A1_thresh <- quantile(clinical_SLC19A1$SLC19A1, 0.75)
    clinical_SLC19A1 <- clinical_SLC19A1 %>% mutate(high_SLC19A1 = ifelse(SLC19A1 > SLC19A1_thresh, "high_SLC19A1", "low_SLC19A1"))
    
    # define high and low ENPP1 groups
    ENPP1_thresh <- quantile(clinical_SLC19A1$ENPP1, 0.75)
    clinical_SLC19A1 <- clinical_SLC19A1 %>% mutate(high_ENPP1 = ifelse(ENPP1 > ENPP1_thresh, "high_ENPP1", "low_ENPP1"))
    
    # define high and low SLCENPP groups
    clinical_SLC19A1 <- clinical_SLC19A1 %>% mutate(high_SLCENPP = ifelse((high_ENPP1 == "high_ENPP1" & high_SLC19A1 == "high_SLC19A1"), 
                                                    "high_SLCENPP1", "low_SLCENPP1"))
    
    # Cox univariate analysis between binary classification based on SLC19A1
    coxUVABinary <- coxph(Surv(as.numeric(PFI.time), PFI) ~ high_SLC19A1,
                          data=clinical_SLC19A1)
    resultsBinaryCox <- data.frame(tumour = tumourType, CI2.5 = exp(confint(coxUVABinary)[1]), CI97.5 = exp(confint(coxUVABinary)[2]))
    hrPFSSurvivalLowHigh <- rbind(hrPFSSurvivalLowHigh, resultsBinaryCox)
    
    coxUVABinary <- coxph(Surv(as.numeric(OS.time), OS) ~ high_SLC19A1,
                          data=clinical_SLC19A1)
    resultsBinaryCox <- data.frame(tumour = tumourType, CI2.5 = exp(confint(coxUVABinary)[1]), CI97.5 = exp(confint(coxUVABinary)[2]))
    hrOSSurvivalLowHigh <- rbind(hrOSSurvivalLowHigh, resultsBinaryCox)
    
    # Cox univariate analysis between binary classification based on SLC19A1 & ENPP1
    coxUVABinary <- coxph(Surv(as.numeric(PFI.time), PFI) ~ high_SLCENPP,
                          data=clinical_SLC19A1)
    resultsBinaryCox <- data.frame(tumour = tumourType, CI2.5 = exp(confint(coxUVABinary)[1]), CI97.5 = exp(confint(coxUVABinary)[2]))
    hrPFSSurvivalLowHighSLCENPP <- rbind(hrPFSSurvivalLowHighSLCENPP, resultsBinaryCox)
    
    coxUVABinary <- coxph(Surv(as.numeric(OS.time), OS) ~ high_SLCENPP,
                          data=clinical_SLC19A1)
    resultsBinaryCox <- data.frame(tumour = tumourType, CI2.5 = exp(confint(coxUVABinary)[1]), CI97.5 = exp(confint(coxUVABinary)[2]))
    hrOSSurvivalLowHighSLCENPP <- rbind(hrOSSurvivalLowHighSLCENPP, resultsBinaryCox)
    
    # Cox univariate analysis with continuous SLC19A1 expression
    coxUVAContinuous <- coxph(Surv(as.numeric(PFI.time), PFI) ~ SLC19A1,
                          data=clinical_SLC19A1)
    resultsContinuousCox <- data.frame(tumour = tumourType, CI2.5 = exp(confint(coxUVAContinuous)[1]), CI97.5 = exp(confint(coxUVAContinuous)[2]))
    hrPFSSurvivalSLC19A1 <- rbind(hrPFSSurvivalSLC19A1, resultsContinuousCox)
    
    coxUVAContinuous <- coxph(Surv(as.numeric(OS.time), OS) ~ SLC19A1,
                          data=clinical_SLC19A1)
    resultsContinuousCox <- data.frame(tumour = tumourType, CI2.5 = exp(confint(coxUVAContinuous)[1]), CI97.5 = exp(confint(coxUVAContinuous)[2]))
    hrOSSurvivalSLC19A1 <- rbind(hrOSSurvivalSLC19A1, resultsContinuousCox)
  }
}

p <- hrOSSurvivalSLC19A1 %>% mutate(CI50 = (CI2.5+CI97.5) / 2) %>%
  mutate(col = ifelse(CI2.5 >= 1, "red", "grey")) %>%  
  ggplot(aes(y = reorder(tumour, CI2.5), xmin = CI2.5, xmax = CI97.5, x = CI50, col = col)) + 
  geom_vline(xintercept = 1, lty = "dashed") + 
  scale_color_manual(values = c("lightgrey", "firebrick3")) + 
  geom_pointrange() + 
  theme_classic() + 
  xlim(c(0, 5)) + 
  theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(strip.background = element_rect(colour="white")) + 
  theme(legend.position = "none") + 
  labs(x = "Hazard Ratio", y = "") + 
  scale_x_log10(breaks = c(0.1, 1, 10), limits = c(0.1, 10)) 

# save plot
save_ggplot(p, file.path(FIG_DIR, "SupFig16_pancan_os_slc19a1"), w = 100, h = 100)

# plot forest plot for PFS Hazard Ratio with SLC19A1 as continuous variable
p <- hrPFSSurvivalSLC19A1 %>% mutate(CI50 = (CI2.5+CI97.5) / 2) %>%
  mutate(col = ifelse(CI2.5 >= 1, "red", "grey")) %>%  
  ggplot(aes(y = reorder(tumour, CI2.5), xmin = CI2.5, xmax = CI97.5, x = CI50, col = col)) + 
  geom_vline(xintercept = 1, lty = "dashed") + 
  scale_color_manual(values = c("lightgrey", "firebrick3")) + 
  geom_pointrange() + 
  theme_classic() + 
  theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(strip.background = element_rect(colour="white")) + 
  theme(legend.position = "none") + 
  labs(x = "Hazard Ratio", y = "") + 
  scale_x_log10(breaks = c(0.1, 1, 10), limits = c(0.1, 10)) 

# save plot
save_ggplot(p, file.path(FIG_DIR, "SupFig16_pancan_pfs_slc19a1"), w = 100, h = 100)
