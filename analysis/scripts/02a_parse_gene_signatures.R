### Get list of gene signatures to run ssGSEA analysis ###

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(msigdbr)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "data", "meta")

# DEFINE GENE SIGNATURES --------------------------------------------------

# 50 MSIGDB HALLMARK SIGNATURES
# get hallmark gene sets using msigdbr
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

# Add also some pathways of interest, such as RNA-Splicing, methylation and epigenetic regulation of gene expression
bp_gene_sets <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  filter(gs_name %in% c("GOBP_RNA_SPLICING", "GOBP_METHYLATION", "GOBP_REGULATION_OF_GENE_EXPRESSION_EPIGENETIC"))

hallmark_gene_sets <- rbind(hallmark_gene_sets, bp_gene_sets)

# wrangle hallmark gene sets for fsgea to work with
msigdbr_H <- split(x = hallmark_gene_sets$gene_symbol, f = hallmark_gene_sets$gs_name)

# MOTZER ET AL (2020) - MOLECULAR SUBSETS
angiogenesis <- "VEGFA, KDR, ESM1, PECAM1, ANGPTL4, CD34"
T_effector <- "CD8A, EOMES, PRF1, IFNG, CD274"
FAO_AMPK <- "CPT2, PPARA, CPT1A, PRKAA2, PDK2, PRKAB1"
cell_cycle <- "CDK2, CDK4, CDK6, BUB1B, CCNE1, POLQ, AURKA, MKI67, CCNB2"
FAS_pentose_phosphate <- "FASN, PARP1, ACACA, G6PD, TKT, TALDO1, PGD"
stroma <- "FAP, FN1, COL5A1, COL5A2, POSTN, COL1A1, COL1A2, MMP2"
myeloid_inflammation <- "CXCL1, CXCL2, CXCL3, CXCL8, IL6, PTGS2"
complement_cascade <- "F2, C1S, C1R, CFB, C3"
omega_oxidation <- "CYP4F3, CYP8B1, NNMT, MGST1, MAOA, CYP4F11, CYP4F2, CYP4F12"
sno_rna <- "SNORD38A, SNORD104, SNORD32A, SNORD68, SNORD66, SNORD100"

motzer_et_al <- list(
  motzer_angiogenesis = unlist(str_split(angiogenesis, ", ")),
  motzer_t_effector = unlist(str_split(T_effector, ", ")),
  motzer_FAO_AMPK = unlist(str_split(FAO_AMPK, ", ")),
  motzer_cell_cycle = unlist(str_split(cell_cycle, ", ")),
  motzer_FAS_pentose_phosphate = unlist(str_split(FAS_pentose_phosphate, ", ")),
  motzer_stroma = unlist(str_split(stroma, ", ")),
  motzer_myeloid_inflammation = unlist(str_split(myeloid_inflammation, ", ")),
  motzer_complement_cascade = unlist(str_split(complement_cascade, ", ")),
  motzer_omega_oxidation = unlist(str_split(omega_oxidation, ", ")),
  motzer_sno_rna = unlist(str_split(sno_rna, ", "))
)

# MCDERMOTT ET AL (2020) - CLINICAL ACTIVITY AND MOLECULAR CORRELATES OF RESPONSE
angiogenesis <- "VEGFA, KDR, ESM1, PECAM1, ANGPTL4, CD34"
T_effector <- "CD8A, EOMES, PRF1, IFNG, CD274"
myeloid_inflammation <- "IL-6, CXCL1, CXCL2, CXCL3, CXCL8, PTGS2"

mc_dermott_et_al <- list(
  mc_dermott_angiogenesis = unlist(str_split(angiogenesis, ", ")),
  mc_dermott_t_effector = unlist(str_split(T_effector, ", ")),
  mc_dermott_myeloid_inflammation = unlist(str_split(myeloid_inflammation, ", "))
)

# MOTZER ET AL (2020) - JAVELIN 101 TRIAL
# immune
t_cell_receptor_signature <- "CD3G, CD3E, CD8B, THEMIS, TRAT1, GRAP2, CD247"
t_cell_activation_proliferation_differentiation <- "CD2, CD96, PRF1, CD6, IL7R, ITK, GPR18, EOMES, SIT1, NLRC3"
nk_cell_mediated_cytotoxicy <- "CD2, CD96, PRF1, CD244, KLRD1, SH2D1A"
chemokine <- "CCL5, XCL2"
other_immune_response_genes <- "CST7, GFI1, KCNA3, PSTPIP1"
# angiogenesis
angiogenesis <- "NRARP	RAMP2	ARHGEF15	VIP NRXN3	KDR	SMAD6	KCNAB1 CALCRL	NOTCH4	AQP1	RAMP3 TEK	FLT1 GATA2 CACNB2 ECSCR GJA5 ENPP2	CASQ2 PTPRB	TBX2 ATP1A2 CD34	HEY2 EDNRB"

javelin <- list(
  J101_immune_26genes = unlist(str_split(c(
    t_cell_receptor_signature,
    t_cell_activation_proliferation_differentiation,
    nk_cell_mediated_cytotoxicy, chemokine,
    other_immune_response_genes
  ), ", ")),
  J101_angiogenesis = unlist(str_split(angiogenesis, "\\s"))
)

# HAKIMI ET AL 2019 COMPARZ ANGIO
angio <- c(
  "CD93", "CDH5", "CLEC14A", "ECSCR", "ELTD1", "EMCN", "ENG", "ESAM", "GPR116", "KDR",
  "LDB2", "MYCT1", "PTPRB", "RHQJ", "ROBO4", "S1PR1", "SPARCL1", "TEK", "TIE1", "VWF"
)
comparz <- list(comparz_angio = angio)

# CIN70 signature
cin70 <- c(
  "TPX2", "PRC1", "FOXM1", "CDC2", "TGIF2", "MCM2", "H2AFZ", "TOOP2A", "PCNA", "UBE2C", "MELK", "TRIP13", "CNAP1", "MCM7", "RNASEH2A",
  "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "TOG", "NUP205", "CDC20",
  "CKS2", "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5", "CDCA8", "PTTG1", "CEP55", "H2AFX",
  "CMAS", "BRRN1", "MCM10", "LSM4", "MTB", "ASF1B", "ZWINT", "TOPK", "FLJ10036", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21",
  "ACTL6A", "GplandMGC13096", "SFRS2", "HDGF", "NXT1", "NEK2", "DHCR7", "STK6", "NDUFAB1", "KIAA0286", "KIF4A"
)
cin70 <- list(cin70 = cin70)


# T-CELL EXHAUSTION SIGNATURES
braun_exhaustion <- c("CD3E", "PDCD1", "IL2RA", "CD69", "ICOS", "CD24", "TNFRSF4", "TNFRSF9", "CD28", "CD27")

annelaure_exhaustion <- c("HAVCR2", "ENTPD1", "PCPD1", "TIGIT", "LAG3", "TOX")
exhaustion <- list(braun_exhaustion = braun_exhaustion, annelaure_exhaustion = annelaure_exhaustion)

# concat to run everything together with gsva function.
all_sign <- c(msigdbr_H, motzer_et_al, mc_dermott_et_al, javelin, comparz, cin70, exhaustion)

saveRDS(all_sign, file.path(OUT_DIR, "gene_signatures.rds"))
