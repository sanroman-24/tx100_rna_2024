### Generate an oncoprint representing TRACERx Renal RNA-Seq cohort

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(this.path)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
# Clinical annotation from TRACERx Renal, Cell 2018 (in Sup Table 1)
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
# Spatial annotation from Zhao paper
SPATIAL_INFO_FPATH <- file.path(META_DIR, "Zhao_2021_spatial_annotation.tsv")


# FUNCTIONS ---------------------------------------------------------------

get_n_subclonal <- function(annotation, driver) {
  alt <- ifelse(annotation[[driver]] == 0, 0, 1)
  tab <- table(alt, annotation$Patient)
  n_subclonal <- sum(tab[1, ] > 0 & tab[2, ] > 0)
  return(n_subclonal)
}

source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "plotting_theme.R"))
# LOAD DATA ---------------------------------------------------------------

annotation <- read_delim(ANNOTATION_PATH)
annotation$sample <- str_replace(annotation$sample, "-", "_")
annotation <- annotation[annotation$sample != "K207-R4", ]

# add clinical information
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]
clinical_data$Patient <- clinical_data$Subject
annotation <- left_join(annotation, clinical_data, by = "Patient")
annotation$stage <- annotation$`Overall Stage`

# add spatial information
spatial_annotation <- read_delim(SPATIAL_INFO_FPATH)
spatial_annotation$sample <- paste0(spatial_annotation$Patient, "_", spatial_annotation$Region)
spatial_annotation <- spatial_annotation %>% dplyr::select(sample, X, Y, `distance to boundary (mm)`)
annotation <- left_join(annotation, spatial_annotation, by = "sample")
annotation$distance_to_boundary <- annotation$`distance to boundary (mm)`

# WRANGLE DATA TO PLOT ----------------------------------------------------

annotation$type_collapsed <- str_replace(str_to_title(annotation$type_collapsed), "_", " ")

annotation$wgii <- as.numeric(annotation$wgii)
annotation$ITH <- as.numeric(annotation$ITH)
annotation$ITH <- log(annotation$ITH + 1)
annotation$purity <- as.numeric(annotation$purity)
annotation$ploidy <- as.numeric(annotation$ploidy)

# pull out relevant features for plotting the oncoprint
info_oncoprint <- annotation %>% dplyr::select(
  Patient, sample, # to have IDs
  VHL, PBRM1, SETD2, BAP1,
  loss_3p, loss_9p, loss_14q,
  purity, ploidy, wgii, stage, ITH,
  type_collapsed
) %>% # annotation OncoPrint
  dplyr::rename(Sample = sample)

# transpose data so that columns are patients and rows features
oncoprint_mat <- info_oncoprint %>% t()
# store sample IDs as the column names of the matrix that will be plot
colnames(oncoprint_mat) <- oncoprint_mat[2, ]
# store in a vector the patient IDs to use them to split oncoprint for easier visualization
patients <- oncoprint_mat[1, ]
# remove sample and patient IDs from the matrix data
oncoprint_mat <- oncoprint_mat[-c(1, 2), ]

# modify accordingly "1" and "0" with more representative labels for mutations, copy number
# gains and losses
# for mutations:
drivers <- oncoprint_mat[1:4, ]
drivers[drivers == "3"] <- "methylation"
drivers[drivers == "2"] <- "mutation"
drivers[drivers == "1"] <- "mutation"
drivers[drivers == "0"] <- ""
# for copy number losses:
loss_cnv <- oncoprint_mat[5:7, ]
loss_cnv[loss_cnv == "1"] <- "copy_number_loss"
loss_cnv[loss_cnv == "0"] <- ""

# merge all the relabelled data, which will be what we will plot in the main body of the oncoprint
alterations <- rbind(drivers, loss_cnv)
# Oncoprint ---------------------------------------------------------------
# define the colors for each of the alterations
col <- c(
  "mutation" = "sandybrown",
  "methylation" = "lightgreen",
  "copy_number_loss" = "skyblue3"
)

# define function for transforming values of the matrix
# into the actual oncoprint values
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h,
      gp = gpar(fill = "grey80", col = "white")
    )
  },
  mutation = function(x, y, w, h) {
    grid.rect(x, y, w, h,
      gp = gpar(fill = col["mutation"], col = "white")
    )
  },
  methylation = function(x, y, w, h) {
    grid.rect(x, y, w, h,
      gp = gpar(fill = col["methylation"], col = "white")
    )
  },
  copy_number_loss = function(x, y, w, h) {
    grid.rect(x, y, w, h,
      gp = gpar(fill = col["copy_number_loss"], col = "white")
    )
  }
)

# define the colours for each of the annotation features that will be included
purity_col_fun <- colorRamp2(c(0, 1), c("white", "green4"))
ploidy_col_fun <- colorRamp2(c(1, 6), c("white", "salmon4"))
wgii_col <- colorRamp2(c(0, 1), c("white", "indianred4"))
ITH_col <- colorRamp2(c(0, 15), c("white", "navyblue"))
SampleType_col <- c(
  "Metastasis" = "tomato3", "Primary" = "green",
  "Lymph node" = "mediumorchid2", "Thrombus" = "cyan2",
  "Renal met" = "gold3", "#N/A" = "grey60"
)

stage_col <- c(
  "IV" = "#ff0000", "I" = "#ffff66",
  "III" = "#ff9900", "II" = "#ffcc00",
  "NA" = "grey60"
)

# create annotation for the oncoprint
col_ha <- HeatmapAnnotation(
  `Tumour purity` = info_oncoprint$purity,
  `Tumour ploidy` = info_oncoprint$ploidy,
  `wGII` = info_oncoprint$wgii,
  `stage` = info_oncoprint$stage,
  `Genetic ITH` = info_oncoprint$ITH,
  `Type of sample` = info_oncoprint$type_collapsed,
  annotation_legend_param = list(
    `Tumour purity` = list(direction = "horizontal"),
    `Tumour ploidy` = list(direction = "horizontal"),
    `wGII` = list(direction = "horizontal"),
    `Genetic ITH` = list(direction = "horizontal"),
    `Type of sample` = list(direction = "horizontal")
  ),
  gp = gpar(col = "white"),
  col = list(
    `Tumour purity` = purity_col_fun, `Tumour ploidy` = ploidy_col_fun,
    `wGII` = wgii_col, `Genetic ITH` = ITH_col, `Type of sample` = SampleType_col,
    `stage` = stage_col
  )
)


# Plot oncoprint and save it ----------------------------------------------

# pdf(file.path(OUT_DIR, "SupFig1_tx_rna_oncoprint.pdf"), width = 16, height = 4.5)
ht <- oncoPrint(alterations,
  alter_fun = alter_fun, col = col,
  top_annotation = NULL,
  right_annotation = NULL,
  row_order = 1:nrow(alterations),
  column_order = 1:ncol(alterations),
  bottom_annotation = col_ha,
  column_split = annotation$Patient
)
ht <- draw(ht, annotation_legend_side = "bottom")
# dev.off()
save_baseplot(ht,
  file.path(OUT_DIR, "SupFig1_tx_rna_oncoprint"),
  w = 16 * 25.4, h = 4.5 * 25.4
)

# ESTIMATE NUMBER PATIENTS WITH SUBCLONAL MUTATION ------------------------
n_pats <- length(unique(annotation$Patient))
get_n_subclonal(annotation, "SETD2")
get_n_subclonal(annotation, "PBRM1")
get_n_subclonal(annotation, "loss_9p")
get_n_subclonal(annotation, "loss_14q")
