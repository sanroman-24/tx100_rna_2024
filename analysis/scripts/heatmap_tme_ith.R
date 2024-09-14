# TME ITH groups in TRACERx Renal

rm(list = ls(all = TRUE))

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(circlize)
library(ComplexHeatmap)


# PATHS -------------------------------------------------------------------

BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
PLOT_DIR <- file.path(BASE, "analysis", "figures")
META_DIR <- file.path(BASE, "data", "meta")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
TME_PATH <- file.path(BASE, "data", "processed", "tum_consensustme.rds")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))


# LOAD DATA ---------------------------------------------------------------
tme <- read_rds(TME_PATH)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
annotation$ITH <- as.numeric(annotation$ITH)
annotation$wgii <- as.numeric(annotation$wgii)
annotation$Regions <- as.numeric(annotation$Regions)
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH)
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]
clinical_data <- clinical_data %>%
    dplyr::rename("Patient" = "Subject", "Stage" = "Overall Stage")
annotation <- left_join(annotation, clinical_data, by = "Patient")

annotation <- annotation[annotation$sample != "K207-R4", ]
tme <- tme[, colnames(tme) != "K207-R4"]

# check that it is good to go
all(annotation$sample == colnames(tme))

# Subset to patients with more than one tumour region (regardless of loc)
keep_pts <- unique(annotation$Patient[annotation$Regions > 1])
annotation <- annotation[annotation$Patient %in% keep_pts, ]
tme <- t(tme)
tme <- tme[rownames(tme) %in% annotation$sample, ]
# check that same values in annotation and consensustme data
all(rownames(tme) == annotation$sample)


# TRANSFORM CELL ABUNDANCE TO Z-SCORE -----------------------------------
tme <- tme %>%
    as.data.frame() %>%
    dplyr::select(-c(Fibroblasts, Endothelial, Immune_Score)) %>%
    as.matrix() %>%
    scale()

colnames(tme) <- str_replace_all(colnames(tme), "_"," ")
# CONSTRUCT HEATMAP -----------------------------------------------------
ht <- Heatmap(t(tme),
    name = "Row Z-Scores",
    col = circlize::colorRamp2(c(-2, 0, 2), c(tx_palette[["darkblue"]], "white", tx_palette[["darkred"]])),
    column_dend_reorder = FALSE,
    clustering_distance_columns = "manhattan",
    clustering_method_columns = "complete", clustering_method_rows = "complete",
    column_split = 3,
    show_column_dend = TRUE, show_row_dend = FALSE,
    show_column_names = TRUE, show_heatmap_legend = TRUE,
    column_names_gp = gpar(fontsize = 6),
    row_names_gp = gpar(fontsize = 10), rect_gp = gpar(col = "white", lwd = 0.5)
)

ht <- draw(ht)


# GET ANNOTATION -----------------------------------------------------

# this returns a list split with numbers that would reorder the matrix
# columns as we see it in the heatmap, with values split into 3 elements
# of the list corresponding to the 3 clusters: cold, intermediate, hot
order_columns <- column_order(ht)

# define vector that labels each column of the matrix(sample) accordingly
category_cluster <- vector(length = nrow(tme))
category_cluster[order_columns[[1]]] <- "hot"
category_cluster[order_columns[[2]]] <- "cold"
category_cluster[order_columns[[3]]] <- "intermediate"

# get all the sample IDs
samples <- rownames(tme)
# get column order in the heatmap to know in which cluster each patient
order_columns <- column_order(ht)
# get the number of samples in each cluster
n_samples_cluster1 <- length(order_columns[[1]])
n_samples_cluster2 <- length(order_columns[[2]])
n_samples_cluster3 <- length(order_columns[[3]])
# get samples of each cluster in the order in which they appear in the heatmap
samples_cluster1 <- samples[order_columns[[1]]]
samples_cluster2 <- samples[order_columns[[2]]]
samples_cluster3 <- samples[order_columns[[3]]]
# combine all the samples to get all IDs in the order in the clustering
all_samples <- c(samples_cluster1, samples_cluster2, samples_cluster3)
# get the patient IDs. These are the first 4 digits always
all_patients <- substr(all_samples, 1, 4)
# placeholder matrix
matrix_annotation <- matrix(
    nrow = length(unique(all_patients)),
    ncol = length(all_patients)
)

# set values to "0" initially to faciliate posterior modification
matrix_annotation[is.na(matrix_annotation)] <- "0"

# add column and rownames
colnames(matrix_annotation) <- all_patients
rownames(matrix_annotation) <- unique(all_patients)

# when there is a match in patient between columns and rows
# (which means different samples from the same patient by
# how we define the matrix), include the proper label

for (i in 1:nrow(matrix_annotation)) {
    for (j in 1:ncol(matrix_annotation)) {
        if (colnames(matrix_annotation)[j] == rownames(matrix_annotation)[i]) {
            if (j <= n_samples_cluster1) {
                matrix_annotation[i, j] <- "hot"
            } else if ((j > n_samples_cluster1) & (j <= (n_samples_cluster1 + n_samples_cluster2))) {
                matrix_annotation[i, j] <- "cold"
            } else if (j > (n_samples_cluster1 + n_samples_cluster2)) {
                matrix_annotation[i, j] <- "intermediate"
            }
        }
    }
}

# now, include intermediate values. If the samples from different regions, mark it.
for (i in 1:nrow(matrix_annotation)) {
    # get index of the cells from the same patient, which were labelled before accordingly
    indixes <- which(!str_detect(matrix_annotation[i, ], "0"))
    interval_to_fill <- range(indixes)
    fill_from <- interval_to_fill[1]
    fill_to <- interval_to_fill[2]
    # decipher where to colour with the lines of all immune_low, all immune_hot, all immune_inter
    # or immune_heterogeneous
    if (all(matrix_annotation[i, ][indixes] == "cold")) {
        connecting_values <- "homogeneous_cold"
    } else if (all(matrix_annotation[i, ][indixes] == "hot")) {
        connecting_values <- "homogeneous_hot"
    } else if (all(matrix_annotation[i, ][indixes] == "intermediate")) {
        connecting_values <- "homogeneous_intermediate"
    } else {
        connecting_values <- "heterogeneous"
    }
    # modify cells without matches without altering already written values
    matrix_annotation[i, fill_from:fill_to][-c(indixes - fill_from + 1)] <- connecting_values
}

# change 0 values to "", so that oncoprint works perfectly with them
matrix_annotation <- str_replace(matrix_annotation, "0", "") %>%
    matrix(nrow = length(unique(all_patients)), ncol = length(all_patients))

# add column and rownames again to the recreated matrix of the prior step
colnames(matrix_annotation) <- all_samples
rownames(matrix_annotation) <- unique(all_patients)

# define colours to use in oncoprint
col <- c(
    "hot" = tx_palette[["darkred"]], "cold" = tx_palette[["darkblue"]], 
    "intermediate" = tx_palette[["lightpurple"]],
    "heterogeneous" = tx_palette[["darkpurple"]], "homogeneous_cold" = tx_palette[["darkblue"]],
    "homogeneous_hot" = tx_palette[["darkred"]], "homogeneous_intermediate" = tx_palette[["lightpurple"]]
)

# define alter function for oncoprint. Cells full size, except for connection
# which will be with reduced height so that they look like connecting lines
alter_fun <- list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w, h,
            gp = gpar(fill = "grey99", col = NA)
        )
    },
    # hot
    hot = function(x, y, w, h) {
        grid.rect(x, y, w, h,
            gp = gpar(fill = col["hot"], col = NA)
        )
    },
    # intermediate
    intermediate = function(x, y, w, h) {
        grid.rect(x, y, w, h,
            gp = gpar(fill = col["intermediate"], col = NA)
        )
    },
    # cold
    cold = function(x, y, w, h) {
        grid.rect(x, y, w, h,
            gp = gpar(fill = col["cold"], col = NA)
        )
    },
    # heterogeneous
    heterogeneous = function(x, y, w, h) {
        grid.rect(x, y, w, h * 0.1,
            gp = gpar(fill = col["heterogeneous"], col = NA)
        )
    },
    # homogeneous_cold
    homogeneous_cold = function(x, y, w, h) {
        grid.rect(x, y, w, h * 0.1,
            gp = gpar(fill = col["homogeneous_cold"], col = NA)
        )
    },
    # homogeneous_intermediate
    homogeneous_intermediate = function(x, y, w, h) {
        grid.rect(x, y, w, h * 0.1,
            gp = gpar(fill = col["homogeneous_intermediate"], col = NA)
        )
    },
    # homogeneous_hot
    homogeneous_hot = function(x, y, w, h) {
        grid.rect(x, y, w, h * 0.1,
            gp = gpar(fill = col["homogeneous_hot"], col = NA)
        )
    }
)

# uncomment to visualize oncoPrint
#
# oncoPrint(matrix_annotation, alter_fun = alter_fun,
# column_order = 1:ncol(matrix_annotation),
# row_order = 1:nrow(matrix_annotation),
# show_pct = FALSE, name = "annotation")


# reorder matrix of all cases
matrix_annotation2 <- matrix(
    nrow = length(unique(all_patients)),
    ncol = nrow(tme)
)
colnames(matrix_annotation2) <- rownames(tme)
rownames(matrix_annotation2) <- unique(all_patients)
for (name_col in colnames(matrix_annotation2)) {
    tmp <- matrix_annotation[, colnames(matrix_annotation) == name_col]
    matrix_annotation2[, colnames(matrix_annotation2) == name_col] <- tmp
}

stage <- annotation[unlist(order_columns), ]$Stage
type_sample <- annotation[unlist(order_columns), ]$type_collapsed

# set colours to use
pal <- RColorBrewer::brewer.pal(n = 4, name = "RdBu")
stage_col <- c("I" = pal[4], "II" = pal[3], "III" = pal[2], "IV" = pal[1])
type_sample_col <- c(
    "METASTASIS" = "tomato3", "PRIMARY" = "lightgreen",
    "LYMPH_NODE" = "mediumorchid2", "THROMBUS" = "cyan2",
    "RENAL_MET" = "lightcoral", "#N/A" = "grey50"
)
immune_col <- c(
    "cold" = "darkblue",
    "intermediate" = "salmon1",
    "hot" = "darkred"
)

ha <- HeatmapAnnotation(
    stage = stage,
    type_sample = type_sample,
    immune_cluster = category_cluster,
    col = list(
        stage = stage_col,
        type_sample = type_sample_col,
        immune_cluster = immune_col
    ),
    gp = gpar(col = "white"),
    height = unit(5, "mm"),
    width = unit(50, "mm")
)


# CONCAT HEATMAP AND ANNOTATION -----------------------------------------------
ht <- Heatmap(t(tme),
    name = "Row Z-Scores",
    col = circlize::colorRamp2(c(-2, 0, 2), c("darkblue", "white", "darkred")),
    column_dend_reorder = FALSE,
    clustering_distance_columns = "manhattan",
    clustering_method_columns = "complete", clustering_method_rows = "complete",
    show_column_dend = TRUE, show_row_dend = FALSE,
    show_column_names = FALSE, show_heatmap_legend = TRUE,
    # column_names_gp =  gpar(fontsize = 6),
    row_names_gp = gpar(fontsize = 4),
    height = unit(25, "mm"),
    width = unit(50, "mm")
)

# legend for parameters

oncoprint_legend_param <- list(
    title = "",
    at = c(
        "hot", "intermediate", "cold", "heterogeneous",
        "homogeneous_cold", "homogeneous_intermediate",
        "homogeneous_hot"
    ),
    labels = c(
        "Highly immune infiltrated region",
        "Intermediately immune infiltrated region",
        "Lowly immune infiltrated region",
        "Regions from patient with varied immune infiltration",
        "All regions of patient with low immune infiltration",
        "All regions of patient with intermediate immune infiltration",
        "All regions of patient with high immune infiltration"
    )
)

ht2 <- oncoPrint(matrix_annotation2,
    alter_fun = alter_fun,
    column_order = column_order(ht),
    row_order = 1:nrow(matrix_annotation), show_pct = FALSE,
    row_names_gp = gpar(fontsize = 0),
    name = "annotation", height = unit(30, "mm"), width = unit(50, "mm")
)

png(file.path(PLOT_DIR, "Fig4A_ITH_immune_heatmap.png"),
    width = 180, height = 100, units = "mm", res = 400
)
ht %v% ht2
dev.off()


pdf(file.path(PLOT_DIR, "Fig4A_ITH_immune_heatmap.pdf"),
    width = 180 / 25.4, height = 80 / 25.4
)
ht %v% ht2
dev.off()


svg(file.path(PLOT_DIR, "Fig4A_ITH_immune_heatmap.svg"),
    width = 180 / 25.4, height = 100 / 25.4
)
ht %v% ht2
dev.off()

# SURVIVAL BY SUBGROUP -----------------------------------------------

# lbls = sapply(1:nrow(matrix_annotation), function(i) {
#     r = matrix_annotation[i, ]
#     region_lbl = c("hot", "intermediate", "cold", "")
#     if (sum(!r %in% region_lbl) == 0){
#         lbl = paste0("homogeneous_", unique(r[!r %in% ""]))
#     } else {
#         lbl = unique(r[!r %in% c("hot", "cold", "intermediate", "")])
#     }
# })

# names(lbls) = rownames(matrix_annotation)

# df = data.frame(Patient = names(lbls), lbl = lbls)
# write_delim(df, file.path(OUT_DIR, "TME_het_group.tsv"), delim = "\t")

# clinical_data = merge(clinical_data, df, by = "Patient") %>%
#   mutate(pfs = ifelse(`PFS (months)` == "-", 0, 1)) %>%
#   mutate(pfs_time = as.numeric(ifelse(`PFS (months)` == "-",
#     `Total follow up (months)`, `PFS (months)`
#   ))) %>%
#   mutate(OS = ifelse(Outcome == "Death", 1, 0)) %>%
#   mutate(OS_time = as.numeric(`Total follow up (months)`))

# # only 1 homogeneous cold patient, so remove
# clinical_data = clinical_data[clinical_data$lbl != "homogeneous_cold", ]

# plt = "lancet"
# # PFS
# km_fit <- survfit(Surv(pfs_time, pfs) ~ lbl,
#   data = clinical_data
# )

# km_plot <- plot_km(clinical_data, km_fit, plt = plt)

# save_baseplot(km_plot, file.path(PLOT_DIR, "SuppFig12_PFS_TME_groups"), w = 70, h = 90)

# coxph(Surv(pfs_time, pfs) ~ lbl, data = clinical_data)

# # OS
# km_fit <- survfit(Surv(OS_time, OS) ~ lbl,
#   data = clinical_data
# )

# km_plot <- plot_km(clinical_data, km_fit, plt = plt)

# save_baseplot(km_plot, file.path(PLOT_DIR, "SuppFig12_OS_TME_groups"), w = 70, h = 90)

# coxph(Surv(OS_time, OS) ~ lbl, data = clinical_data)
