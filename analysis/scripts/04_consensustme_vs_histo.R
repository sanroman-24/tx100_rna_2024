### Correlation between estimated cell type abundance consensusTME
### & histo estimation

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(plyr)
library(here)

# PATHS -------------------------------------------------------------------
BASE <- here::here()
HISTO_DIR <- file.path(BASE, "data", "meta", "histo")
HISTO_PATHS <- grep(".tsv", list.files(list.dirs(HISTO_DIR), full.names = TRUE), value = TRUE)
OUT_DIR <- file.path(BASE, "analysis", "figures")
TME_PATH <- file.path(BASE, "data", "processed", "tum_consensustme.rds")


# FUNCTIONS ---------------------------------------------------------------
# Returns merged dataframe with value of IHC marker and
# estimated abundance of specific cell type in TME from RNA-Seq
#' @param histo_df dataframe with expression of different markers
#' @param tme matrix with estimated abundance of cell counts (cells x samples)
#' @param ihc_marker label of IHC marker that we want to analyse
#' @param tme_cell_type label of cell population we want to associate with IHC marker
get_marker_df <- function(histo_df, tme, ihc_marker, tme_cell_types) {
  df <- tme %>%
    t() %>%
    as.data.frame() %>% # matrix to df
    rownames_to_column("sample") %>%
    dplyr::select("sample", tme_cell_types)
  df <- data.frame(sample = df$sample, cell_abundance = rowMeans(df[, -1]))
  histo_df <- histo_df %>%
    filter(IHC_marker == ihc_marker) %>%
    dplyr::select("sample", "Num Positive per mm^2")
  colnames(histo_df) <- c("sample", "ihc_estimate")
  return(merge(df, histo_df))
}

source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

# LOAD DATA ---------------------------------------------------------------
tme <- readRDS(TME_PATH)
histo <- do.call("rbind.fill", lapply(HISTO_PATHS, function(x) read_delim(x, delim = "\t")))

# reformat histo data
histo <- histo %>%
  # Patient ID is before first underscore always
  mutate(Patient = str_extract(Image, "^[^_]+")) %>%
  # Sample ID is before second underscore always
  mutate(sample = str_extract(Image, "^[^_]+_[^_]+")) %>%
  # Remove P that is not in my sample names
  mutate(sample = str_replace(sample, "_P", "_")) %>%
  # Marker analysed is after the last underscore always
  mutate(IHC_marker = str_extract(Image, "[^_]+$")) %>%
  # only consider the total tissue area column for correlations
  filter(Class == "Tissue_area")

# wrangle to make sample IDs match in both datasets
# remove trailing 0s in TME dataset
histo <- mutate(histo, sample = str_replace(sample, "R0+", "R"))
# change the - for _ in colnames in TRACERx Renal
tme <- tme[, colnames(tme) != "K207-R4"]
colnames(tme) <- str_replace_all(colnames(tme), "-", "_")

keep_smps <- dplyr::intersect(colnames(tme), histo$sample) # 53 total samples
# length(keep_smps) # 53 total samples
tme <- tme[, colnames(tme) %in% keep_smps]
histo <- histo[histo$sample %in% keep_smps, ]



# ASSOCIATION T CELL ABUNDANCE AND IHC CD3 QUANTIFICATION -----------------
cell_types <- c("T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "T_regulatory_cells")
marker <- "CD3"
df <- get_marker_df(histo, tme, marker, cell_types)

p <- scatter_plot(
  df = df[!is.na(df$ihc_estimate), ],
  y_str = "cell_abundance", x_str = "ihc_estimate",
  x_title = "T cell quantification IHC",
  y_title = "T cell estimation consensusTME"
) + geom_smooth(method = "lm", se = FALSE, col = "red") +
  ggpubr::stat_cor(method = "spearman", size = 2) + xlim(c(0, 4000))

save_ggplot(p = p, fp = file.path(OUT_DIR, "SupFigXX_Tcell_benchmark"), w = 50, h = 50)



# ASSOCIATION MACROPHAGE ABUNDANCE AND IHC CD68 QUANTIFICATION ------------
cell_types <- c("Macrophages_M1", "Macrophages_M2")
marker <- "PGM1"
df <- get_marker_df(histo, tme, marker, cell_types)

p <- scatter_plot(
  df = df,
  y_str = "cell_abundance", x_str = "ihc_estimate",
  x_title = "Macrophage quantification IHC",
  y_title = "Macrophage estimation consensusTME"
) + geom_smooth(method = "lm", se = FALSE, col = "red") +
  ggpubr::stat_cor(method = "spearman", size = 2) +
  xlim(c(0, 4500))

save_ggplot(p = p, fp = file.path(OUT_DIR, "SupFigXX_Macrophage_benchmark"), w = 50, h = 50)
