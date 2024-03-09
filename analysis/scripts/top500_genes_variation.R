# Calculate top 500 genes with highest variation in TRACERx Renal
# (Used to calculate I-TED as in Martinez-Ruiz et al, Nature 2023)

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(here)


# PATHS -------------------------------------------------------------------
BASE <- here::here()
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")

# LOAD DATA ---------------------------------------------------------------
vst <- readRDS(VST_PATH)

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

# FIND TOP 500 WITH HIGHEST VARIATION -------------------------------------

# Filter to genes with at least 5 counts in at least 20% of the cohort
# The filter here is more strict to avoid high variability in genes with
# problematic detection only in a subset of samples
vst <- assay(vst)
keep_genes <- (rowMeans(vst > 5) > .2)
vst <- vst[keep_genes, ] # 16704 genes passed filtering
top500_genes <- apply(vst, 1, sd) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500) %>%
    {
        names(.)
    }
saveRDS(top500_genes, file.path(OUT_DIR, "top500_variable_genes.rds"))

## Plot variance per gene
v <- apply(vst, 1, sd) %>% sort(decreasing = TRUE)
df <- data.frame(variance = v, genes = names(v))
df$hghl <- ifelse(df$genes %in% top500_genes, "yes", "no")

p <- ggplot(
    df, aes(x = reorder(genes, v), y = v, fill = hghl)
    ) +
    geom_col() +
    theme(legend.position = "none") +
    scale_fill_manual(
        values = c("no" = "grey90", "yes" = "coral")
    ) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

save_ggplot(p,file.path(FIG_DIR, "tr_var_pergene"), w = 100, h = 50)
