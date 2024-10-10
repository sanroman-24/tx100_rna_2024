# Spatial changes in ENPP1 and SLC19A1 correlation with immune infiltrate
rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(survival)
library(survminer)
library(here)

# PATHS
BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
TME_PATH <- file.path(BASE, "data", "processed", "tum_consensustme.rds")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
SSGSEA_TUMOURS_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")
VST_PATH <- file.path(BASE, "data", "processed", "tum_vst_exp_filtered.rds")

# LOAD DATA
tme <- read_rds(TME_PATH)
cells <- rownames(tme)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
vst <- readRDS(VST_PATH)
ssgsea_tumour <- read_delim(SSGSEA_TUMOURS_PATH)
ssgsea_tumour$patient <- str_remove(ssgsea_tumour$sample, "[-_].*$")


# FUNCTIONS
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "plotting_theme.R"))

# WRANGLE DATA ---------------------------------------------------------------
# filter to genes of interest
vst <- assay(vst[rownames(vst) %in% c("ENPP1", "SLC19A1"), ]) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
rownames(vst) <- str_sub(rownames(vst), 1, 12)

vst <- rownames_to_column(vst, "sample") %>%
    mutate(Patient = str_sub(sample, 1, 4))

# add information and classification of interest to sample annotation
annotation <- merge(vst, annotation)
annotation <- annotation %>%
    mutate(
        high_SLCENPP =
            case_when(
                ENPP1 > quantile(ENPP1, .75, na.rm = T) |
                    SLC19A1 > quantile(SLC19A1, .75, na.rm = TRUE) ~ 1,
                TRUE ~ 0
            )
    ) %>%
    mutate(high_wgii = case_when(
        wgii > quantile(wgii, .5, na.rm = T) ~ 1,
        TRUE ~ 0
    )) %>%
    mutate(SLCENPP_wgii = case_when(
        high_SLCENPP & high_wgii ~ "rep_wgii",
        !high_SLCENPP & high_wgii ~ "norep_wgii",
        high_SLCENPP & !high_wgii ~ "rep_nowgii",
        !high_SLCENPP & !high_wgii ~ "norep_nowgii"
    ))

summary(nlme::lme(high_SLCENPP ~ high_wgii, random = ~1|Patient, data = annotation))
summary(nlme::lme(high_SLCENPP ~ loss_9p, random = ~1|Patient, data = annotation[annotation$high_wgii == 1,]))
summary(nlme::lme(high_SLCENPP ~ loss_9p, random = ~1|Patient, data = annotation[annotation$high_wgii == 0,]))

tme <- tme[, colnames(tme) != "K207-R4"]
tme <- t(tme) %>%
    as.data.frame() %>%
    rownames_to_column("sample")
df <- merge(tme, annotation, by = "sample")

df <- merge(df, ssgsea_tumour)

# Check association with TME across samples

df$t2myel <- scale(df$motzer_myeloid_inflammation) - 
             scale(df$motzer_t_effector)


summary(run_lme("motzer_t_effector", "high_SLCENPP", "Patient", df))
summary(run_lme("motzer_myeloid_inflammation", "high_SLCENPP", "Patient", df))
summary(run_lme("t2myel", "high_SLCENPP", "Patient", df))


cgas_plt <- c("wgii_high_slcenpp_high" = "tomato3", 
              "wgii_low_slcenpp_high" = "grey25", 
              "wgii_high_slcenpp_low" = "tomato1", 
              "wgii_low_slcenpp_low" = "grey50")


df <- mutate(df, group = case_when(
        df$high_SLCENPP & df$high_wgii ~ "wgii_high_slcenpp_high",
        !df$high_SLCENPP & df$high_wgii ~ "wgii_high_slcenpp_low",
        df$high_SLCENPP & !df$high_wgii ~ "wgii_low_slcenpp_high",
        !df$high_SLCENPP & !df$high_wgii ~ "wgii_low_slcenpp_low"
    )
)

p <- ggplot(df, aes(
    y = motzer_t_effector, x = group
)) +
    geom_boxplot(aes(fill = group)) +
    theme(legend.position = "none") + 
    scale_fill_manual(values = cgas_plt) + 
    scale_x_discrete(limits = c(
        "wgii_low_slcenpp_low", 
        "wgii_high_slcenpp_low",
         "wgii_low_slcenpp_high", 
         "wgii_high_slcenpp_high"
    )) + 
    ggpubr::stat_compare_means(
        comparisons = list(
            c("wgii_high_slcenpp_high", "wgii_high_slcenpp_low"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_high"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_low")
        )
    )

save_ggplot(p, file.path(FIG_DIR, "slcenpp_teff"), w = 50, h = 50)


p <- ggplot(df, aes(
    y = motzer_myeloid_inflammation, x = group
)) +
    geom_boxplot(aes(fill = group)) +
    theme(legend.position = "none") + 
    scale_fill_manual(values = cgas_plt) + 
    scale_x_discrete(limits = c(
        "wgii_low_slcenpp_low", 
        "wgii_high_slcenpp_low",
         "wgii_low_slcenpp_high", 
         "wgii_high_slcenpp_high"
    )) + ggpubr::stat_compare_means(
        comparisons = list(
            c("wgii_high_slcenpp_high", "wgii_high_slcenpp_low"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_high"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_low")
        )
    )


save_ggplot(p, file.path(FIG_DIR, "slcenpp_myel"), w = 50, h = 50)



p <- ggplot(df, aes(
    y = t2myel, x = group
)) +
    geom_boxplot(aes(fill = group)) +
    theme(legend.position = "none") + 
    scale_fill_manual(values = cgas_plt) + 
    scale_x_discrete(limits = c(
        "wgii_low_slcenpp_low", 
        "wgii_high_slcenpp_low",
         "wgii_low_slcenpp_high", 
         "wgii_high_slcenpp_high"
    )) + 
    ggpubr::stat_compare_means(
        comparisons = list(
            c("wgii_high_slcenpp_high", "wgii_high_slcenpp_low"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_high"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_low")
        )
    )


save_ggplot(p, file.path(FIG_DIR, "slcenpp_t2myel"), w = 50, h = 50)

p <- ggplot(df, aes(
    y = motzer_t_effector, x = group
)) +
    geom_boxplot(aes(fill = group)) +
    theme(legend.position = "none") + 
    scale_fill_manual(values = cgas_plt) + 
    scale_x_discrete(limits = c(
        "wgii_low_slcenpp_low", 
        "wgii_high_slcenpp_low",
         "wgii_low_slcenpp_high", 
         "wgii_high_slcenpp_high"
    )) + 
    ggpubr::stat_compare_means(
        comparisons = list(
            c("wgii_high_slcenpp_high", "wgii_high_slcenpp_low"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_high"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_low")
        )
    )


save_ggplot(p, file.path(FIG_DIR, "slcenpp_wgii_teff"), w = 50, h = 50)

p <- ggplot(df, aes(
    y = motzer_myeloid_inflammation, x = group
)) +
    geom_boxplot(aes(fill = group)) +
    theme(legend.position = "none") + 
    scale_fill_manual(values = cgas_plt) + 
    scale_x_discrete(limits = c(
        "wgii_low_slcenpp_low", 
        "wgii_high_slcenpp_low",
         "wgii_low_slcenpp_high", 
         "wgii_high_slcenpp_high"
    )) + 
    ggpubr::stat_compare_means(
        comparisons = list(
            c("wgii_high_slcenpp_high", "wgii_high_slcenpp_low"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_high"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_low")
        )
    )


save_ggplot(p, file.path(FIG_DIR, "slcenpp_wgii_myel"), w = 50, h = 50)


p <- ggplot(df, aes(
    y = t2myel, x = group
)) +
    geom_boxplot(aes(fill = group)) +
    theme(legend.position = "none") + 
    scale_fill_manual(values = cgas_plt) + 
    scale_x_discrete(limits = c(
        "wgii_low_slcenpp_low", 
        "wgii_high_slcenpp_low",
         "wgii_low_slcenpp_high", 
         "wgii_high_slcenpp_high"
    )) + 
    ggpubr::stat_compare_means(
        comparisons = list(
            c("wgii_high_slcenpp_high", "wgii_high_slcenpp_low"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_high"),
            c("wgii_high_slcenpp_high", "wgii_low_slcenpp_low")
        )
    )


save_ggplot(p, file.path(FIG_DIR, "slcenpp_wgii_t2myel"), w = 50, h = 50)

# Subset to paired analysis
tb <- table(df$high_SLCENPP, df$Patient)
pts <- colnames(tb)[tb[1, ] & tb[2, ]]

df_sub <- df[df$Patient %in% pts, ]
summary(run_lme("motzer_t_effector", "high_SLCENPP", "Patient", df_sub))
summary(run_lme("t2myel", "high_SLCENPP", "Patient", df_sub))
summary(run_lme("Macrophages", "high_SLCENPP", "Patient", df_sub))

p <- df_sub %>%
    mutate(high_SLCENPP = ifelse(
        high_SLCENPP == 1, "cGAS repressed", "WT"
    )) %>%
    df_to_ggpaired(
        pat_var = "Patient",
        cond = "high_SLCENPP", var_str = "t2myel"
    ) %>%
    plot_paired_boxplot(
        cond1 = "cGAS repressed", cond2 = "WT",
        ylab = "Myeloid Inflammation - T effector",
        xlab = "", ylim = c(NA)
    )

save_ggplot(p,
    file.path(FIG_DIR, paste0("paired_slcenpp_t2myel")),
    w = 45, h = 45
)
