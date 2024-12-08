# Calculate prognostic association of proliferation + 9p loss in TRACERx Renal

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------


library(tidyverse)
library(survminer)
library(survival)
library(ggthemes)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
DATA_DIR <- file.path(BASE, "data", "raw")
META_DIR <- file.path(BASE, "data", "meta")
OUT_DIR <- file.path(BASE, "analysis", "outputs")
FIG_DIR <- file.path(BASE, "analysis", "figures")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
CLINICAL_ANNOTATION_PATH <- file.path(META_DIR, "TRACERx_s1_1_clinical.txt")
SSGSEA_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")
HALLMARK_GROUPS_PATH <- file.path(META_DIR, "martinez_ruiz_2023_hallmark_gs_groups.txt")

# FUNCTIONS ----------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

get_p_comp <- function(grp1, grp2, surv_end, df, gr_var) {
    df <- df[df[[gr_var]] %in% c(grp1, grp2), ]
    df$gr <- df[[gr_var]]
    if (surv_end == "OS") {
        return(summary(coxph(Surv(os_time, os) ~ gr, data = df)))
    } else if (surv_end == "PFS") {
        return(summary(coxph(Surv(pfs_time, pfs) ~ gr, data = df)))
    } else {
        stop("Please provide valid clinical endpoint")
    }
}

# LOAD DATA ---------------------------------------------------------------
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
gene_groups <- read_delim(HALLMARK_GROUPS_PATH)
ssgsea_scores <- read_delim(SSGSEA_PATH)
# clinical annotation of the TRACERx Renal cohort
clinical_data <- read_delim(CLINICAL_ANNOTATION_PATH, delim = "\t")
colnames(clinical_data) <- clinical_data[1, ]
clinical_data <- clinical_data[-1, ]
annotation <- annotation[annotation$sample != "K207-R4", ]

# only primary samples
annotation <- annotation[annotation$type_collapsed == "PRIMARY", ]
# get a "proliferation" score
annotation <- merge(annotation, ssgsea_scores, "sample")

gene_groups$pathway <- paste0("HALLMARK_", str_to_upper(gene_groups$Hallmark))
pf_paths <- gene_groups$pathway[gene_groups$Functional_group == "proliferation"]
ifn_paths <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")

scores <- annotation %>%
    dplyr::select(Patient, pf_paths, ifn_paths) %>%
    group_by(Patient) %>%
    summarise_all(max) %>%
    {
        cbind(
            .[, 1],
            data.frame(proliferation = rowMeans(.[, pf_paths])),
            ifn = rowMeans(.[, ifn_paths])
        )
    }

# Add information on status 9p loss
loss_9p <- annotation %>%
    group_by(Patient) %>%
    summarise(loss_9p = sum(loss_9p)) %>%
    mutate(loss_9p = ifelse(loss_9p > 0, 1, 0))

scores <- merge(scores, loss_9p)

scores$Subject <- scores$Patient
scores <- merge(scores, clinical_data, by = "Subject") %>%
    mutate(pfs = ifelse(`PFS (months)` == "-", 0, 1)) %>%
    mutate(pfs_time = as.numeric(ifelse(`PFS (months)` == "-",
        `Total follow up (months)`, `PFS (months)`
    ))) %>%
    mutate(os = ifelse(Outcome == "Death", 1, 0)) %>%
    mutate(os_time = as.numeric(`Total follow up (months)`))

prol_thresh <- quantile(scores$proliferation, .75)
scores$high_proliferation <- scores$proliferation > prol_thresh
scores <- scores %>% mutate(group = case_when(
    high_proliferation & loss_9p ~ "9p loss + high proliferation",
    !high_proliferation & loss_9p ~ "9p loss + low proliferation",
    high_proliferation & !loss_9p ~ "9p WT + high proliferation",
    !high_proliferation & !loss_9p ~ "9p WT + low proliferation",
))

write_delim(scores, file.path(OUT_DIR, "surv_9p_ssgsea_scores.tsv"), delim = "\t")

### SURVIVAL ANALYSIS
# PFS
plt <- "lancet"

km_fit <- survfit(Surv(pfs_time, pfs) ~ group,
    data = scores
)

km_plot <- plot_km(scores, km_fit, plt = plt)

save_baseplot(km_plot, file.path(FIG_DIR, "SupFig13_9p_loss_prolif_pfs"),
    h = 90, w = 70
)

# OS
km_fit <- survfit(Surv(os_time, os) ~ group,
    data = scores
)

km_plot <- plot_km(scores, km_fit, plt = plt)

save_baseplot(km_plot, file.path(FIG_DIR, "SupFig13_9p_loss_prolif_os"),
    h = 90, w = 70
)


get_p_comp(
    "9p loss + high proliferation", "9p loss + low proliferation",
    "OS", scores, "group"
)

get_p_comp(
    "9p loss + high proliferation", "9p loss + low proliferation",
    "PFS", scores, "group"
)

get_p_comp(
    "9p loss + low proliferation", "9p WT + low proliferation",
    "OS", scores, "group"
)

get_p_comp(
    "9p loss + low proliferation", "9p WT + low proliferation",
    "PFS", scores, "group"
)


# Check upon correction by stage
scores$stage_simple <- ifelse(scores$`Overall Stage` %in% c("I", "II"), "I-II", "III-IV")
coxph(Surv(pfs_time, pfs) ~ loss_9p:high_proliferation + stage_simple, data = scores)
coxph(Surv(os_time, pfs) ~ proliferation + stage_simple, data = scores[scores$loss_9p > 0, ])
coxph(Surv(os_time, pfs) ~ proliferation + stage_simple, data = scores[scores$loss_9p == 0, ])


# # 9p loss status
# hd_9p <- read_delim("tmp/9p_status.csv")
# hd_9p$Patient <- str_remove(hd_9p$Patient, "_.+")
# pts <- scores$Patient[scores$group == "9p loss + high proliferation"] 
# sapply(pts, function(pt){
#     if (any(hd_9p$status_9p[hd_9p$Patient == pt] == "HD")){
#         return("HD")
#     } else if (any(hd_9p$status_9p[hd_9p$Patient == pt] == "LOH")){
#         return("LOH")
#     } else {
#         return("WT")
#     }
# })

