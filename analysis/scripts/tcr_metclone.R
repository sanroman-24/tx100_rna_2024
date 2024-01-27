# Compare TCR similarity primary-met in met clone vs not

rm(list = ls(all = TRUE))


# LIBRARIES ---------------------------------------------------------------

library(here)
library(tidyverse)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
OUT_DIR <- file.path(BASE, "analysis", "outputs")
PLOT_DIR <- file.path(BASE, "analysis", "figures")
TCR_SIM_PATH <- file.path(OUT_DIR, "tcr_repertoire_similarity.RDS")
ANNOTATION_PATH <- file.path(BASE, "data", "meta", "tx_annotation.tsv")
CLONES_PATH <- file.path(META_DIR, "clones_annotation.txt")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "get_clonal_dist.R"))
source(file.path(BASE, "src", "utils.R"))

# LOAD DATA ---------------------------------------------------------------
tcr_sim <- readRDS(TCR_SIM_PATH)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
clones_annotation <- read_delim(CLONES_PATH, delim = "\t")
clones_annotation <- clones_annotation[!is.na(clones_annotation$parent_clone), ]


# Calculate TCR similarity between mets and primaries
prim_lbl <- "PRIMARY"
met_lbl <- c("RENAL_MET", "THROMBUS", "LYMPH_NODE", "METASTASIS")
# binary variable that indicates if the sample is a metastasis
annotation$is_met <- ifelse(annotation$type_collapsed %in% prim_lbl, 1, 0)
# binary variable that indicates if the sample is from primary tumour
annotation$is_primary <- ifelse(annotation$type_collapsed %in% met_lbl, 1, 0)

# filter to patients with at least one met and primary
keep_pts <- annotation %>%
    group_by(Patient) %>%
    summarise(num_mets = sum(is_met), num_prims = sum(is_primary)) %>%
    filter(num_mets >= 1 & num_prims > 1) %>%
    dplyr::select(Patient) %>%
    as_vector()

annotation_sub <- annotation[annotation$Patient %in% keep_pts, ]

# dataframe where to store the final result
tcr_sim_sub <- data.frame(Patient = c(), s1 = c(), s2 = c(), sim = c())
for (patientID in unique(annotation_sub$Patient)) {
    samples <- annotation_sub$sample[annotation_sub$Patient == patientID]
    samples <- str_replace_all(samples, "-", "_")
    samples <- str_replace_all(samples, "_T(\\d+$)", "_THR\\1")
    pairs <- combn(samples, 2)
    for (pair in 1:ncol(pairs)) {
        s1 <- pairs[1, pair]
        s2 <- pairs[2, pair]
        sim <- tcr_sim[rownames(tcr_sim) == s1, colnames(tcr_sim) == s2]
        tcr_sim_sub <- rbind(
            tcr_sim_sub,
            data.frame(patient = patientID, s1 = s1, s2 = s2, sim = sim)
        )
    }
}


# subset only to primary-met pairs
mets <- annotation %>%
    filter(type_collapsed %in% met_lbl) %>%
    dplyr::select(sample, Patient, type_collapsed)

# make sure that always one of the sample is a met and the other is primary
tcr_sim_sub <- tcr_sim_sub[tcr_sim_sub$s1 %in% mets$sample |
    tcr_sim_sub$s2 %in% mets$sample, ]

tcr_sim_sub <- tcr_sim_sub[!tcr_sim_sub$s1 %in% mets$sample |
    !tcr_sim_sub$s2 %in% mets$sample, ]

tcr_sim_sub <- tcr_sim_sub %>%
    mutate(
        prim_sample = ifelse(s1 %in% mets$sample, s2, s1),
        met_sample = ifelse(s1 %in% mets$sample, s1, s2)
    ) %>%
    dplyr::select(met_sample, prim_sample, sim)
tcr_sim_sub <- merge(tcr_sim_sub, annotation, by.x = "met_sample", by.y = "sample")


# CALCULATE NUMBER OF CLONES SEPARATING SAMPLES ---------------------------
clonal_dists <- sapply(1:nrow(tcr_sim_sub), function(i) {
    prim_smp <- str_replace(tcr_sim_sub$prim_sample[i], "-", "_")
    met_smp <- str_replace(tcr_sim_sub$met_sample[i], "-", "_")
    get_met_prim_dist(met_smp, prim_smp, clones_annotation)
})


# CORRELATION CLONAL DISTANCE AND TRANSCRIPTIONAL DIVERGENCE --------------
tcr_sim_sub$clonal_dists <- clonal_dists
tcr_sim_sub$metastasising_clone <- ifelse(clonal_dists == 0,
    "yes_met_clone", "no_met_clone"
)

# summary plot across all the patients

tcr_sim_sub <- tcr_sim_sub[!is.na(tcr_sim_sub$clonal_dists), ]

p <- violin_plot(
    df = tcr_sim_sub, x_str = "metastasising_clone",
    y_str = "sim", x_title = "Primary region",
    y_title = "Transcriptional distance",
    labs_x = c("No metastasising clone", "Metastasising clone"),
    fun.y = "mean",
    ylim = c(0, 1)
)

p <- p + scale_x_discrete(limits = c("yes_met_clone", "no_met_clone"))

save_ggplot(p, file.path(FIG_DIR, "Fig6E_tcr_sim_metclone"), w = 40, h = 40)

summary(run_lme("sim", "metastasising_clone", "Patient", tcr_sim_sub))
# p_value = 0.047 for greater similarity to primary regions with met clone