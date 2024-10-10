# Analysis of the transitions between anti-tumor and immunosuppressive TME

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(here)

# PATHS
BASE <- here::here()
META_DIR <- file.path(BASE, "data", "meta")
FIG_DIR <- file.path(BASE, "analysis", "figures", "review")
TME_PATH <- file.path(BASE, "data", "processed", "tum_consensustme.rds")
ANNOTATION_PATH <- file.path(META_DIR, "tx_annotation.tsv")
SSGSEA_TUMOURS_PATH <- file.path(BASE, "data", "processed", "tx_ssGSEA.tsv")
# from ssgsea_early_late.R
CLONES_ANNOTATION_PATH <- file.path(
    BASE, "analysis", "outputs",
    "clones_annotation_with_distance_mrca.tsv"
)

# LOAD DATA
clones_annotation <- read_delim(CLONES_ANNOTATION_PATH, delim = "\t")
tme <- read_rds(TME_PATH)
cells <- rownames(tme)
annotation <- read_delim(ANNOTATION_PATH, delim = "\t")
ssgsea_tumour <- read_delim(SSGSEA_TUMOURS_PATH)
ssgsea_tumour$patient <- str_remove(ssgsea_tumour$sample, "[-_].*$")


# FUNCTIONS
source(file.path(BASE, "src", "utils.R"))
source(file.path(BASE, "src", "plotting_functions.R"))
source(file.path(BASE, "src", "plotting_theme.R"))

get_coocurrence_transition <- function(transitions_df, df, v, is_driver) {
    m_s1 <- match(transitions_df$s1, df$sample)
    m_s2 <- match(transitions_df$s2, df$sample)
    if (is_driver) {
        df[[v]] <- ifelse(df[[v]] > 0, 1, df[[v]])
    }
    transitions_df[[paste0(v, "_s1")]] <- df[[v]][m_s1]
    transitions_df[[paste0(v, "_s2")]] <- df[[v]][m_s2]
    transitions_df <- transitions_df[complete.cases(transitions_df), ]
    transitions_df[[paste0(v, "_change")]] <- paste0(
        transitions_df[[paste0(v, "_s1")]],
        "_",
        transitions_df[[paste0(v, "_s2")]]
    )
    return(table(
        transitions_df[[paste0(v, "_change")]],
        transitions_df$transition
    ))
}

# WRANGLE DATA ---------------------------------------------------------------
tme <- tme[, colnames(tme) != "K207-R4"] 
tme <- t(tme) %>%
    as.data.frame() %>%
    rownames_to_column("sample")
annotation <- annotation[annotation$type_collapsed == "PRIMARY", ]
df <- merge(tme, annotation, by = "sample")

df <- merge(df, ssgsea_tumour)

# Only focus in patients with more than one tumor sample
patients <- table(df$patient) %>%
    {
        names(.)[. > 1]
    }

df <- df[df$patient %in% patients, ]
# Define difference between myeloid inflammation and T effector
df$t2myel <- scale(df$motzer_myeloid_inflammation) -
    scale(df$motzer_t_effector)

p <- ggplot(df, aes(x = t2myel)) + 
    geom_histogram(aes(y = ..density..), fill = tx_palette[["lightblue"]], alpha = .8) + 
    geom_vline(xintercept = 0, col = tx_palette[["lightred"]], linetype = "dashed") + 
    labs(y = "Number of patients", x = "Myeloid inflammation - T effector ssGSEA scores")

ggsave(file.path(FIG_DIR, "histogram_teff2myel.pdf"), width = 55, height = 55, units = "mm")

df$tme_category <- ifelse(df$t2myel > 0, "Immunosuppressive",
    ifelse(df$t2myel < 0, "T effector", "Intermediate")
)

table(df$tme_category, df$patient)

transitions_df <- data.frame(
    patient = c(),
    s1 = c(), s2 = c(),
    tme1 = c(), tme2 = c(),
    transition = c()
)

for (patient in unique(df$patient)) {
    samples <- df$sample[df$patient == patient]
    pairs <- combn(samples, 2)
    for (i in 1:ncol(pairs)) {
        s1 <- pairs[1, i]
        s2 <- pairs[2, i]
        tme1 <- df$tme_category[df$sample == s1]
        tme2 <- df$tme_category[df$sample == s2]
        transition <- paste0(tme1, "-", tme2)
        if (transition == "Immunosuppressive-T effector") {
            tmp <- s1
            s1 <- s2
            s2 <- tmp
            transition <- "T effector-Immunosuppressive"
        }
        transitions_df <- rbind(
            transitions_df,
            data.frame(
                patient = patient, s1 = s1, s2 = s2,
                tme1 = tme1, tme2 = tme2, transition = transition
            )
        )
    }
}

# Clean all IDs to ensure clean matching of IDs across dataframes
annotation$sample <- clean_ids(annotation$sample)
clones_annotation$sample <- clean_ids(clones_annotation$sample)
transitions_df$s1 <- clean_ids(transitions_df$s1)
transitions_df$s2 <- clean_ids(transitions_df$s2)

# Prepare data to add information on whether sample contains terminal clone
clones_annotation <- clones_annotation[!is.na(clones_annotation$pat_terminal), ]
clones_annotation <- clones_annotation %>%
    group_by(sample) %>%
    summarise(pat_terminal = sum(pat_terminal == "Terminal") > 0) %>%
    mutate(pat_terminal = ifelse(pat_terminal, "Terminal", "No terminal"))

# COMPARE TRANSITIONS
table(transitions_df$transition)

# 1) See if transitions more frequently co-occurring with given event
# 2) Check directionality of the change

# BY TERMINAL CLONE -----------
get_coocurrence_transition(
    transitions_df, clones_annotation,
    "pat_terminal",
    is_driver = F
)

# No terminal > Terminal has more transitions (?)
fisher.test(matrix(c(54 + 104, 71, 3+5+21+5, 12+3), nrow = 2, byrow = T))

# No terminal > Terminal vs Terminal > No terminal
# fisher.test(matrix(c(12, 3, 7.5, 7.5), nrow = 2, byrow = T)) # 0.13
chisq.test(c(12,3), p = c(0.5, 0.5)) # 0.02


# BY SETD2 -------------
get_coocurrence_transition(
    transitions_df, annotation,
    "SETD2",
    is_driver = T
)

# SETD2 WT > WT
fisher.test(matrix(c(54 + 104, 71, 10 + 1 + 8 + 13, 15), nrow = 2, byrow = T))

# SETD2 WT > MT vs SETD2 MT > WT
# fisher.test(matrix(c(15, 0, 7.5, 7.5), nrow = 2, byrow = T)) # 0.002
chisq.test(c(15, 0), p = c(0.5, 0.5)) # 0.0001075

# By 9p loss -----------
get_coocurrence_transition(transitions_df, annotation, "loss_9p", TRUE)

# 9p WT > 9p loss
fisher.test(matrix(c(54 + 104, 71, 8+8+7+5, 25), nrow = 2, byrow = T)) # 0.03

# 9p WT > 9p loss vs 9p loss > WT
# fisher.test(matrix(c(24, 1, 12.5, 12.5), nrow = 2, byrow = T)) # 0.0002
chisq.test(c(24, 1), p = c(0.5, 0.5)) # 4.2-06 (p < 0.0001)
