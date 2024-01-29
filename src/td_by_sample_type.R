get_samples <- function(pt, annotation) {
    # get sample IDs of primaries and noprim from that patient
    prims <- annotation %>%
        filter(Patient == pt & is_primary == 1) %>%
        dplyr::select(sample) %>%
        as_vector()
    no_prims <- annotation %>%
        filter(Patient == pt & is_noprim == 1) %>%
        dplyr::select(sample) %>%
        as_vector()
    return(list(prims = prims, no_prims = no_prims))
}

# gets the I-TED prim-prim and prim-noprim from a given patient
# with prim samples `prims` and no prim samples `no_prims`
get_ith <- function(d_mat, prims, no_prims) {
    prim_prim_ITH <- c()
    prim_noprim_ITH <- c()
    for (sample1 in c(prims, no_prims)) {
        for (sample2 in prims) {
            if (sample1 != sample2) {
                ITH <- d_mat[rownames(d_mat) == sample1, colnames(d_mat) == sample2]
                if (sample1 %in% prims & (!paste0(sample1, sample2) %in% names(prim_prim_ITH)) & (!paste0(sample2, sample1) %in% names(prim_prim_ITH))) {
                    prim_prim_ITH[paste0(sample1, sample2)] <- ITH
                } else if ((!paste0(sample1, sample2) %in% names(prim_prim_ITH)) & (!paste0(sample2, sample1) %in% names(prim_prim_ITH))) {
                    prim_noprim_ITH[paste0(sample1, sample2)] <- ITH
                }
            }
        }
    }
    return(list(prim_prim = max(prim_prim_ITH), prim_noprim = max(prim_noprim_ITH)))
}


get_pair_dist <- function(d_mat, prims, no_prims) {
    out <- data.frame(sample1 = c(), sample2 = c(), t_d = c())
    for (sample1 in c(prims, no_prims)) {
        for (sample2 in prims) {
            if (sample1 != sample2) {
                ITH <- d_mat[rownames(d_mat) == sample1, colnames(d_mat) == sample2]
                out <- rbind(out, data.frame(sample1 = sample1, sample2 = sample2, t_d = ITH))
            }
        }
    }
    return(out)
}


# Returns data frame with prim-prim and prim-noprim heterogeneities for all the patients in the annotation
# patient IDs in `Patient` column annotation
# all the patient have ge 1 no prim sample and gt 1 prim sample
get_ited_by_tissue <- function(d_mat, annotation, prim_labs, noprim_labs) {
    # define binary variable that indicates if the sample is primary
    annotation$is_noprim <- ifelse(annotation$type_collapsed %in% c(noprim_labs), 1, 0)
    # define binary variable that indicates if the sample is from primary tumour
    annotation$is_primary <- ifelse(annotation$type_collapsed %in% c(prim_labs), 1, 0)

    # filter to patients with ge 1 no prim region and gt 1 prim region
    keep_patients <- annotation %>%
        dplyr::group_by(Patient) %>%
        # count how many primary and noprim samples
        dplyr::summarise(num_noprims = sum(is_noprim), num_prims = sum(is_primary)) %>%
        # get IDs of patients with ge 1 no prim samples and gt 1 prim samples
        filter(num_noprims >= 1 & num_prims > 1) %>%
        dplyr::select(Patient) %>%
        as_vector()

    annotation <- annotation[annotation$Patient %in% keep_patients, ]

    # dataframe where to store the final result
    out <- data.frame(Patient = c(), ITH_prim_prim = c(), ITH_prim_noprim = c())
    for (pt in unique(annotation$Patient)) {
        samples <- get_samples(pt, annotation)
        tmp_ith <- get_ith(d_mat, samples$prims, samples$no_prims)
        out <- rbind(
            out,
            data.frame(
                Patient = pt,
                ITH_prim_prim = tmp_ith$prim_prim,
                ITH_prim_noprim = tmp_ith$prim_noprim
            )
        )
    }

    return(out)
}

# Returns the transcriptional distance between all the possible pairs of
# samples where one of the samples is primary and the other is not
get_all_td_pairs_prim_noprim <- function(d_mat, annotation, prim_labs, noprim_labs) {
    # define binary variable that indicates if the sample is not a primary sample
    annotation$is_noprim <- ifelse(annotation$type_collapsed %in% c(noprim_labs), 1, 0)
    # define binary variable that indicates if the sample is from primary tumour (do not yet consider primary tumour here)
    annotation$is_primary <- ifelse(annotation$type_collapsed %in% c(prim_labs), 1, 0)

    # filter to patients with ge 1 no primary sample and gt 1 primary sample
    keep_patients <- annotation %>%
        dplyr::group_by(Patient) %>%
        dplyr::summarise(num_noprims = sum(is_noprim), num_prims = sum(is_primary)) %>%
        filter(num_noprims >= 1 & num_prims > 1) %>%
        dplyr::select(Patient) %>%
        as_vector()

    annotation <- annotation[annotation$Patient %in% keep_patients, ]

    # dataframe where to store the final result
    out <- data.frame(sample1 = c(), sample2 = c(), t_d = c())
    for (pt in unique(annotation$Patient)) {
        samples <- get_samples(pt, annotation)
        tmp <- get_pair_dist(d_mat, samples$prims, samples$no_prims)
        out <- rbind(out, tmp)
    }
    return(out)
}
