get_monoclonal_regions_df <- function(clones_annotation) {
    smps <- unique(clones_annotation$sample)
    out_df <- data.frame(smp = c(), clone = c())
    for (smp in smps) {
        term_clones <- get_terminal_clones(clones_annotation, smp)
        if (length(term_clones) == 1) {
            clone <- paste0(str_remove(smp, "_.*$"), "-", term_clones)
            out_df <- rbind(out_df, data.frame(smp = smp, clone = clone))
        }
    }
    return(out_df)
}

get_terminal_clones <- function(clones_annotation, smp) {
    smp_clones_annotation <- clones_annotation[clones_annotation$sample == smp, ]
    term_clones <- smp_clones_annotation$clone[!smp_clones_annotation$clone %in% smp_clones_annotation$parent_clone]
    return(term_clones)
}

get_dist <- function(smp1, smp2, clones_annotation) {
    if (!smp2 %in% clones_annotation$sample | !smp1 %in% clones_annotation$sample) {
        return(NA)
    }
    smp2_clones <- get_terminal_clones(clones_annotation, smp2)

    dists <- c()
    for (clone in smp2_clones) {
        dists <- c(dists, get_clonal_dist(clones_annotation, clone, smp1, smp2, 0))
    }
    return(min(dists))
}


get_clonal_dist <- function(clones_annotation, clone, smp1, smp2, curr_dist) {
    smp1_clones <- clones_annotation[clones_annotation$sample == smp1, ]$clone
    if (clone %in% smp1_clones | clone == "GL") {
        return(curr_dist)
        # dist_term_smp1 = get_dist_to_terminal(clones_annotation, smp1, c(clone), c(0))
        # dist_term_smp2 = get_dist_to_terminal(clones_annotation, smp2, c(clone), curr_dist)
    } else {
        par_clone <- clones_annotation[clones_annotation$clone == clone & clones_annotation$sample == smp2, ]$parent_clone
        curr_dist <- curr_dist + 1
        get_clonal_dist(clones_annotation, par_clone, smp1, smp2, curr_dist)
    }
}

get_dist_to_terminal <- function(clones_annotation, smp, curr_clones, curr_dists) {
    smp_clones <- clones_annotation[clones_annotation$sample == smp, ]
    clonal_dists <- c(0)
    clonal_dist <- traverse_clones(smp_clones, curr_clones, curr_dists, clonal_dists)
    return(clonal_dist + curr_dists)
}

get_met_prim_dist <- function(met_smp, prim_smp, clones_annotation) {
    if (!prim_smp %in% clones_annotation$sample | !met_smp %in% clones_annotation$sample) {
        return(NA)
    }
    prim_clones <- get_terminal_clones(clones_annotation, prim_smp)
    dists <- c()
    for (clone in prim_clones) {
        dists <- c(dists, get_clonal_dist(clones_annotation, clone, met_smp, prim_smp, 0))
    }
    return(max(dists))
}

calculate_clonal_age <- function(clones_annotation, pat, curr_clone, age) {
    if (str_detect(curr_clone, "GL")) {
        return(age)
    }
    if (!curr_clone %in% clones_annotation$clone_code) {
        return(NA)
    }
    parent_clone <- unique(clones_annotation[clones_annotation$clone_code == curr_clone, ]$parent_clone)
    curr_clone <- paste0(pat, "-", parent_clone)
    age <- age + 1
    return(calculate_clonal_age(clones_annotation, pat, curr_clone, age))
}

get_terminal_clones_pt <- function(clones_annotation, pat) {
    pat_clones_annotation <- clones_annotation[clones_annotation$patient == pat, ]
    term_clones <- pat_clones_annotation$clone[!pat_clones_annotation$clone %in% pat_clones_annotation$parent_clone]
    return(term_clones)
}

get_dist2mrca <- function(smp, clones_annotation) {
    if (!smp %in% clones_annotation$sample) {
        return(NA)
    }
    smp_clones <- get_terminal_clones(clones_annotation, smp)

    dists <- c()
    for (clone in smp_clones) {
        dists <- c(dists, get_branch_length(clones_annotation, clone, smp, 0))
    }
    return(min(dists))
}

get_branch_length <- function(clones_annotation, clone, smp, curr_dist) {
    if (clone == "GL") {
        return(curr_dist)
    } else {
        par_clone <- clones_annotation$parent_clone[
            clones_annotation$clone == clone &
                clones_annotation$sample == smp
        ]
        curr_dist <- curr_dist + 1
        get_branch_length(clones_annotation, par_clone, smp, curr_dist)
    }
}
