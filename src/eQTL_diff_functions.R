library(data.table)

parse_pair <- function(pair) {
    match <- str_match(pair, "([^:]+):(.*$)")
    return(c(match[, 2], match[, 3]))
}

# Gets Es1 - Es2 in all pairs from same patient
# returns a Genes x pairs matrix
get_exp_dif <- function(vst) {
    vst <- vst[, !duplicated(colnames(vst))]
    genes <- rownames(vst)
    samples <- colnames(vst)
    pts <- get_pt(samples)
    # filter to patients with more than one sample
    pts <- table(pts) %>%
        {
            names(.)[. > 1]
        }
    # get differences in expression between pairs from same patient
    edif_matrix <- as.data.frame(rbindlist(lapply(unique(pts), function(pt) {
        print(pt)
        pairs <- combn(samples[grepl(pt, samples)], 2)
        pairs_edifs <- rbindlist(lapply(1:ncol(pairs), function(j) {
            s1 <- pairs[1, j]
            s2 <- pairs[2, j]
            pair <- paste0(s1, ":", s2)
            edifs <- sapply(genes, function(gene) {
                vst[rownames(vst) == gene, colnames(vst) == s1] -
                    vst[rownames(vst) == gene, colnames(vst) == s2]
            })
            o <- as.data.frame(t(edifs))
            o$pair <- pair
            return(o)
        }))
        return(pairs_edifs)
    })))
    pairs <- edif_matrix$pair
    edif_matrix <- edif_matrix[, !colnames(edif_matrix) %in% "pair"]
    edif_matrix <- as.matrix(edif_matrix)
    rownames(edif_matrix) <- pairs
    # output as gene x pairs
    return(t(edif_matrix))
}

# Gets CNs1 - CNs2 in all pairs from same patient
# where CN is cn at locus of given gene
# returns a Genes x pairs matrix
get_cn_dif <- function(cnmat, pairs) {
    genes <- rownames(cnmat)
    cnmat <- cnmat[, !duplicated(colnames(cnmat))]
    cn_dif_mat <- as.data.frame(rbindlist(lapply(genes, function(gene) {
        cn_difs <- sapply(pairs, function(pair) {
            samples <- parse_pair(pair)
            s1 <- samples[1]
            s2 <- samples[2]
            cn_dif <- cnmat[rownames(cnmat) == gene, colnames(cnmat) == s1] -
                cnmat[rownames(cnmat) == gene, colnames(cnmat) == s2]
            return(cn_dif)
        })
        o <- as.data.frame(t(cn_difs))
        return(o)
    })))
    rownames(cn_dif_mat) <- genes
    return(cn_dif_mat)
}

# Returns a vector of the patient to which each pair corresponds
get_pat_vector <- function(pairs) {
    return(get_pt(pairs))
}

# Gets Ds1 - Ds2 in all pairs from same patient
# where D is driver
# returns a vector with differences in driver for each sample
get_driver_dif <- function(annotation, driver, pairs) {
    annotation <- annotation[!duplicated(annotation$sample), ]
    annotation[[driver]] <- ifelse(annotation[[driver]] > 0, 1, 0)
    return(sapply(pairs, function(pair) {
        samples <- parse_pair(pair)
        s1 <- samples[1]
        s2 <- samples[2]
        return(annotation[[driver]][annotation$sample == s1] -
            annotation[[driver]][annotation$sample == s2])
    }))
}

# Gets Ps1 - Ps2 in all pairs from same patient
# where P is purity
get_pur_dif <- function(annotation, purity_v, pairs) {
    annotation <- annotation[!duplicated(annotation$sample), ]
    return(sapply(pairs, function(pair) {
        samples <- parse_pair(pair)
        s1 <- samples[1]
        s2 <- samples[2]
        return(annotation[[purity_v]][annotation$sample == s1] -
            annotation[[purity_v]][annotation$sample == s2])
    }))
}

# Gets distance between clones in two samples
get_clonal_dist <- function(clonal_dist, pairs) {
    clonal_dist$pair1 <- paste0(clonal_dist$smp1, ":", clonal_dist$smp2)
    clonal_dist$pair2 <- paste0(clonal_dist$smp2, ":", clonal_dist$smp1)
    return(sapply(pairs, function(pair) {
        samples <- parse_pair(pair)
        s1 <- samples[1]
        s2 <- samples[2]
        if (pair %in% clonal_dist$pair1 | pair %in% clonal_dist$pair2) {
            return(unique(c(
                # with unique we simply omit lack of match in one of entries
                clonal_dist$clonal_dist[clonal_dist$pair1 == pair],
                clonal_dist$clonal_dist[clonal_dist$pair2 == pair]
            )))
        }
        return(NA)
    }))
}

# Gets Ns1 - Ns2 in all pairs from same patient
# where N is annotation of whether sample is normal
# returns a vector with 1 for normal-prim, 0 if same, -1 if prim-normal
get_normal_dif <- function(annotation, is_normal_v, pairs) {
    annotation <- annotation[!duplicated(annotation$sample), ]
    return(sapply(pairs, function(pair) {
        samples <- parse_pair(pair)
        s1 <- samples[1]
        s2 <- samples[2]
        return(annotation[[is_normal_v]][annotation$sample == s1] -
            annotation[[is_normal_v]][annotation$sample == s2])
    }))
}

# Gets Ms1 - Ms2 in all pairs from same patient
# where N is annotation of whether sample is normal
# returns a vector with 1 for met-prim, 0 if same, -1 if prim-met
get_met_dif <- function(annotation, is_met_v, pairs) {
    annotation <- annotation[!duplicated(annotation$sample), ]
    return(sapply(pairs, function(pair) {
        samples <- parse_pair(pair)
        s1 <- samples[1]
        s2 <- samples[2]
        return(annotation[[is_met_v]][annotation$sample == s1] -
            annotation[[is_met_v]][annotation$sample == s2])
    }))
}
