# helper functions to run eQTL analysis

get_cntab <- function(smp, cndir, as.gr = TRUE, chr_v = "Chromosome", start_v = "Start.bp", end_v = "End.bp", cn_v = "expected_cn") {
    if (str_detect(smp, "-")) {
        smp <- str_replace(smp, "-", "_")
    }
    pt <- str_extract(smp, "K\\d{3}")
    fp <- list.files(cndir, pattern = pt, full = T)
    if (length(fp) < 1) {
        print(paste0("No CN info for patient ", pt))
        return(NA)
    }
    cntab <- read_delim(fp, show_col_types = FALSE)
    cntab$sample <- cnid2smp(pt, cntab$sample)
    if (!any(str_detect(cntab$sample, smp))) {
        print(paste0("Check sample ", smp))
    }
    if (!smp %in% cntab$sample) {
        return(NA)
    }
    cntab <- cntab[cntab$sample == smp, ]
    if (as.gr) {
        o <- GRanges(
            seqnames = cntab[[chr_v]],
            IRanges(
                start = cntab[[start_v]],
                end = cntab[[end_v]]
            )
        )
        o$cn <- cntab[[cn_v]]
        return(o)
    }
    return(cntab)
}

get_cngenes <- function(genes.gr, cntab.gr) {
    ovs <- findOverlaps(genes.gr, cntab.gr)
    cns <- sapply(1:length(genes.gr), function(i) {
        j <- ovs@to[ovs@from == i]
        if (length(j) < 1) {
            return(NA)
        }
        return(mean(cntab.gr$cn[j]))
    })
    names(cns) <- genes.gr$id
    return(cns)
}

get_cnmat <- function(smps, cndir, genes.gr) {
    o <- matrix(NA, nrow = length(genes.gr), ncol = length(smps))
    colnames(o) <- smps
    rownames(o) <- genes.gr$id
    for (i in 1:ncol(o)) {
        print(i)
        smp <- smps[i]
        cntab.gr <- get_cntab(smp, cndir)
        if (length(cntab.gr) < 2) {
            o[, i] <- cn
            next()
        }
        cn <- get_cngenes(genes.gr, cntab.gr)
        o[, i] <- cn
    }
    return(o)
}

get_pur_vector <- function(smps, annotation) {
    annotation$sample <- clean_ids(annotation$sample)
    return(sapply(smps, function(smp) {
        annotation$purity[annotation$sample == smp][1]
    }))
}

get_pat_vector <- function(smps) {
    return(str_extract(smps, "K\\d{3}"))
}

filter_low <- function(counts, thr = 0.75) {
    idx <- which((rowSums(counts != 0) / ncol(counts)) > thr)
    return(rownames(counts)[idx])
}

run_lme_eqtl <- function(exp, ..., fixed = c("cn"), random = c("pt")) {
    l <- list(...)
    df <- as.data.frame(l)
    df$e <- exp
    is_ok <- complete.cases(df)
    if (mean(is_ok) < 0.5) {
        return(list(model = NA, p = NA))
    }
    df <- df[is_ok, ]
    random_formula <- as.formula(glue::glue("~ 1 | {random}"))
    random_formula <- lapply(random, function(v) ~ 1)
    names(random_formula) <- random
    covariates <- paste(fixed, collapse = " + ")
    fixed_formula <- as.formula(glue::glue("e ~ {covariates}"))
    full_model <- nlme::lme(fixed_formula, random = random_formula, data = df)
    null_model <- nlme::lme(e ~ 1, random = random_formula, data = df)
    a <- anova(null_model, full_model)
    p <- a[2, 9]
    return(list(model = full_model, p = p))
}

# pass genes only if we want it in output table
pull_eqtl_stats <- function(lfits, x, genes = NULL, na.rm = T) {
    # Assumes null model should be exp[i,] ~ 1 + (1|Patient)
    l <- lapply(1:length(lfits), function(i) {
        # handle cases where we could not fit model
        if (length(lfits[[i]]) == 1){
            return(data.frame(gene = genes[i], t = NA, p = NA, model_p = NA))
        }
        model <- lfits[[i]]$model 
        if (length(model) == 1) {
            return(data.frame(gene = genes[i], t = NA, p = NA, model_p = NA))
        }
        
        # get stats of LME
        r <- summary(lfits[[i]]$model) %>%
            {
                .$tTable
            }
        idx <- na.omit(match(x, rownames(r)))
        t <- NA; t[x %in% rownames(r)] = r[idx, "t-value"]
        p <- NA; p[x %in% rownames(r)] = r[idx, "p-value"]
        model_p <- lfits[[i]]$p
        return(data.frame(gene = genes[i], t = t, p = p, model_p = model_p))
    })
    o <- as.data.frame(data.table::rbindlist(l))
    if (na.rm) {
        o <- o[complete.cases(o), ]
    }
    return(o)
}


pull_all_eqtl_stats <- function(lfits, variables, genes, na.rm = T){
    o <- rbindlist(lapply(variables, function(x){
        df_stats <- pull_eqtl_stats(lfits, x, genes, na.rm)
        df_stats$fdr <- p.adjust(df_stats$p, method = "fdr")
        df_stats$variable <- x
        return(df_stats)
    }))
    return(as.data.frame(o))
}

get_proportion_explained <- function(lfits) {
    o <- sapply(1:length(lfits), function(i) {
        model <- lfits[[i]]$model
        # handle cases where we could not fit model
        if (length(model) == 1) {
            return(NA)
        }
        r_squared <- MuMIn::r.squaredGLMM(model)
        return(r_squared[2])
    })
    return(o)
}
