library(nlme)
library(dplyr)

# little wrapper to quickly clean IDs in the annotation to match IDs when they have only _
clean_ids <- function(ids){
  ids %>%
  str_replace("K328R19", "K328-R19") %>%
  str_replace("-", "_") %>%
  str_replace("K328_T1", "K328_THR1") %>%
  str_replace("K245_T1", "K245_THR1")
}

get_pt <- function(smp){
  return(str_extract(smp, "K\\d{3}"))
}

get_gene_coordinates <- function(hgnc_genes, filt_chr = T, filt_dup = T){
  m <- useMart('ensembl', dataset='hsapiens_gene_ensembl') 
  df <- getBM(mart=m, attributes=c('hgnc_symbol', 'chromosome_name',
                                   'start_position', 'end_position'),
            filters='hgnc_symbol', values=hgnc_genes) 
  if (filt_chr){
    df <- df[df$chromosome %in% 1:22,]
  }
  if (filt_dup){
    df <- df[!duplicated(df$hgnc_symbol),]
  }
  return(df)
}

# from ID like R1d1xx1 or R4dna1xx1, get only R1 or R4
cnid2smp <- function(pt, ids){
    smps <- sapply(ids, function(id){
        e <- str_locate(id, "\\d+")[1,2]
        smp <- str_sub(id, 1, e)
        return(paste0(pt, "_", smp))
    })
    return(smps)
}

# Returns linear-mixed effect model with variables contained in dataframe
#' @param Y string with name of dependent variable
#' @param X string with the name of independent variable(s)
#' @param random_covariate name of the random variable to correct for in LME (usually patient ID)
#' @param df dataframe containing all the values for variables in the model
run_lme <- function(Y, X, random_covariate, df) {
    df <- df[complete.cases(df[X]), ]
    form <- as.formula(paste0(Y, " ~ ", paste0(X, collapse = " + ")))
    random_form <- as.formula(paste0("~ 1 | ", random_covariate))
    return(lme(form, random = random_form, data = df))
}

# Returns % of variance explained by each of the variables in the model
#' @param fm_str formula of linear regression as string
#' @param df dataframe containing all the values for variables in the model
get_perc_variation <- function(fm_str, df){
  fm = as.formula(fm_str)
  lmod = lm(fm, data = df)
  af = anova(lmod)
  afs = af$`Sum Sq`
  af$perc_var = afs / sum(afs) * 100
  return(af)
}

find_subclonal_alteration <- function(df, gene) {
    tab <- df %>%
        group_by(Patient) %>%
        dplyr::summarise(
            n_wt = sum(!!sym(gene) == 0),
            n_mt = sum(!!sym(gene) > 0)
        )
    pts_subclonal <- tab[tab$n_wt > 0 & tab$n_mt > 0, ]
    return(df[df$Patient %in% pts_subclonal$Patient, ])
}

df_to_ggpaired = function(df, pat_var, cond, var_str){
  df = df %>% 
  group_by(!!sym(pat_var), !!sym(cond)) %>% 
    dplyr::summarise(var = mean(!!sym(var_str)))
  keep_pats = df[[pat_var]][duplicated(df[[pat_var]])]
  df = df[df[[pat_var]] %in% keep_pats, ] %>% 
    pivot_wider(names_from = cond, values_from = var)
  return(df)
}
