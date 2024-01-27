library(energy)
library(lsa)
library(tidyverse)

# Returns the distance between two vectors with numerical data using function f
# With previous filtering of elements in both vectors to specific targets
#' @param c_s1 counts in the first sample
#' @param c_s2 counts in the second sample
#' @param t targets on which to calculate the distance
#' @param f function to use to calculate distance
get_dist <- function(c_s1, c_s2, t = c(), f) {
  f <- match.fun(f)
  if (length(t) > 0) {
    c_s1 <- c_s1[names(c_s1) %in% t]
    c_s2 <- c_s2[names(c_s2) %in% t]
  }
  return(1 - f(c_s1, c_s2)[1])
}

# Returns a pairwise distance matrix from a numerical matrix of features x samples
# subsetting the features to specific targets
#' @param counts matrix with numerical features
#' @param t target features on which to calculate distance
#' @param f function to use to calculate pairwise distances
get_dist_matrix <- function(counts, t, f = "dcor") {
  samps <- colnames(counts)
  all_pairs <- combn(colnames(counts), 2)
  out <- matrix(NA, nrow = length(samps), ncol = length(samps))
  rownames(out) <- samps
  colnames(out) <- samps
  for (pair in 1:ncol(all_pairs)) {
    s1 <- all_pairs[1, pair]
    s2 <- all_pairs[2, pair]
    c_s1 <- counts[, colnames(counts) == s1]
    c_s2 <- counts[, colnames(counts) == s2]
    t_dist <- get_dist(c_s1, c_s2, t, f)
    out[rownames(out) == s1, colnames(out) == s2] <- t_dist
    out[rownames(out) == s2, colnames(out) == s1] <- t_dist
  }
  return(out)
}

# Returns a per-patient estimate of intratumour heterogeneity
# by applying function summ_f to the distances between samples
# from given patient
#' @param d_mat pairwise distance matrix
#' @param pat_id patient for whom to estimate ITH
#' @param summ_f function to estimate patient ITH
summarise_ited_per_patient <- function(d_mat, pat_id, summ_f = "median") {
  d_mat <- d_mat[str_detect(rownames(d_mat), pat_id), str_detect(colnames(d_mat), pat_id)]
  s_pairs <- combn(colnames(d_mat), 2)
  t_d <- c()
  for (pair in 1:ncol(s_pairs)) {
    s1 <- s_pairs[1, pair]
    s2 <- s_pairs[2, pair]
    t_d <- c(t_d, d_mat[rownames(d_mat) == s1, colnames(d_mat) == s2])
  }
  f <- match.fun(summ_f)
  return(f(t_d))
}

# Wrapper to apply summarise_ited function to several patients
#' @param d_mat pairwise distance matrix
#' @param pat_id patients for whom to estimate ITH
#' @param summ_f function to estimate patient ITH
summarise_ited <- function(d_mat, pat_ids, summ_f) {
  return(sapply(pat_ids, function(pat_id) summarise_ited_per_patient(d_mat, pat_id, summ_f)))
}


# Same as 2 functions above but with similarity instead of distance matrix
# Just for readability
summarise_sim_per_patient <- function(sim_matrix, pat_id, summ_f){
  sim_matrix = sim_matrix[str_detect(rownames(sim_matrix), pat_id), str_detect(colnames(sim_matrix), pat_id)]
  s_pairs = combn(colnames(sim_matrix), 2)
  t_d = c()
  for (pair in 1:ncol(s_pairs)){
    s1 = s_pairs[1, pair]
    s2 = s_pairs[2, pair]
    t_d = c(t_d, sim_matrix[rownames(sim_matrix) == s1, colnames(sim_matrix) == s2])
  }
  f = match.fun(summ_f)
  return(f(t_d))
}

summarise_sim <- function(sim_matrix, pat_ids, summ_f){
  return(sapply(pat_ids, function(pat_id) summarise_sim_per_patient(sim_matrix, pat_id, summ_f)))
}
