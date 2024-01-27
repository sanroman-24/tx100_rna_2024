assign_expression_to_clones = function(dat, regions_per_clone_df){
  if (!all(regions_per_clone_df$smp %in% colnames(dat))){
    stop("all the regions of the annotation should be in the data from which to get the data")
  }
  out_exp = matrix(NA, nrow = nrow(dat))
  rownames(out_exp) = rownames(dat)
  clones = unique(regions_per_clone_df$clone)
  for (clone in clones){
    clone_regions = regions_per_clone_df[regions_per_clone_df$clone == clone, ]$smp
    if (length(clone_regions) > 1){
      clone_exp = rowMeans(dat[,which(colnames(dat) %in% clone_regions)])
    } else {
      clone_exp = dat[,colnames(dat) == clone_regions]
    }
    out_exp = cbind(out_exp, clone_exp)
    colnames(out_exp) = c(colnames(out_exp)[-ncol(out_exp)], clone)
  }
  out_exp = as.matrix(out_exp[,-1])
}

assign_purity_to_clones = function(annotation, regions_per_clone_df){
  if (!all(regions_per_clone_df$smp %in% annotation$sample)){
    stop("all the regions of the annotation should be in the data from which to get the data")
  }
  out_pur = c()
  clones = unique(regions_per_clone_df$clone)
  for (clone in clones){
    clone_regions = regions_per_clone_df[regions_per_clone_df$clone == clone, ]$smp
    if (length(clone_regions) > 1){
      clone_pur = mean(annotation[annotation$sample %in% clone_regions, ]$purity, na.rm = TRUE)
    } else {
      clone_pur = annotation[annotation$sample %in% clone_regions, ]$purity
    }
    out_pur = c(out_pur, clone_pur)
    names(out_pur) = c(names(out_pur)[-length(out_pur)], clone)
  }
  out_pur
}