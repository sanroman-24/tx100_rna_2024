# Get subclonal subclonal SCNA 

get_pat_subSCNA <- function(cndist_pairs) {
  cndist_pairs %>%
    mutate(patient = str_extract(s1, "[A-Z]*_K\\d{3}")) %>%
    dplyr::group_by(patient) %>%
    dplyr::summarise(subclonal_scna = median(fga_subclonal_wgd_corrected)) %>%
    mutate(patient = str_remove(patient, "^[A-Z]_"))
}