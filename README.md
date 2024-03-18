# Tracking non-genetic evolution from primary to metastatic ccRCC: TRACERx Renal 100

## Structure and content of this repository

This repository contains all the code and data necessariy to reproduce the results, figures (main and supplementary) and tables of the manuscript "*Tracking non-genetic evolution from primary to metastatic ccRCC: TRACERx Renal*".

Running the scripts in [analysis/scripts](https://github.com/sanroman-24/tx100_rna_2024/tree/main/analysis/scripts) reproduces all the analyses in the manuscript, with the exception of single-cell RNA-Seq results, for which data availability is controlled by corresponding authors. Auxiliar functions are located in [src](https://github.com/sanroman-24/tx100_rna_2024/tree/main/src). Input, intermediate and output data are provided, so that each script is runnable without the need to run previous scripts. 

We provide below an index of what scripts reproduce the different figures and tables in the manuscript.

## Linking code to manuscript figures (main and supplementary) and tables

The table below explicitly lists which code reproduces each figure in the manuscript:

| Figure type   | Figure number  | Code location |
| ------------- | -------------- | --------------| 
| Main          | 1              | [analysis/scripts/generate_oncoprint.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/generate_oncoprint.R)|
| Main          | 2a              | [analysis/scripts/umap.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 2b              | [analysis/scripts/estimate_transcriptional_ited.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/estimate_transcriptional_ited.R)|
| Main          | 2c              | [analysis/scripts/correlates_primary_ited.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/correlates_primary_ited.R)|
| Main          | 2d              | [analysis/scripts/survival_transcriptional_ited.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/survival_transcriptional_ited.R)|
| Main          | 3a              | [analysis/scripts/td_prim_norm.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_prim_norm.R)|
| Main          | 3b              | [analysis/scripts/td_prim_met.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_prim_met.R)|
| Main          | 3c              | [analysis/scripts/td_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_vs_clonal_dist.R)|
| Main          | 3d              | [analysis/scripts/td_metclone.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_metclone.R)|
| Main          | 3e              | [analysis/scripts/ssgsea_early_late.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/ssgsea_early_late.R)|
| Main          | 3f             | None, manually using graphical design software|
| Main          | 4a             | None, manually using graphical design software|
| Main          | 4b             | [analysis/scripts/paired_ssgsea_driver.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/paired_ssgsea_driver.R)|
| Main          | 4c             | [analysis/scripts/cgas_dea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_dea.R)|
| Main          | 4d             | [analysis/scripts/cgas_surv_tx.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_surv_tx.R)|
| Main          | 4e             | [analysis/scripts/cgas_surv_tcga.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_surv_tcga.R)|
| Main          | 5a             | [analysis/scripts/heatmap_tme_ith.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/heatmap_tme_ith.R)|
| Main          | 5b             | [analysis/scripts/correlates_primary_tme_ith.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/correlates_primary_tme_ith.R)|
| Main          | 5c             | [analysis/scripts/survival_tme_ith.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/survival_tme_ith.R)|
| Main          | 5d             | [analysis/scripts/ssgsea_early_late.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/ssgsea_early_late.R)|
| Main          | 5e             | [analysis/scripts/tme_by_evotype.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_by_evotype.R)|
| Main          | 6a             | [analysis/scripts/get_tcr_bcr_sim.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/get_tcr_bcr_sim.R)|
| Main          | 6b             | [analysis/scripts/TCR_BCR_similarity_survival.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/TCR_BCR_similarity_survival.R)|
| Main          | 6c             | [analysis/scripts/tcr_bcr_matches.Rmd](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_matches.Rmd)|
| Main          | 6d             | [analysis/scripts/tcr_bcr_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_vs_clonal_dist.R)|
| Main          | 6e             | [analysis/scripts/tcr_metclone.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_metclone.R)|
| Main          | 7a             | [analysis/scripts/fig7_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig7_HERV_analysis.R)|
| Main          | 7b             | [analysis/scripts/fig7_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig7_HERV_analysis.R)|
| Main          | 7c             | [analysis/scripts/fig7_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig7_HERV_analysis.R)|
| Main          | 7d             | [analysis/scripts/fig7_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig7_HERV_analysis.R)|
| Supplementary          | 1             | [analysis/scripts/umap.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Supplementary          | 2             | [analysis/scripts/lme_transcriptional_pc.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/lme_transcriptional_pc.R)|
| Supplementary          | 3             | [analysis/scripts/lme_transcriptional_pc.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/lme_transcriptional_pc.R)|
| Supplementary          | 4             | [analysis/scripts/umap.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Supplementary          | 5             | None, manually using graphical design software. Mock scatter plots generated in [src/td_mock_example.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/src/td_mock_example.R)|
| Supplementary          | 6             | [analysis/scripts/survival_transcriptional_ited.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/survival_transcriptional_ited.R)|
| Supplementary          | 7             | [analysis/scripts/prim_norm_td_clonal.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/prim_norm_td_clonal.R)|
| Supplementary          | 8             | None, manually using graphical design software|
| Supplementary          | 9             | [analysis/scripts/td_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_vs_clonal_dist.R)|
| Supplementary          | 10             | None, manually using graphical design software|
| Supplementary          | 11             | scRNA-Seq data is not publically available, hence code is not available|
| Supplementary          | 12             | [analysis/scripts/surv_9p_prolif.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/surv_9p_prolif.R)|
| Supplementary          | 13a             | [analysis/scripts/cgas_surv_tx.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_surv_tx.R)|
| Supplementary          | 13b             | [analysis/scripts/cgas_surv_tcga.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_surv_tcga.R)|
| Supplementary          | 14             | [analysis/scripts/surv_pan_cancer_slc19a1.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/surv_pan_cancer_slc19a1.R)|
| Supplementary          | 15             | [analysis/scripts/04_consensustme_vs_histo.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/04_consensustme_vs_histo.R)|
| Supplementary          | 16             | [analysis/scripts/tme_dist_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_dist_vs_clonal_dist.R)|
| Supplementary          | 17             | [analysis/scripts/TME_early_late.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/TME_early_late.R)|
| Supplementary          | 18             | [analysis/scripts/CD8_paired_drivers.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/CD8_paired_drivers.R)|
| Supplementary          | 19             | [analysis/scripts/TCR_BCR_similarity_survival.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/TCR_BCR_similarity_survival.R)|
| Supplementary          | 20             | [analysis/scripts/tcr_bcr_matches.Rmd](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_matches.Rmd)|
| Supplementary          | 21             | [analysis/scripts/tcr_bcr_matches.Rmd](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_matches.Rmd)|
| Supplementary          | 22             | [analysis/scripts/tcr_bcr_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_vs_clonal_dist.R)|
| Supplementary          | 23             | [analysis/scripts/fig7_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig7_HERV_analysis.R)|
| Supplementary          | 24             | [analysis/scripts/fig7_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig7_HERV_analysis.R)|
| Supplementary          | Table 1        | [analysis/scripts/normal_prim_ssgsea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/normal_prim_ssgsea.R)|
| Supplementary          | Table 1        | [analysis/scripts/prim_met_ssgsea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/prim_met_ssgsea.R)|
| Supplementary          | Table 3        | [data/meta/martinez_ruiz_2023_hallmark_gs_groups.txt](https://github.com/sanroman-24/tx100_rna_2024/blob/main/data/meta/martinez_ruiz_2023_hallmark_gs_groups.txt). Obtained from [Martínez-Ruiz et al, Nature 2023](https://www.nature.com/articles/s41586-023-05706-4#Sec8) |
| Supplementary          | Table 4        | [analysis/scripts/cgas_dea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_dea.R)|
| Supplementary          | Table 5        | [analysis/scripts/ssgsea_early_late.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/ssgsea_early_late.R)|


