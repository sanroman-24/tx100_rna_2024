# Tracking non-genetic evolution from primary to metastatic ccRCC: TRACERx Renal 100

## Structure and content of this repository

This repository contains all the code and data necessary to reproduce the results, figures (main and supplementary) and tables of the manuscript "*Tracking non-genetic evolution from primary to metastatic ccRCC: TRACERx Renal*".

Running the scripts in [analysis/scripts](https://github.com/sanroman-24/tx100_rna_2024/tree/main/analysis/scripts) reproduces all the analyses in the manuscript, with the exception of single-cell RNA-Seq results, for which data availability is controlled by corresponding authors. Helper functions are located in [src](https://github.com/sanroman-24/tx100_rna_2024/tree/main/src). Input, intermediate and output data are provided, so that each script is runnable without the need to run previous scripts. 

*Extensive analyses included during the last round of review are included in [analysis/scripts]*

We provide below an index of what scripts reproduce the different figures and tables in the manuscript.

## Linking code to manuscript figures (main and supplementary) and tables

The table below explicitly lists which code reproduces each figure in the manuscript:

| Figure type   | Figure number  | Code location |
| ------------- | -------------- | --------------| 
| Main          | 1a              | [analysis/scripts/umap.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Main          | 1b             | None, manually using graphical design software. Mock scatter plots generated in [src/td_mock_example.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/src/td_mock_example.R)|
| Main          | 1c              | [analysis/scripts/estimate_transcriptional_ited.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/estimate_transcriptional_ited.R)|
| Main          | 1d              | [analysis/scripts/correlates_primary_ited.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/correlates_primary_ited.R)|
| Main          | 1e              | None, manually using graphical design software |
| Main          | 1f              | [analysis/scripts/eqtl_td.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/eqtl_td.R)|
| Main          | 2a              | [analysis/scripts/td_prim_norm.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_prim_norm.R)|
| Main          | 2b              | [analysis/scripts/eqtl_prim_normal.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/eqtl_prim_normal.R)|
| Main          | 2c              | [analysis/scripts/td_prim_met.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_prim_met.R)|
| Main          | 2d              | [analysis/scripts/eqtl_prim_met.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/eqtl_prim_met.R)|
| Main          | 2e              | [analysis/scripts/td_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_vs_clonal_dist.R)|
| Main          | 2f              | [analysis/scripts/td_metclone.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_metclone.R)|
| Main          | 2g              | [analysis/scripts/ssgsea_early_late_purity.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/ssgsea_early_late_purity.R)|
| Main          | 2h             | None, manually using graphical design software |
| Main          | 3a             | None, manually using graphical design software |
| Main          | 3b             | [analysis/scripts/paired_ssgsea_driver.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/paired_ssgsea_driver.R)|
| Main          | 3c             | [analysis/scripts/het_hd_9p.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/het_hd_9p.R)|
| Main          | 3d             | [analysis/scripts/cgas_dea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_dea.R)|
| Main          | 3e             | None, from mIF results |
| Main          | 3f             | None, from mIF results |
| Main          | 3g             | [analysis/scripts/cgas_surv_tx.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_surv_tx.R)|
| Main          | 4a             | [analysis/scripts/tme_by_evotype.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_by_evotype.R)|
| Main          | 4b             | [analysis/scripts/heatmap_tme_ith.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/heatmap_tme_ith.R)|
| Main          | 4b             | [analysis/scripts/heatmap_tme_ith.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/heatmap_tme_ith.R)|
| Main          | 4c             | [analysis/scripts/tme_transitions.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_transitions.R)|
| Main          | 4d             | [analysis/scripts/tme_transitions.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_transitions.R)|
| Main          | 4e             | [analysis/scripts/tme_transitions.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_transitions.R)|
| Main          | 4f             | [analysis/scripts/tme_transitions.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_transitions.R)|
| Main          | 5a             | [analysis/scripts/get_tcr_bcr_sim.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/get_tcr_bcr_sim.R)|
| Main          | 5b             | [analysis/scripts/TCR_BCR_similarity_survival.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/TCR_BCR_similarity_survival.R)|
| Main          | 5c             | [analysis/scripts/tcr_bcr_matches.Rmd](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_matches.Rmd)|
| Main          | 5d             | [analysis/scripts/tcr_bcr_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_vs_clonal_dist.R)|
| Main          | 5e             | [analysis/scripts/tcr_metclone.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_metclone.R)|
| Main          | 6a             | [analysis/scripts/fig6_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig6_HERV_analysis.R)|
| Main          | 6b             | [analysis/scripts/fig6_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig6_HERV_analysis.R)|
| Main          | 6c             | [analysis/scripts/fig6_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig6_HERV_analysis.R)|
| Main          | 6d             | [analysis/scripts/fig6_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig6_HERV_analysis.R)|
| Supplementary          | 1              | [analysis/scripts/generate_oncoprint.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/generate_oncoprint.R)|
| Supplementary          | 2a,b             | [analysis/scripts/umap.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/umap.R)|
| Supplementary          | 2c             | [analysis/scripts/intra_vs_inter.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/intra_vs_inter.R)|
| Supplementary          | 3a             | None, manually using graphical design software. Mock scatter plots generated in [src/td_mock_example.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/src/td_mock_example.R)|
| Supplementary          | 3b             |[analysis/scripts/ited_top500_vs_all.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/ited_top500_vs_all.R) |
| Supplementary          | 4             | [analysis/scripts/survival_transcriptional_ited.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/survival_transcriptional_ited.R)|
| Supplementary          | 5a             | [analysis/scripts/correlates_primary_ited.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/correlates_primary_ited.R)|
| Supplementary          | 5b             | [analysis/scripts/eqtl_td.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/eqtl_td.R)|
| Supplementary          | 6             | [analysis/scripts/eqtl_td.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/eqtl_td.R)|
| Supplementary          | 7             | [analysis/scripts/td_prim_norm.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_prim_norm.R)|
| Supplementary          | 8             | None, manually using graphical design software|
| Supplementary          | 9             | [analysis/scripts/td_metclone.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_metclone.R)|
| Supplementary          | 10             | [analysis/scripts/td_metclone.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/td_metclone.R)|
| Supplementary          | 11             | None, manually using graphical design software|
| Supplementary          | 12             | [analysis/scripts/paired_ssgsea_driver.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/paired_ssgsea_driver.R)|
| Supplementary          | 13             | [analysis/scripts/het_hd_9p.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/het_hd_9p.R)|
| Supplementary          | 14             | scRNA-Seq data is not publically available, hence code is not re-runnable|
| Supplementary          | 15             | [analysis/scripts/surv_9p_prolif.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/surv_9p_prolif.R)|
| Supplementary          | 16             | None, from mIF results |
| Supplementary          | 17             | [analysis/scripts/enpp1_slc19a1_tme.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/enpp1_slc19a1_tme.R)|
| Supplementary          | 18a             | [analysis/scripts/cgas_surv_tx.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_surv_tx.R)|
| Supplementary          | 18b,18c             | [analysis/scripts/cgas_surv_tcga.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_surv_tcga.R)|
| Supplementary          | 19             | [analysis/scripts/surv_pan_cancer_slc19a1.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/surv_pan_cancer_slc19a1.R)|
| Supplementary          | 20             | [analysis/scripts/04_consensustme_vs_histo.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/04_consensustme_vs_histo.R)|
| Supplementary          | 21a             | [analysis/scripts/estimate_tme_ith.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/estimate_tme_ith.R)|
| Supplementary          | 21b             | [analysis/scripts/tme_dist_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_dist_vs_clonal_dist.R)|
| Supplementary          | 21c,d             | [analysis/scripts/survival_tme_ith.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/survival_tme_ith.R)|
| Supplementary          | 22a             | [analysis/scripts/tme_transitions.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tme_transitions.R)|
| Supplementary          | 22b             | [analysis/scripts/ssgsea_early_late_purity.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/ssgsea_early_late_purity.R)|
| Supplementary          | 22c,d             | [analysis/scripts/paired_ssgsea_driver.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/paired_ssgsea_driver.R)|
| Supplementary          | 23             | [analysis/scripts/CD8_paired_drivers.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/CD8_paired_drivers.R)|
| Supplementary          | 24             | [analysis/scripts/TCR_BCR_similarity_survival.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/TCR_BCR_similarity_survival.R)|
| Supplementary          | 25             | [analysis/scripts/tcr_bcr_matches.Rmd](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_matches.Rmd)|
| Supplementary          | 26             | [analysis/scripts/tcr_bcr_matches.Rmd](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_matches.Rmd)|
| Supplementary          | 27             | [analysis/scripts/tcr_bcr_vs_clonal_dist.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_vs_clonal_dist.R)|
| Supplementary          | 28             | [analysis/scripts/tcr_bcr_cor.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_cor.R)|
| Supplementary          | 29             | [analysis/scripts/fig6_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig6_HERV_analysis.R)|
| Supplementary          | 30             | [analysis/scripts/herv_cnas.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/herv_cnas.R)|
| Supplementary          | 31             | [analysis/scripts/fig6_HERV_analysis.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/fig6_HERV_analysis.R)|
| Supplementary          | Table 1        | [analysis/scripts/normal_prim_ssgsea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/normal_prim_ssgsea.R)|
| Supplementary          | Table 2        | [analysis/scripts/normal_prim_ssgsea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/normal_prim_ssgsea.R)|
| Supplementary          | Table 3        | [analysis/scripts/prim_met_ssgsea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/prim_met_ssgsea.R)|
| Supplementary          | Table 4        | [data/meta/martinez_ruiz_2023_hallmark_gs_groups.txt](https://github.com/sanroman-24/tx100_rna_2024/blob/main/data/meta/martinez_ruiz_2023_hallmark_gs_groups.txt). Obtained from [Martínez-Ruiz et al, Nature 2023](https://www.nature.com/articles/s41586-023-05706-4#Sec8) |
| Supplementary          | Table 5        | [analysis/scripts/ssgsea_early_late_purity.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/ssgsea_early_late_purity.R)|
| Supplementary          | Table 6        | [analysis/scripts/cgas_dea.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/cgas_dea.R)|
| Supplementary          | Table 7        | [analysis/scripts/tcr_bcr_matches.Rmd](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_matches.Rmd)|
| Supplementary          | Table 8        | [analysis/scripts/tcr_bcr_matches.Rmd](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/tcr_bcr_matches.Rmd)|
| Supplementary          | Table 9        | [analysis/scripts/top500_genes_variation.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/top500_genes_variation.R)|