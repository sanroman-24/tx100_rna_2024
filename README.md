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
| Main          | 2a              | [analysis/scripts/estimate_transcriptional_ited.R.R](https://github.com/sanroman-24/tx100_rna_2024/blob/main/analysis/scripts/estimate_transcriptional_ited.R.R)|




