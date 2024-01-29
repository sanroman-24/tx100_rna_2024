# Differential expression analysis 9p and 14q cells vs wildtype

rm(list = ls(all = TRUE))

# PACKAGES ----------------------------------------------------------------
library(Seurat)
library(MAST)
library(tidyverse)
library(here)
library(nlme)
library(ggpubr)
library(ggbeeswarm)
library(msigdbr)
library(fgsea)

# PATHS -------------------------------------------------------------------

BASE = here::here()
OUT_DIR = file.path(BASE, "analysis", "outputs", "scrna")
FIG_DIR = file.path(BASE, "analysis", "figures", "scrna") # TODO tmp
SEURAT_PATH = file.path(OUT_DIR, "chr3p_subset_qc_exp_filt_merged_clustered.rds")
META_DIR = file.path(BASE, "data", "meta")
HALLMARK_GROUPS_PATH <- file.path(META_DIR, "martinez_ruiz_2023_hallmark_gs_groups.txt")


# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

GSEADESeq2 <- function(DESeq2_results_df, pathways_info){
  # rank genes by -log FDR * sign(logFC). This allows to go from the most 
  # statistically significant upregulated to the most statistically significant downregulated
  
  # some modifications to input df
  DESeq2_results_df$padj = DESeq2_results_df$p_val_adj
  DESeq2_results_df$log2FoldChange = DESeq2_results_df$avg_log2FC
  DESeq2_results_df = rownames_to_column(DESeq2_results_df, "gene_name")
  DESeq2_results_df = DESeq2_results_df[!is.na(DESeq2_results_df$p_val) & 
                                          !is.na(DESeq2_results_df$log2FoldChange),]
  # to avoid issues getting values = Inf
  DESeq2_results_df$p_val[DESeq2_results_df$p_val == 0] = 1e-22
  ranking <- DESeq2_results_df %>% mutate(rank = -log10(p_val) * sign(log2FoldChange)) %>% # add rank
    arrange(desc(rank)) %>% # go from most significant upregulated to downregulated
    dplyr::select(gene_name, rank) 
  
  # create a vector with genes ordered by rank
  ranked_genes <- ranking$rank
  names(ranked_genes) <- ranking$gene_name
  # run and output GSEA
  fgsea(pathways_info, ranked_genes, nPermSimple = 10000)
}

subset_samples = function(srat, v, thresh = c(0.1,0.9, 5)){
  # get the fraction of cells within a sample that has alteration
  tb = srat@meta.data %>% group_by(Patient_id) %>% summarise(freq = mean(!!sym(v)), num = sum(!!sym(v)))
  keep_pts = tb$Patient_id[tb$freq > thresh[1] & tb$freq < thresh[2] & tb$num > thresh[3]]
  # only keep the single-cells from the patients that have subclonal event
  return(srat[,srat@meta.data$Patient_id %in% keep_pts])
}

seurat_dea = function(srat, v, ident.1, ident.2, outname, pathways){
  srat = NormalizeData(srat)
  # Find DE features using FindMarkers analysis
  Idents(srat) <- v
  de.markers <- FindMarkers(srat, ident.1 = "TRUE", ident.2 = "FALSE")
  # view results
  # head(de.markers)
  # write results
  saveRDS(de.markers, file.path(OUT_DIR, paste0(outname, "_seurat_dea.rds")))
  gsea = GSEADESeq2(de.markers, pathways)
  saveRDS(gsea, file.path(OUT_DIR, paste0(outname, "_seurat_gsea.rds")))
  return(gsea)
}

pseudobulk_dea = function(srat, v, thresh = c(0.1,0.9, 5), ident.1, ident.2, outname, pathways){
  # subset to only samples with certain fraction and total number of cells to study subclonal event
  # and avoid excessive weight into noisy inferCNV detected cells
  pseudo_srat = subset_samples(srat, v, thresh = thresh)
  pseudo_srat = AggregateExpression(pseudo_srat, assays = "RNA", return.seurat = T, group.by = c("Patient_id", v))
  
  # Differential expression analysis using DESeq2 on pseudobulk counts
  Idents(pseudo_srat) <- v
  de.markers <- FindMarkers(pseudo_srat, ident.1 = "TRUE", ident.2 = "FALSE", test.use = "DESeq2")
  # view results
  # head(de.markers)
  saveRDS(de.markers, file.path(OUT_DIR, paste0(outname, "_pseudo_deseq.rds")))
  
  gsea = GSEADESeq2(de.markers, pathways)
  saveRDS(gsea, file.path(OUT_DIR, paste0(outname, "_pseudo_gsea.rds")))
  return(gsea)
}

# LOAD DATA ---------------------------------------------------------------
gene_groups <- read_delim(HALLMARK_GROUPS_PATH, delim = "\t")
srat = readRDS(SEURAT_PATH)
srat = srat[,!is.na(srat$Disease_type) & srat$Disease_type == "ccRCC_sporadic" | srat$Publication == "Obradovic-etal", ]
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
pathways = split(x = hallmark_gene_sets$gene_symbol, f = hallmark_gene_sets$gs_name)

srat$is_9p = srat$Chr9p21.3_Status == "Loss"
srat$is_14q = srat$Chr14q23_31_Status == "Loss"


# DIFFERENTIAL EXPRESSION ANALYSIS 9P LOSS --------------------------------
# Seurat @ single-cell level

seurat_9p_gsea = seurat_dea(srat, "is_9p", "TRUE", "FALSE", "loss_9p_wt", pathways)
# Pseudobulk DEA when subclonal 9p loss
pseudo_9p_gsea = pseudobulk_dea(srat, "is_9p", c(0.1, 0.9, 10), "TRUE", "FALSE", "loss_9p_wt", pathways)
# DIFFERENTIAL EXPRESSION ANALYSIS 14Q LOSS -------------------------------
# Seurat @ single-cell level
seurat_14q_gsea = seurat_dea(srat, "is_14q", "TRUE", "FALSE", "loss_14q_wt", pathways)
# Pseudobulk DEA when subclonal 14q loss
pseudo_14q_gsea = pseudobulk_dea(srat, "is_14q", c(0.1, 0.9, 10), "TRUE", "FALSE", "loss_14q_wt", pathways)




# PLOT RESULTS ------------------------------------------------------------
# 9p loss
pseudo_9p_gsea$Hallmark <- str_remove(str_to_lower(pseudo_9p_gsea$pathway), "hallmark_")
pseudo_9p_gsea <- merge(pseudo_9p_gsea, gene_groups, by = "Hallmark")

pseudo_9p_gsea = pseudo_9p_gsea %>% as.data.frame %>% 
  mutate(pathway = ifelse(padj < 0.05, pathway, "")) %>% 
  mutate(alpha = ifelse(padj < 0.05, "dark", "soft")) %>% 
  mutate(logp = -log10(padj))  

p = plot_volcano(df = pseudo_9p_gsea, x_str = "NES", 
                 y_str = "logp", fill_str = "Functional_group", 
                 lab_str = "pathway", alpha_str = "alpha", 
                 lgd = "yes", xlim = c(-2,2))

save_ggplot(p, file.path(FIG_DIR, "dea_9p_scrna"), w = 100, h = 100)

# 14q loss
pseudo_14q_gsea$Hallmark <- str_remove(str_to_lower(pseudo_14q_gsea$pathway), "hallmark_")
pseudo_14q_gsea <- merge(pseudo_14q_gsea, gene_groups, by = "Hallmark")

pseudo_14q_gsea = pseudo_14q_gsea %>% as.data.frame %>% 
  mutate(pathway = ifelse(padj < 0.05, pathway, "")) %>% 
  mutate(alpha = ifelse(padj < 0.05, "dark", "soft")) %>% 
  mutate(logp = -log10(padj))  

p = plot_volcano(df = pseudo_14q_gsea, x_str = "NES", 
                 y_str = "logp", fill_str = "Functional_group", 
                 lab_str = "pathway", alpha_str = "alpha", 
                 lgd = "yes", xlim = c(-2,2))

save_ggplot(p, file.path(FIG_DIR, "dea_14q_scrna"), w = 100, h = 100)


# TSNE --------------------------------------------------------------------

srat_sub = subset_samples(srat, "is_9p", c(0.1, 0.9, 10))
srat_sub$cat = paste0(srat_sub$Patient_id, ":", srat_sub$is_9p)
Idents(srat_sub) = "is_9p"
srat_sub = RunTSNE(srat_sub)
TSNEPlot(srat_sub)

srat_sub = subset_samples(srat, "is_9p", c(0.1, 0.9, 10))
srat_sub$cat = paste0(srat_sub$Patient_id, ":", srat_sub$is_9p)
Idents(srat_sub) = "is_14q"
srat_sub = RunTSNE(srat_sub)
TSNEPlot(srat_sub)