## module load r/4.5.0
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(glue)
library(bigreadr)

## set path
proj_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(proj_path)
ref_path <- "/ix1/wchen/xiangyu/Ref_data/NicheNet/"

source("code/FUNCTION/multinichenet/multinichnetr.R")
## load in LR reference
organism = "Human"
ligand_target_matrix <- readRDS(glue("{ref_path}{organism}/ligand_target_matrix.rds"))
lr_network <- readRDS(glue("{ref_path}{organism}/lr_network.rds"))
weighted_networks <- readRDS(glue("{ref_path}{organism}/weighted_networks.rds"))

## 0. set parameters
data_path <- "03_output/02_clean/"
out_path <- "03_output/05_Interaction/multinichenetr/"
# 0.1 set data info
sample_id <- "Sample_ID_exp"
group_id <- "Condition"
celltype_id <- "ann_level5_final"
covariates <- NA # provide DEGs directly 
batches <- NA

# 0.2 set cell type of interest
ct = "CD4_CD103_TRM_Th17"
message(ct)
#
ct_level <- readRDS("03_output/03_clustering/ct_levels_refine.rds")
all_ct <- ct_level[[celltype_id]]
receivers_oi <- senders_oi <- all_ct[grep("^CD4_CD103_TRM_Th17", all_ct)]
##
de_method <- "DESeq2"
message(de_method)

# for (type in c("TI", "Colon")) {
lapply(c("TI", "Colon"), function(type){
  # type = "TI"
  message(type)
  # # 0.3 set contrast
  # for (contrast_name in c("II_vs_NN", "NU_vs_NN")) {
  lapply(c("II_vs_NN", "NU_vs_NN"), function(contrast_name){
    
    # contrast_name = "II_vs_NN"
    message(glue("{ct}: {contrast_name}"))
    conditions_keep <- strsplit(contrast_name, "_vs_")[[1]]
    contrasts_oi <- c(paste(conditions_keep, collapse = "-"),
                      paste(rev(conditions_keep), collapse = "-"))
    contrast_tbl <- tibble(contrast = contrasts_oi, 
                           group = conditions_keep)
    ## 1. load data
    meta_path <- glue("{data_path}{type}/{ct}/metadata.tsv")
    tenX_path <- glue("{data_path}{type}/{ct}/RNA/")
    sce_obj <- create.from.tenX(tenX_path = tenX_path,
                                meta_path = meta_path,
                                sample_id = sample_id,
                                group_id = group_id,
                                celltype_id = celltype_id,
                                covariates = covariates,
                                batch = batch)
    ## 2. set DE info
    de_path <- glue("03_output/04_Diff/DEG/{type}/{ct}/{de_method}/")
    if (de_method == "MAST") {
      defile <- glue("{de_path}/lmm_{celltype_id}_{contrast_name}.rds")
    } else {
      defile <- glue("{de_path}/pseudo_bulk_{celltype_id}_{contrast_name}.rds")
    }
    celltype_de <- format.deg(defile_list = defile,
                              gene_col = "Term",
                              ct_col = "cluster",
                              fc_col = "log2FC",
                              stat_col = "stat",
                              p_col = "pvalue",
                              padj = "padj",
                              contrasts_oi = contrasts_oi)
    ## 3. run 
    multinichenet_output <- run.multinichenet(sce = sce_obj,
                                              sample_id = sample_id,
                                              group_id = group_id,
                                              celltype_id = celltype_id,
                                              covariates = covariates,
                                              batches = batches,
                                              senders_oi = senders_oi,
                                              receivers_oi = receivers_oi,
                                              contrast_tbl = contrast_tbl,
                                              celltype_de = celltype_de,
                                              ligand_target_matrix = ligand_target_matrix,
                                              lr_network = lr_network,
                                              weighted_networks = weighted_networks,
                                              min_cells = 10,
                                              min_sample_prop = 0.5,
                                              fraction_cutoff = 0.05,
                                              logFC_threshold = 0.25,
                                              p_val_threshold = 0.05,
                                              p_val_adj = T,
                                              top_n_target = 250,
                                              n_cores = 10)
    
    system(glue("mkdir -p {out_path}{type}/{ct}/{de_method}/{celltype_id}/"))
    saveRDS(multinichenet_output, glue("{out_path}{type}/{ct}/{de_method}/{celltype_id}/multinichenet_output_{contrast_name}.rds"))
    
    return(glue("{out_path}{type}/{ct}/{de_method}/{celltype_id}/multinichenet_output_{contrast_name}.rds"))
  })
  
})

# }
# }
