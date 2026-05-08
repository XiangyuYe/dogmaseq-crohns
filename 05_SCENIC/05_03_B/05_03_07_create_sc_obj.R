#####################################################
load.AUC <- function(h5mu_file = NULL, dir = NULL){
  
  auc_list <- h5read(h5mu_file, dir)
  auc_df <- as.data.frame(auc_list$X)
  dimnames(auc_df) <- list(auc_list$var$`_index`,
                           gsub("\\___cisTopic", "", auc_list$obs$Cell))
  return(auc_df)
}
#####################################################
library(cols4all)
library(rhdf5)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(glue)
library(stringr)
library(bigreadr)
##
proj_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(proj_path)
##
ann_col <- "ann_level4_final"
ct = "B"
meta_path <- "03_output/03_clustering/B/WNN_ADT_RNA/sc_meta.rds"

## Load the .h5mu file
sp_outpath <- glue("03_output/07_SCENIC/{ct}/scplus_pipeline/Snakemake/")
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds")) 
h5mu_file <- glue("{sp_outpath}scplusmdata.h5mu")
h5mu <- h5ls(h5mu_file)
#
gene_based_AUC <- rbind(load.AUC(h5mu_file, "/mod/direct_gene_based_AUC"),
                        load.AUC(h5mu_file, "/mod/extended_gene_based_AUC")) %>%
  as.data.frame()
region_based_AUC <- rbind(load.AUC(h5mu_file, "/mod/direct_region_based_AUC"),
                          load.AUC(h5mu_file, "/mod/extended_region_based_AUC")) %>%
  as.data.frame()
##
regulon_metadatata <- rbind(fread2(glue("{sp_outpath}eRegulon_direct.tsv")),
                            fread2(glue("{sp_outpath}eRegulons_extended.tsv")))
##
AUC_mat <- rbind(gene_based_AUC, region_based_AUC)[e_regulon_name_trans$signature,]
rownames(AUC_mat) <- e_regulon_name_trans$signature_name_simplify
sc_meta <- readRDS(meta_path)
sc_meta <- sc_meta[colnames(AUC_mat),]

####
scAUC_obj <- CreateSeuratObject(counts = as.data.frame(AUC_mat),
                                assay = "AUC",
                                meta.data = sc_meta)
# fake normalize and scale (FindMarkers will auto-detect)
scAUC_obj <- NormalizeData(scAUC_obj)
scAUC_obj[["AUC"]]$data <- scAUC_obj[["AUC"]]$counts
scAUC_obj[["AUC"]]$scale.data <- scAUC_obj[["AUC"]]$counts
Idents(scAUC_obj) <- scAUC_obj@meta.data$Condition

saveRDS(scAUC_obj, file = glue("{sp_outpath}scAUC_obj.rds"))
saveRDS(AUC_mat, file = glue("{sp_outpath}AUC_mat.rds"))

