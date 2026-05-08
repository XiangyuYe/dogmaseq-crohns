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
ct <- "CD4T"
sp_outpath <- glue("03_output/07_SCENIC/{ct}/scplus_pipeline/Snakemake/")
h5mu_file <- glue("{sp_outpath}scplusmdata.h5mu")
h5mu <- h5ls(h5mu_file)

gene_based_AUC <- rbind(load.AUC(h5mu_file, "/mod/direct_gene_based_AUC"),
                        load.AUC(h5mu_file, "/mod/extended_gene_based_AUC")) %>% t
region_based_AUC <- rbind(load.AUC(h5mu_file, "/mod/direct_region_based_AUC"),
                          load.AUC(h5mu_file, "/mod/extended_region_based_AUC")) %>% t

##
regulon_metadatata <- rbind(fread2(glue("{sp_outpath}eRegulon_direct.tsv")),
                            fread2(glue("{sp_outpath}eRegulons_extended.tsv")))
##
colnames(gene_based_AUC) <- str_split_i(colnames(gene_based_AUC), "\\(", 1) %>%
  gsub("\\_$", "", .)
colnames(region_based_AUC) <- str_split_i(colnames(region_based_AUC), "\\(", 1) %>%
  gsub("\\_$", "", .)
region_based_AUC <- region_based_AUC[,colnames(gene_based_AUC)]

## calculate correlation
Rho_gr_eregulon <- lapply(1:ncol(gene_based_AUC), function(nn){
  
  c(cor.test(gene_based_AUC[,nn], region_based_AUC[,nn])$estimate,
    cor.test(gene_based_AUC[,nn], region_based_AUC[,nn], method = "spearman")$estimate)
}) %>% Reduce("rbind", .) %>% as.data.frame()
dimnames(Rho_gr_eregulon) <- list(colnames(gene_based_AUC),
                                  c("Rho_pearson", "Rho_spearman"))
regulon_metadatata$Rho_pearson <- Rho_gr_eregulon$Rho_pearson[match(regulon_metadatata$eRegulon_name,
                                                                    colnames(gene_based_AUC))]
regulon_metadatata$Rho_spearman <- Rho_gr_eregulon$Rho_spearman[match(regulon_metadatata$eRegulon_name,
                                                                      colnames(gene_based_AUC))]
saveRDS(regulon_metadatata, file = glue("{sp_outpath}regulon_metadatata_rho_add.rds")) 

## add correlation to filtered metadata
e_regulon_filter <- fread2(glue("{sp_outpath}e_regulon_metadata_filtered.txt"))
regulon_metadatata <- readRDS(glue("{sp_outpath}regulon_metadatata_rho_add.rds")) 
e_regulon_filter$Rho_pearson <- Rho_gr_eregulon$Rho_pearson[match(e_regulon_filter$eRegulon_name,
                                                                  colnames(gene_based_AUC))]
e_regulon_filter$Rho_spearman <- Rho_gr_eregulon$Rho_spearman[match(e_regulon_filter$eRegulon_name,
                                                                    colnames(gene_based_AUC))]
## add simplify name
e_regulon_filter$Region_signature_name_simplify <- 
  e_regulon_filter$Region_signature_name %>%
  gsub("_extended_|_direct_", "", .) %>%
  gsub("\\+\\/\\+", "\\(\\+\\)", .) %>%
  gsub("\\-\\/\\+", "\\(\\-\\)", .)
e_regulon_filter$Gene_signature_name_simplify <- 
  e_regulon_filter$Gene_signature_name %>%
  gsub("_extended_|_direct_", "", .) %>%
  gsub("\\+\\/\\+", "\\(\\+\\)", .) %>%
  gsub("\\-\\/\\+", "\\(\\-\\)", .)
## add name for match
e_regulon_filter$Region_signature_mkname <- 
  e_regulon_filter$Region_signature_name %>%
  gsub("\\_", "\\-", .)
e_regulon_filter$Gene_signature_mkname <- 
  e_regulon_filter$Gene_signature_name %>%
  gsub("\\_", "\\-", .)
saveRDS(e_regulon_filter, file = glue("{sp_outpath}e_regulon_filter_rho_add.rds"))

## additional filteration by RNA-ATAC correlation
e_regulon_ff <- subset(e_regulon_filter, e_regulon_filter$Rho_spearman > 0.2)
fwrite2(e_regulon_ff, file = glue("{sp_outpath}e_regulon_metadata_ff.txt"))
## e_regulon_name_trans for match in subsequent analysis
e_regulon_ff_unique <- e_regulon_ff[!duplicated(e_regulon_ff$eRegulon_name),]
e_regulon_name_trans <- data.frame(signature = c(e_regulon_ff_unique$Gene_signature_name,
                                                 e_regulon_ff_unique$Region_signature_name),
                                   signature_mkname = c(e_regulon_ff_unique$Gene_signature_mkname,
                                                        e_regulon_ff_unique$Region_signature_mkname),
                                   signature_name_simplify = c(e_regulon_ff_unique$Gene_signature_name_simplify,
                                                               e_regulon_ff_unique$Region_signature_name_simplify))
fwrite2(e_regulon_name_trans, file = glue("{sp_outpath}e_regulon_name_trans.txt"))
saveRDS(e_regulon_name_trans, file = glue("{sp_outpath}e_regulon_name_trans.rds"))
## direct regulon only
e_regulon_name_trans_dir_pos <- e_regulon_name_trans[grepl("\\+", e_regulon_name_trans$signature_name_simplify) &
                                                       grepl("direct", e_regulon_name_trans$signature),]
fwrite2(e_regulon_name_trans_dir_pos, file = glue("{sp_outpath}e_regulon_name_trans_dir_pos.txt"))
##
e_regulon_ff_pos <- subset(e_regulon_ff, !e_regulon_ff$is_extended & grepl("\\+", e_regulon_ff$Gene_signature_name_simplify))
fwrite2(e_regulon_ff_pos, file = glue("{sp_outpath}e_regulon_ff_pos.txt"))

