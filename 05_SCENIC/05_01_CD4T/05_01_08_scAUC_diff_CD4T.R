## module load r/4.5.0
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
source("code/FUNCTION/diff/pseudo_pair_diff.R")
##
ct = "CD4T"
sp_outpath <- glue("03_output/07_SCENIC/{ct}/scplus_pipeline/Snakemake/")
scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
scAUC_obj$ann_level2_final <- "CD4T"
#
all_group <- c("II", "NU", "NN")
comb_group <- combn(all_group, 2)
#
out_path_ti <- glue("03_output/04_Diff/AUC/TI/{ct}/")
out_path_cr <- glue("03_output/04_Diff/AUC/Colon/{ct}/")
##
lapply(c("ann_level2_final", "ann_level4_final"), function(ann_col){
  apply(comb_group[, 2:3], 2, function(comb_groupx){
    g1 <- comb_groupx[1] # case
    g2 <- comb_groupx[2] # control
    scAUC_obj_sub <- subset(scAUC_obj, Condition %in% c(g1, g2))
    scAUC_obj_sub$Condition <- factor(scAUC_obj_sub$Condition, 
                                      levels = c(g2, g1))
    out_prefix <- glue("scAUC_pseudo_pair_{ann_col}_{g1}_vs_{g2}.rds")
    ## terminal intestine
    diff_auc_ti <- pseudo.pair.diff(seurat_obj = subset(scAUC_obj_sub,
                                                        Section == "TI"),
                                    assay_use = "AUC",
                                    ident_col = ann_col,
                                    group_col = "Condition",
                                    sample_col = "Sample_ID_exp",
                                    block_col = "Sample_exp",
                                    test_use = "wilcox") %>% 
      Reduce("rbind", .) %>% as.data.frame()
    diff_auc_ti$contrast <- glue("{g1}-{g2}")
    saveRDS(diff_auc_ti, file = glue("{out_path_ti}/{out_prefix}"))
    ## colon
    diff_auc_cr <- pseudo.pair.diff(seurat_obj = subset(scAUC_obj_sub,
                                                        Section == "Colon"),
                                    assay_use = "AUC",
                                    ident_col = ann_col,
                                    group_col = "Condition",
                                    sample_col = "Sample_ID_exp",
                                    block_col = "Sample_exp",
                                    test_use = "wilcox") %>% 
      Reduce("rbind", .) %>% as.data.frame()
    diff_auc_cr$contrast <- glue("{g1}-{g2}")
    saveRDS(diff_auc_cr, file = glue("{out_path_cr}/{out_prefix}"))
    #
    return(out_prefix)
  })
  
})

