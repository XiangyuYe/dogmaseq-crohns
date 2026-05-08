## module load r/4.5.0
library(bigreadr)
library(Seurat)
library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(glue)
##### data input and select parameters
Project_path="/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(Project_path)
##
ct <- "CD4_CD103_TRM_Th17"
cov_use = c("Age", "Gender", 
            "use_Aminosalic", "use_Immunomodu", "use_MTX",
            "use_AntiTNF", "use_AntiIL23", "use_AntiIntegrin")
block_col = "Sample_exp"
sample_col = "Sample_ID_exp"
ann_col <- "ann_level5_final"
group_col <- "Condition"
subgroup_col <- "Section"

##
data_path <- "03_output/03_clustering/"
out_path <- "03_output/04_Diff/DAR/"
sc_meta_cd4cd103trm17 <- readRDS(glue("{data_path}{ct}/WNN_RNA_ATAC/sc_meta.rds"))
sc_meta_cd4cd103trm17$Age <- scale(sc_meta_cd4cd103trm17$Age, center = T, scale = T) %>% as.numeric()
#
scATAC_obj <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_scATAC_obj.rds")
scATAC_obj <- subset(scATAC_obj, cells = rownames(sc_meta_cd4cd103trm17))
DefaultAssay(scATAC_obj)<-"peaks"
##
sc_meta <- sc_meta_cd4cd103trm17[, c(block_col, sample_col, ann_col, 
                         group_col, subgroup_col, cov_use)]
scATAC_obj@meta.data <- sc_meta
Idents(scATAC_obj) <- scATAC_obj@meta.data[[group_col]]
##
all_group <- c("II", "NU", "NN")
comb_group <- combn(all_group, 2)

source("code/FUNCTION/diff/pseudoBulk.R")
filter_ByCPM <- grepl("level2", ann_col)
out_path_ti <- glue("{out_path}TI/{ct}/")
out_path_cr <- glue("{out_path}Colon/{ct}/")

##### DESeq2 #####
de_method <- "DESeq2"
system(glue("mkdir -p {out_path_ti}/{de_method}/"))
system(glue("mkdir -p {out_path_cr}{de_method}/"))
for (ngx in c(2, 3, 1)) {
  
  g1 <- comb_group[1, ngx]
  g2 <- comb_group[2, ngx]
  if (g1 == "II" & g2 == "NU") {
    block_use <- NULL
  } else {
    block_use <- block_col
  }
  scATAC_obj_sub <- subset(scATAC_obj, Condition %in% c(g1, g2))
  scATAC_obj_sub$Condition <- factor(scATAC_obj_sub$Condition, 
                                     levels = c(g2, g1))
  out_prefix <- glue("pseudo_bulk_{ann_col}_{g1}_vs_{g2}.rds")
  ## terminal intestine
  dar_list_ti <- pseudo.bulk.deg(seurat_obj = subset(scATAC_obj_sub,
                                                     Section == "TI"),
                                 use_assay = "peaks",
                                 ident_col = ann_col,
                                 group_col = group_col,
                                 cov_col = cov_use,
                                 agg_strategy = "sum",
                                 vif_thre = 3,
                                 sample_col = sample_col,
                                 block_col = block_use,
                                 paired_only = T,
                                 do_norm = F,
                                 min_cell_per_sample = 10,
                                 min_percent = 0.3,
                                 filter_ByExpr = filter_ByCPM,
                                 test_use = "DESeq2-LRT",
                                 n_core = 10)
  saveRDS(dar_list_ti, file = glue("{out_path_ti}/{de_method}/{out_prefix}"))
  ## colon
  dar_list_cr <- pseudo.bulk.deg(seurat_obj = subset(scATAC_obj_sub,
                                                     Section == "Colon"),
                                 use_assay = "peaks",
                                 ident_col = ann_col,
                                 group_col = group_col,
                                 cov_col = cov_use,
                                 agg_strategy = "sum",
                                 vif_thre = 3,
                                 sample_col = sample_col,
                                 block_col = block_use,
                                 paired_only = T,
                                 do_norm = F,
                                 min_cell_per_sample = 10,
                                 min_percent = 0.3,
                                 filter_ByExpr = filter_ByCPM,
                                 test_use = "DESeq2-LRT",
                                 n_core = 10)
  saveRDS(dar_list_cr, file = glue("{out_path_cr}/{de_method}/{out_prefix}"))
  
}
