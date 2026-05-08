##
library(bigreadr)
library(Seurat)
library(Seurat)
library(dplyr)
library(tibble)
library(stringr)
library(glue)
library(SCOPfunctions)
##### data input and select parameters
Project_path="/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(Project_path)
##
ct <- "CD4T"
cov_use = c("Age", "Gender", 
            "use_Aminosalic", "use_Immunomodu", "use_MTX",
            "use_AntiTNF", "use_AntiIL23", "use_AntiIntegrin")
block_col = "Sample_exp"
sample_col = "Sample_ID_exp"
ann_col <- "ann_level4_final"
group_col <- "Condition"
subgroup_col <- "Section"

##
data_path <- "03_output/03_clustering/"
out_path <- "03_output/04_Diff/DEG/"
sc_meta_cd4 <- readRDS(glue("{data_path}{ct}/WNN_ADT_RNA/sc_meta.rds"))
sc_meta_cd4$Age <- scale(sc_meta_cd4$Age, center = T, scale = T) %>% as.numeric()
#
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
scRNA_obj <- subset(scRNA_obj, cells = rownames(sc_meta_cd4))
DefaultAssay(scRNA_obj)<-"RNA"
##
sc_meta <- sc_meta_cd4[, c(block_col, sample_col, ann_col, 
                           group_col, subgroup_col, cov_use)]
scRNA_obj@meta.data <- sc_meta
Idents(scRNA_obj) <- scRNA_obj@meta.data[[group_col]]
##
all_group <- c("II", "NU", "NN")
comb_group <- combn(all_group, 2)

source("code/FUNCTION/diff/pseudoBulk.R")
filter_ByCPM <- grepl("level2", ann_col)
out_path_ti <- glue("{out_path}TI/{ct}/")
out_path_rc <- glue("{out_path}Colon/{ct}/")

##### DESeq2 #####
de_method <- "DESeq2"
system(glue("mkdir -p {out_path_ti}/{de_method}/"))
system(glue("mkdir -p {out_path_rc}{de_method}/"))
for (ngx in c(2, 3, 1)) {
  
  g1 <- comb_group[1, ngx]
  g2 <- comb_group[2, ngx]
  if (g1 == "II" & g2 == "NU") {
    block_use <- NULL
  } else {
    block_use <- block_col
  }
  scRNA_obj_sub <- subset(scRNA_obj, Condition %in% c(g1, g2))
  scRNA_obj_sub$Condition <- factor(scRNA_obj_sub$Condition, 
                                     levels = c(g2, g1))
  out_prefix <- glue("pseudo_bulk_{ann_col}_{g1}_vs_{g2}.rds")
  ## small intestine
  deg_list_ti1 <- pseudo.bulk.deg(seurat_obj = subset(scRNA_obj_sub,
                                                     Section == "TI" & Gender == 1),
                                 use_assay = "RNA",
                                 ident_col = ann_col,
                                 group_col = group_col,
                                 cov_col = setdiff(cov_use, "Gender"),
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
  saveRDS(deg_list_ti1, file = glue("{out_path_ti}/{de_method}/{out_prefix}"))
  ## small intestine
  deg_list_ti2 <- pseudo.bulk.deg(seurat_obj = subset(scRNA_obj_sub,
                                                     Section == "TI" & Gender == 2),
                                 use_assay = "RNA",
                                 ident_col = ann_col,
                                 group_col = group_col,
                                 cov_col = setdiff(cov_use, "Gender"),
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
  saveRDS(deg_list_ti2, file = glue("{out_path_ti}/{de_method}/{out_prefix}"))
  ## colorect
  deg_list_rc1 <- pseudo.bulk.deg(seurat_obj = subset(scRNA_obj_sub,
                                                     Section == "Colon" & Gender == 1),
                                 use_assay = "RNA",
                                 ident_col = ann_col,
                                 group_col = group_col,
                                 cov_col = setdiff(cov_use, "Gender"),
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
  saveRDS(deg_list_rc1, file = glue("{out_path_rc}/{de_method}/{out_prefix}_F"))
  ## colorect
  deg_list_rc2 <- pseudo.bulk.deg(seurat_obj = subset(scRNA_obj_sub,
                                                     Section == "Colon" & Gender == 2),
                                 use_assay = "RNA",
                                 ident_col = ann_col,
                                 group_col = group_col,
                                 cov_col = setdiff(cov_use, "Gender"),
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
  saveRDS(deg_list_rc2, file = glue("{out_path_rc}/{de_method}/{out_prefix}_M"))
}
