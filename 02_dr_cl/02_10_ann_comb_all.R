## module load r/4.5.0
library(dplyr)
library(stringr)
library(tibble)
library(bigreadr)
library(glue)

## set path
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

##### combine annotation #####
sc_meta_all <- readRDS("03_output/03_clustering/all/WNN_ADT_RNA/sc_meta_ann_level2_refine.rds")
sc_meta_CD4T <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta.rds")
sc_meta_CD8T <- readRDS("03_output/03_clustering/CD8T/WNN_ADT_RNA/sc_meta.rds")
sc_meta_B <- readRDS("03_output/03_clustering/B/WNN_ADT_RNA/sc_meta.rds")
##
sc_meta_all$ann_level3_final <- sc_meta_all$ann_level2_refine %>% as.character()
sc_meta_all[rownames(sc_meta_CD4T), "ann_level3_final"] <- sc_meta_CD4T$ann_level3_final %>% as.character()
sc_meta_all[rownames(sc_meta_CD8T), "ann_level3_final"] <- sc_meta_CD8T$ann_level3_final %>% as.character()
sc_meta_all[rownames(sc_meta_B), "ann_level3_final"] <- sc_meta_B$ann_level3_final %>% as.character()
sc_meta_all$ann_level3_final <- factor(sc_meta_all$ann_level3_final,
                                       levels = c(levels(sc_meta_CD4T$ann_level3_final),
                                                  levels(sc_meta_CD8T$ann_level3_final),
                                                  "gdT",
                                                  levels(sc_meta_B$ann_level3_final),
                                                  "NK", "ILC", "Macrophage", "cDC", "pDC"))
##
sc_meta_all$ann_level4_final <- sc_meta_all$ann_level3_final %>% as.character()
sc_meta_all[rownames(sc_meta_CD4T), "ann_level4_final"] <- sc_meta_CD4T$ann_level4_final %>% as.character()
sc_meta_all[rownames(sc_meta_CD8T), "ann_level4_final"] <- sc_meta_CD8T$ann_level4_final %>% as.character()
sc_meta_all[rownames(sc_meta_B), "ann_level4_final"] <- sc_meta_B$ann_level4_final %>% as.character()
sc_meta_all$ann_level4_final <- factor(sc_meta_all$ann_level4_final,
                                       levels = c(levels(sc_meta_CD4T$ann_level4_final),
                                                  levels(sc_meta_CD8T$ann_level4_final),
                                                  "gdT",
                                                  levels(sc_meta_B$ann_level4_final),
                                                  "NK", "ILC", "Macrophage", "cDC", "pDC"))
## 
sc_meta_cd103trm <- readRDS("03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/sc_meta.rds")
sc_meta_all$ann_level5_final <- sc_meta_all$ann_level4_final %>% as.character()
sc_meta_all[rownames(sc_meta_cd103trm), "ann_level5_final"] <- sc_meta_cd103trm$ann_level5_final %>% as.character()
sc_meta_all$ann_level5_final <- factor(sc_meta_all$ann_level5_final,
                                       levels = c(levels(sc_meta_CD4T$ann_level4_final)[!grepl("CD103", levels(sc_meta_CD4T$ann_level4_final))],
                                                  levels(sc_meta_cd103trm$ann_level5_final),
                                                  levels(sc_meta_CD8T$ann_level4_final),
                                                  "gdT",
                                                  levels(sc_meta_B$ann_level4_final),
                                                  "NK", "ILC", "Macrophage", "cDC", "pDC"))
##
saveRDS(sc_meta_all, file = "03_output/03_clustering/sc_meta_ann_comb_all1205.rds")
#
ct_levels_refine <- list("ann_level2_refine" = levels(sc_meta_all$ann_level2_refine),
                         "ann_level3_final" = levels(sc_meta_all$ann_level3_final), 
                         "ann_level4_final" = levels(sc_meta_all$ann_level4_final),
                         "ann_level5_final" = levels(sc_meta_all$ann_level5_final))
saveRDS(ct_levels_refine, "03_output/03_clustering/ct_levels_refine.rds")

##### pair sample set #####
sc_meta <- readRDS(glue("03_output/03_clustering/sc_meta_ann_comb_all1205.rds"))
cov_use <- c("Age", "Gender", 
             "use_Aminosalic", "use_Immunomodu", "use_MTX",
             "use_AntiTNF", "use_AntiIL23", "use_AntiIntegrin")
##
sc_meta_bulk <- sc_meta[!duplicated(sc_meta$Sample_ID_exp),]
pair_set <- list("Colon_II_NN" = c("Colon", "II",  "NN"),
                 "Colon_NU_NN" = c("Colon", "NU",  "NN"),
                 "TI_II_NN" = c("TI", "II",  "NN"),
                 "TI_NU_NN" = c("TI", "NU",  "NN"),
                 "TI_II_NU" = c("TI", "II",  "NU"),
                 "TI_Colon_NN" = c("TI", "Colon",  "NN"))
pair_sample_set <- lapply(seq_along(pair_set), function(xx){
  
  pair_setx <- pair_set[[xx]]
  sc_meta_bulkx <- subset(sc_meta_bulk, Section %in% pair_setx &
                            Condition %in% pair_setx)
  paired_samplex <- sc_meta_bulkx$Sample_exp[duplicated(sc_meta_bulkx$Sample_exp)]
  sc_meta_pairx <- subset(sc_meta_bulkx, 
                          Sample_exp %in% paired_samplex,
                          select = c("Sample", "Sample_exp", "Sample_ID_exp", 
                                     "Condition", "Section",
                                     cov_use))
  if (grepl("NU|II", names(pair_set)[xx])) {
    sc_meta_pairx$group_pair <- sc_meta_pairx$Condition
  } else {
    sc_meta_pairx$group_pair <- sc_meta_pairx$Section
  }
  sc_meta_pairx$group_pair <- paste0(sc_meta_pairx$group_pair,
                                     "(",
                                     paste(pair_setx, collapse = "_"),
                                     ")")
  return(sc_meta_pairx)
  
}) %>% Reduce("rbind", .) %>% as.data.frame()
saveRDS(pair_sample_set, file = "03_output/03_clustering/pair_sample_set_1205.rds")

