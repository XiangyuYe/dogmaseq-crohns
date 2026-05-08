##
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(reshape2)
library(cols4all)
library(glue)
library(tibble)
library(ggpubr)
library(rstatix)

###data input and select parameters
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
# set path
data_path = "03_output/03_clustering/"
plt_path = "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/04_plot/Fig1/"

##### Supplementary Tables #####
##### Table S2 #####
# load data
sc_meta_all <- readRDS("03_output/03_clustering/sc_meta_ann_comb_all1205.rds")
sc_meta_all$Condition <- factor(sc_meta_all$Condition, levels = c("NN", "NU", "II"))
sc_meta_all$Section <- factor(sc_meta_all$Section, levels = c("Colon", "TI"))
sc_meta_all$group_comb  <- paste0(sc_meta_all$Condition, ":", sc_meta_all$Section)
sc_meta_all$group_comb <- factor(sc_meta_all$group_comb, 
                                 levels = c("NN:TI", "NU:TI", "II:TI",
                                            "NN:Colon", "NU:Colon", "II:Colon"))
tab_ct <- table(sc_meta_all$ann_level2_refine, sc_meta_all$group_comb)
tab_ct_prop <- round(prop.table(tab_ct, 2)*100, 2) %>% as.data.frame()
tab_ct_comb <- as.data.frame(tab_ct)
tab_ct_comb$comb <- paste0(tab_ct_comb$Freq, " (", tab_ct_prop$Freq, ")")
tab_ct_rp <- reshape2::dcast(tab_ct_comb, Var1 ~ Var2, value.var = "comb")
write.table(tab_ct_rp, file = glue("{plt_path}tab_ct_rp.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S3 #####
magma_df <- fread2("03_output/01_scDRS/01_magma/CD.genes.out")
magma_zdf <- fread2("03_output/01_scDRS/01_magma/zscore_CD.tsv")
magma_df$GENE <- magma_zdf$GENE
# GENE CHR  START   STOP NSNPS NPARAM     N    ZSTAT       P
magma_df <- magma_df[order(magma_df$ZSTAT, decreasing = T), 
                     c("GENE", "NSNPS", "N", "ZSTAT", "P")]

write.table(magma_df[1:1000,], file = glue("{plt_path}magma_df_top1000.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)


