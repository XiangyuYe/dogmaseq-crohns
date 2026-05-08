## module load r/4.5.0
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(reshape2)
library(cols4all)
library(tibble)
library(ggpubr)
library(rstatix)

###data input and select parameters
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
# set path
plt_path = "04_plot/Fig6/"
data_path = "03_output/03_clustering/B/WNN_ADT_RNA/"
sc_meta_b <- readRDS("03_output/03_clustering/B/WNN_ADT_RNA/sc_meta.rds")

##### Table S31 ######s
ct_markers <- readRDS(glue("{data_path}ct_markers.rds"))
ct_markers$cluster <- gsub("_", " ", ct_markers$cluster) %>%
  factor(.,
         levels = c("B", "Plasma",
                    gsub("_", " ", levels(sc_meta_b$ann_level4_final))))
ct_markers_sig <- subset(ct_markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
ct_markers_sig <- ct_markers_sig[order(ct_markers_sig$cluster,
                                       ct_markers_sig$p_val),
                                 c("gene", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
colnames(ct_markers_sig) <- c("Term", "log2FC", "pvalue", "padj", "cluster")
write.table(ct_markers_sig, file = glue("{plt_path}ct_markers_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S32 ######s
use_de <- "DESeq2"
ann_col <- "ann_level4_final"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/B/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/B/{use_de}/")
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
deg_df_it2 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_NU_vs_NN.rds")) %>%
  Reduce("rbind", .)
deg_df_cr1 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
deg_df_cr2 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_NU_vs_NN.rds")) %>%
  Reduce("rbind", .)
deg_df_it1$Section <- "Ileum"
deg_df_it2$Section <- "Ileum"
deg_df_cr1$Section <- "Colon"
deg_df_cr2$Section <- "Colon"

deg_df_it1$cluster <- gsub("_", " ", deg_df_it1$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_b$ann_level4_final)))

deg_df_it2$cluster <- gsub("_", " ", deg_df_it2$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_b$ann_level4_final)))

deg_df_cr1$cluster <- gsub("_", " ", deg_df_cr1$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_b$ann_level4_final)))

deg_df_cr2$cluster <- gsub("_", " ", deg_df_cr2$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_b$ann_level4_final)))
deg_df_it1 <- deg_df_it1[order(deg_df_it1$cluster, deg_df_it1$padj, -abs(deg_df_it1$log2FC)),]
deg_df_it2 <- deg_df_it2[order(deg_df_it2$cluster, deg_df_it2$padj, -abs(deg_df_it2$log2FC)),]
deg_df_cr1 <- deg_df_cr1[order(deg_df_cr1$cluster, deg_df_cr1$padj, -abs(deg_df_cr1$log2FC)),]
deg_df_cr2 <- deg_df_cr2[order(deg_df_cr2$cluster, deg_df_cr2$padj, -abs(deg_df_cr2$log2FC)),]
#
deg_df_sig <- rbind(deg_df_it1[which(deg_df_it1$padj < 0.05 & abs(deg_df_it1$log2FC) > 0.25),],
                    deg_df_it2[which(deg_df_it2$padj < 0.05 & abs(deg_df_it2$log2FC) > 0.25),],
                    deg_df_cr1[which(deg_df_cr1$padj < 0.05 & abs(deg_df_cr1$log2FC) > 0.25),],
                    deg_df_cr2[which(deg_df_cr2$padj < 0.05 & abs(deg_df_cr2$log2FC) > 0.25),])
deg_df_sig$Contrast <- gsub("Condition_", "",
                            deg_df_sig$Contrast) %>%
  gsub("I_I", "II", .) %>%
  gsub("N_U", "NU", .) %>%
  gsub("N_N", "NN", .)  %>%
  gsub("_", " ", .)
write.table(deg_df_sig, file = glue("{plt_path}deg_df_sig_level4.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

