## module load r/4.5.0
library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(reshape2)
library(glue)
library(bigreadr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggrepel)

###data input and select parameters
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
plt_path <- "04_plot/Fig4/"
data_path <- "03_output/03_clustering/CD4_CD103_TRM_Th17/WNN_RNA_ATAC/"
sc_meta_cd4cd103trm17 <- readRDS(glue("{data_path}sc_meta.rds"))

##### Supplementary Tables #####
##### Table S14 #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD4T/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD4T/{use_de}/")
ann_col <- "ann_level4_final"
ctx <- "CD4_CD103_TRM_Th17"
source("code/FUNCTION/diff/Enrichment.R")
enrich_bar_up <- enrich_barplt(prefix = glue("Up_IIvsNN_{ctx}"),
                               df_path = glue("{deg_path_it}Enrichment/"),
                               nn = 5,
                               min_gene = 10,
                               max_nchar = 50)
enrich_bar_down <- enrich_barplt(prefix = glue("Down_IIvsNN_{ctx}"),
                                 df_path = glue("{deg_path_it}Enrichment/"),
                                 nn = 5,
                                 min_gene = 10,
                                 max_nchar = 50)

enrich_df <- rbind(enrich_bar_up$enrich_df,
                   enrich_bar_down$enrich_df)
enrich_df$Direction <- c(rep("Up", nrow(enrich_bar_up$enrich_df)),
                         rep("Down", nrow(enrich_bar_down$enrich_df)))

enrich_dff <- subset(enrich_df, p.adjust < 0.05 & Count >= 10,
                     select = c("ID", "Group", "Description", "Count", "p.adjust", "Direction"))
write.table(enrich_dff, file = glue("{plt_path}enrich_dff_level4.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S15 #####
diff_auc_ti <- readRDS(glue("{diff_path_it}scAUC_pseudo_pair_{ann_col}_II_vs_NN.rds"))
diff_auc_ti <- diff_auc_ti[!is.na(diff_auc_ti$padj),]
diff_auc_ti_trm17 <- subset(diff_auc_ti, cluster == sel_ct)
diff_auc_ti_trm17_sig <- subset(diff_auc_ti_trm17, padj < 0.05&
                                  grepl("\\(\\+\\)", Term))
diff_auc_ti_trm17_sig$TF <- str_split_i(diff_auc_ti_trm17_sig$Term, "\\(", 1)
write.table(diff_auc_ti_trm17_sig, file = glue("{plt_path}diff_auc_ti_trm17_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S16 #####
ct_markers <- readRDS(glue("{data_path}ct_markers.rds"))
ct_markers$cluster <- gsub("_|\\-", " ", ct_markers$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
ct_markers_sig <- subset(ct_markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
ct_markers_sig <- ct_markers_sig[order(ct_markers_sig$cluster,
                                       ct_markers_sig$p_val),
                                 c("gene", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
colnames(ct_markers_sig) <- c("Term", "log2FC", "pvalue", "padj", "cluster")
write.table(ct_markers_sig, file = glue("{plt_path}ct_markers_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

## Table S17 ##
diff_auc <- readRDS(glue("{data_path}diff_auc.rds"))
diff_auc$cluster <- gsub("_|\\-", " ", diff_auc$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
diff_auc_sig <- subset(diff_auc, p_val_adj < 0.05)
diff_auc_sig <- diff_auc_sig[order(diff_auc_sig$cluster,
                                   diff_auc_sig$p_val),
                             c("gene", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
colnames(diff_auc_sig) <- c("Term", "log2FC", "pvalue", "padj", "cluster")
write.table(diff_auc_sig, file = glue("{plt_path}diff_auc_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S18 #####
system(glue("cp -f 03_output/03_clustering/CD4_CD103_TRM_Th17/cNMF/CD4_CD103_TRM_harmony_top_genes_k10.txt {plt_path}CD4_CD103_TRM_harmony_top_genes_k10.txt"))

##### Table S19 #####
sc_meta_cd4cd103trm17$group_comb  <- paste0(sc_meta_cd4cd103trm17$Condition, ":", sc_meta_cd4cd103trm17$Section)
sc_meta_cd4cd103trm17$group_comb <- factor(sc_meta_cd4cd103trm17$group_comb, 
                                           levels = c("NN:TI", "NU:TI", "II:TI",
                                                      "NN:Colon", "NU:Colon", "II:Colon"))
tab_ct <- table(sc_meta_cd4cd103trm17$ann_level5_final, sc_meta_cd4cd103trm17$group_comb)
tab_ct_prop <- round(prop.table(tab_ct, 2)*100, 2) %>% as.data.frame()
tab_ct_comb <- as.data.frame(tab_ct)
tab_ct_comb$comb <- paste0(tab_ct_comb$Freq, " (", tab_ct_prop$Freq, ")")
tab_ct_comb$Var1 <- gsub("_|\\-", " ", tab_ct_comb$Var1) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
tab_ct_rp <- reshape2::dcast(tab_ct_comb, Var1 ~ Var2, value.var = "comb")
write.table(tab_ct_rp, file = glue("{plt_path}tab_ct_rp.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S20 #####
use_de <- "DESeq2"
ann_col <- "ann_level5_final"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD4_CD103_TRM_Th17/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD4_CD103_TRM_Th17/{use_de}/")
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
deg_df_it1$cluster <- factor(deg_df_it1$cluster,
                             levels = levels(sc_meta_cd4cd103trm17$ann_level5_final),
                             labels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
deg_df_it2$cluster <- factor(deg_df_it2$cluster,
                             levels = levels(sc_meta_cd4cd103trm17$ann_level5_final),
                             labels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
deg_df_cr1$cluster <- factor(deg_df_cr1$cluster,
                             levels = levels(sc_meta_cd4cd103trm17$ann_level5_final),
                             labels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
deg_df_cr2$cluster <- factor(deg_df_cr2$cluster,
                             levels = levels(sc_meta_cd4cd103trm17$ann_level5_final),
                             labels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
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
write.table(deg_df_sig, file = glue("{plt_path}deg_df_sig_level5.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S21 #####
use_de <- "DESeq2"
ann_col <- "ann_level5_final"
out_path <- "03_output/04_Diff/DAR/"
dar_path_it <- glue("{out_path}TI/CD4_CD103_TRM_Th17/{use_de}/")
dar_path_cr <- glue("{out_path}Colon/CD4_CD103_TRM_Th17/{use_de}/")

dar_df_it1 <- readRDS(glue("{dar_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
dar_df_it2 <- readRDS(glue("{dar_path_it}/pseudo_bulk_{ann_col}_NU_vs_NN.rds")) %>%
  Reduce("rbind", .)
dar_df_cr1 <- readRDS(glue("{dar_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
dar_df_cr2 <- readRDS(glue("{dar_path_cr}/pseudo_bulk_{ann_col}_NU_vs_NN.rds")) %>%
  Reduce("rbind", .)
dar_df_it1$Section <- "Ileum"
dar_df_it2$Section <- "Ileum"
dar_df_cr1$Section <- "Colon"
dar_df_cr2$Section <- "Colon"
dar_df_it1$cluster <- factor(dar_df_it1$cluster,
                             levels = levels(sc_meta_cd4cd103trm17$ann_level5_final),
                             labels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
dar_df_it2$cluster <- factor(dar_df_it2$cluster,
                             levels = levels(sc_meta_cd4cd103trm17$ann_level5_final),
                             labels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
dar_df_cr1$cluster <- factor(dar_df_cr1$cluster,
                             levels = levels(sc_meta_cd4cd103trm17$ann_level5_final),
                             labels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
dar_df_cr2$cluster <- factor(dar_df_cr2$cluster,
                             levels = levels(sc_meta_cd4cd103trm17$ann_level5_final),
                             labels = gsub("_", " ", levels(sc_meta_cd4cd103trm17$ann_level5_final)))
dar_df_it1 <- dar_df_it1[order(dar_df_it1$cluster, dar_df_it1$padj, -abs(dar_df_it1$log2FC)),]
dar_df_it2 <- dar_df_it2[order(dar_df_it2$cluster, dar_df_it2$padj, -abs(dar_df_it2$log2FC)),]
dar_df_cr1 <- dar_df_cr1[order(dar_df_cr1$cluster, dar_df_cr1$padj, -abs(dar_df_cr1$log2FC)),]
dar_df_cr2 <- dar_df_cr2[order(dar_df_cr2$cluster, dar_df_cr2$padj, -abs(dar_df_cr2$log2FC)),]
#
dar_df_sig <- rbind(dar_df_it1[which(dar_df_it1$padj < 0.05 & abs(dar_df_it1$log2FC) > 0.25),],
                    dar_df_it2[which(dar_df_it2$padj < 0.05 & abs(dar_df_it2$log2FC) > 0.25),],
                    dar_df_cr1[which(dar_df_cr1$padj < 0.05 & abs(dar_df_cr1$log2FC) > 0.25),],
                    dar_df_cr2[which(dar_df_cr2$padj < 0.05 & abs(dar_df_cr2$log2FC) > 0.25),])
dar_df_sig$Contrast <- gsub("Condition_", "",
                            dar_df_sig$Contrast) %>%
  gsub("I_I", "II", .) %>%
  gsub("N_U", "NU", .) %>%
  gsub("N_N", "NN", .)  %>%
  gsub("_", " ", .)
write.table(dar_df_sig, file = glue("{plt_path}dar_df_sig_level5.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S22 #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD4_CD103_TRM_Th17/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD4_CD103_TRM_Th17/{use_de}/")
ann_col <- "ann_level5_final"
source("code/FUNCTION/diff/Enrichment.R")
#
ctx1 <- "CD4_CD103_TRM_Th17_HSP"
enrich_bar_down1 <- enrich_barplt(prefix = glue("Down_IIvsNN_{ctx1}"),
                                  df_path = glue("{deg_path_it}Enrichment/"),
                                  nn = 5,
                                  min_gene = 10,
                                  max_nchar = 50)
ctx2 <- "CD4_CD103_TRM_Th17_CREM"
enrich_bar_down2 <- enrich_barplt(prefix = glue("Down_IIvsNN_{ctx2}"),
                                  df_path = glue("{deg_path_it}Enrichment/"),
                                  nn = 5,
                                  min_gene = 10,
                                  max_nchar = 50)

enrich_df <- rbind(enrich_bar_down1$enrich_df,
                   enrich_bar_down2$enrich_df)
enrich_df$Direction <- c(rep("CD4_CD103_TRM_Th17_HSP", nrow(enrich_bar_down1$enrich_df)),
                         rep("CD4_CD103_TRM_Th17_CREM", nrow(enrich_bar_down2$enrich_df)))

enrich_dff <- subset(enrich_df, p.adjust < 0.05 & Count >= 10,
                     select = c("ID", "Group", "Description", "Count", "p.adjust", "Direction"))
write.table(enrich_dff, file = glue("{plt_path}enrich_dff_level5.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)



