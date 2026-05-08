## model load r/4.5.0
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glue)
library(stringr)
library(reshape2)
library(cols4all)

###data input and select parameters
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
data_path <- "03_output/03_clustering/CD8T/WNN_ADT_RNA/"
plt_path <- "04_plot/Fig5/"
#
pair_list <- list("II_vs_NN_TI" = c("II:TI", "NN:TI"),
                  "NU_vs_NN_TI" = c("NU:TI", "NN:TI"),
                  "II_vs_NN_Colon" = c("II:Colon", "NN:Colon"),
                  "NU_vs_NN_Colon" = c("NU:Colon", "NN:Colon"),
                  "TI_vs_Colon_NN" = c("NN:TI", "NN:Colon"))
##
sc_meta_cd8 <- readRDS(glue("{data_path}sc_meta.rds"))
##### Supplementary Tables #####
##### Table S23 #####
use_de <- "DESeq2"
ann_col <- "ann_level2_refine"
out_path <- "03_output/04_Diff/DEG/"
ctx <- "CD8T"

deg_path_it <- glue("{out_path}TI/all/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/all/{use_de}/")
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
deg_df_it2 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[[ctx]]
deg_df_cr1 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
deg_df_cr2 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[[ctx]]
deg_df_it1$Section <- "Ileum"
deg_df_it2$Section <- "Ileum"
deg_df_cr1$Section <- "Colon"
deg_df_cr2$Section <- "Colon"
deg_df_it1 <- deg_df_it1[order(deg_df_it1$padj, -abs(deg_df_it1$log2FC)),]
deg_df_it2 <- deg_df_it2[order(deg_df_it2$padj, -abs(deg_df_it2$log2FC)),]
deg_df_cr1 <- deg_df_cr1[order(deg_df_cr1$padj, -abs(deg_df_cr1$log2FC)),]
deg_df_cr2 <- deg_df_cr2[order(deg_df_cr2$padj, -abs(deg_df_cr2$log2FC)),]
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
write.table(deg_df_sig, file = glue("{plt_path}deg_df_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S24 #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/all/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/all/{use_de}/")
ctx <- "CD8T"

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
write.table(enrich_dff, file = glue("{plt_path}enrich_dff_level2.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S25 #####
ct_markers <- readRDS(glue("{data_path}ct_markers.rds"))
ct_markers$cluster <- gsub("_", " ", ct_markers$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final)))
ct_markers_sig <- subset(ct_markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
ct_markers_sig <- ct_markers_sig[order(ct_markers_sig$cluster,
                                       ct_markers_sig$p_val),
                                 c("gene", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
colnames(ct_markers_sig) <- c("Term", "log2FC", "pvalue", "padj", "cluster")
write.table(ct_markers_sig, file = glue("{plt_path}ct_markers_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S26 #####
diff_auc <- readRDS(glue("{data_path}diff_auc.rds"))
diff_auc$cluster <- gsub("_", " ", diff_auc$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final)))
diff_auc_sig <- subset(diff_auc, p_val_adj < 0.05)
diff_auc_sig <- diff_auc_sig[order(diff_auc_sig$cluster,
                                   diff_auc_sig$p_val),
                             c("gene", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
colnames(diff_auc_sig) <- c("Term", "log2FC", "pvalue", "padj", "cluster")
write.table(diff_auc_sig, file = glue("{plt_path}diff_auc_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S27 #####
system(glue("cp -f 03_output/03_clustering/CD8_TRM/cNMF/CD8_TRM_harmony_top_genes_k9.txt {plt_path}CD8_TRM_harmony_top_genes_k9.txt"))

##### Table S28 #####
sc_meta_cd8$group_comb  <- paste0(sc_meta_cd8$Condition, ":", sc_meta_cd8$Section)
sc_meta_cd8$group_comb <- factor(sc_meta_cd8$group_comb, 
                                 levels = c("NN:TI", "NU:TI", "II:TI",
                                            "NN:Colon", "NU:Colon", "II:Colon"))
tab_ct <- table(sc_meta_cd8$ann_level4_final, sc_meta_cd8$group_comb)
tab_ct_prop <- round(prop.table(tab_ct, 2)*100, 2) %>% as.data.frame()
tab_ct_comb <- as.data.frame(tab_ct)
tab_ct_comb$comb <- paste0(tab_ct_comb$Freq, " (", tab_ct_prop$Freq, ")")
tab_ct_comb$Var1 <- factor(tab_ct_comb$Var1,
                           levels = levels(sc_meta_cd8$ann_level4_final),
                           labels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final)))
tab_ct_rp <- reshape2::dcast(tab_ct_comb, Var1 ~ Var2, value.var = "comb")
write.table(tab_ct_rp, file = glue("{plt_path}tab_ct_rp.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S29 #####
use_de <- "DESeq2"
ann_col <- "ann_level4_final"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD8T/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD8T/{use_de}/")
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
         levels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final)))

deg_df_it2$cluster <- gsub("_", " ", deg_df_it2$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final)))

deg_df_cr1$cluster <- gsub("_", " ", deg_df_cr1$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final)))

deg_df_cr2$cluster <- gsub("_", " ", deg_df_cr2$cluster) %>%
  factor(.,
         levels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final)))
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

##### Table S30 #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD8T/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD8T/{use_de}/")
ann_col <- "ann_level4_final"
ctx <- "CD8_TRM_Tc17_IL26"
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

enrich_dff <- subset(enrich_df, p.adjust < 0.05 & Count >= 5,
                     select = c("ID", "Group", "Description", "Count", "p.adjust", "Direction"))
write.table(enrich_dff, file = glue("{plt_path}enrich_dff_level4.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)
