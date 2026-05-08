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
plt_path = "04_plot/Fig2/"

##### load data #####
sc_meta_CD4 <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta.rds")
sc_meta_CD4$Condition <- factor(sc_meta_CD4$Condition, levels = c("NN", "NU", "II"))
sc_meta_CD4$Section <- factor(sc_meta_CD4$Section, levels = c("Colon", "TI"), 
                              labels = c("Colon", "Ileum"))
##### Supplementary Tables #####
##### Table S4 #####
use_de <- "DESeq2"
ann_col <- "ann_level2_refine"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/all/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/all/{use_de}/")
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4T"]]
deg_df_it2 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[["CD4T"]]
deg_df_cr1 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4T"]]
deg_df_cr2 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[["CD4T"]]
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
deg_df_sig$Contrast <- gsub("condition_", "",
                            deg_df_sig$Contrast) %>%
  gsub("_", " ", .)
write.table(deg_df_sig, file = glue("{plt_path}deg_df_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S5 #####
use_de <- "DESeq2"
ann_col <- "ann_level2_refine"
out_path <- "03_output/04_Diff/DAR/"
dar_path_it <- glue("{out_path}TI/all/{use_de}/")
dar_path_cr <- glue("{out_path}Colon/all/{use_de}/")
dar_df_it1 <- readRDS(glue("{dar_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4T"]]
dar_df_it2 <- readRDS(glue("{dar_path_it}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[["CD4T"]]
dar_df_cr1 <- readRDS(glue("{dar_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4T"]]
dar_df_cr2 <- readRDS(glue("{dar_path_cr}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[["CD4T"]]
dar_df_it1$Section <- "Ileum"
dar_df_it2$Section <- "Ileum"
dar_df_cr1$Section <- "Colon"
dar_df_cr2$Section <- "Colon"
dar_df_it1 <- dar_df_it1[order(dar_df_it1$padj, -abs(dar_df_it1$log2FC)),]
dar_df_it2 <- dar_df_it2[order(dar_df_it2$padj, -abs(dar_df_it2$log2FC)),]
dar_df_cr1 <- dar_df_cr1[order(dar_df_cr1$padj, -abs(dar_df_cr1$log2FC)),]
dar_df_cr2 <- dar_df_cr2[order(dar_df_cr2$padj, -abs(dar_df_cr2$log2FC)),]
#
dar_df_sig <- rbind(dar_df_it1[which(dar_df_it1$padj < 0.05 & abs(dar_df_it1$log2FC) > 0.25),],
                    dar_df_it2[which(dar_df_it2$padj < 0.05 & abs(dar_df_it2$log2FC) > 0.25),],
                    dar_df_cr1[which(dar_df_cr1$padj < 0.05 & abs(dar_df_cr1$log2FC) > 0.25),],
                    dar_df_cr2[which(dar_df_cr2$padj < 0.05 & abs(dar_df_cr2$log2FC) > 0.25),])
dar_df_sig$Contrast <- gsub("condition_", "",
                            dar_df_sig$Contrast) %>%
  gsub("_", " ", .)
write.table(dar_df_sig, file = glue("{plt_path}dar_df_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S6 #####
motif_enrich_it1 <- readRDS("03_output/04_Diff/DAR/TI/all/DESeq2/motif_enrich_df_CD4T_II_vs_NN.rds")
motif_enrich_it1 <- subset(motif_enrich_it1, fold.enrichment > 1)
motif_enrich_it1$log2FC <- ifelse(motif_enrich_it1$dir == "Up", 
                                  log2(motif_enrich_it1$fold.enrichment), 
                                  -log2(motif_enrich_it1$fold.enrichment))
motif_enrich_it1 <- motif_enrich_it1[order(motif_enrich_it1$p.adjust),]
motif_enrich_it1 <- motif_enrich_it1[!duplicated(motif_enrich_it1$motif),]
motif_enrich_it1$contrast <- "II vs NN"
#
motif_enrich_it2 <- readRDS("03_output/04_Diff/DAR/TI/all/DESeq2/motif_enrich_df_CD4T_NU_vs_NN.rds")
motif_enrich_it2 <- subset(motif_enrich_it2, fold.enrichment > 1)
motif_enrich_it2$log2FC <- ifelse(motif_enrich_it2$dir == "Up", 
                                  log2(motif_enrich_it2$fold.enrichment), 
                                  -log2(motif_enrich_it2$fold.enrichment))
motif_enrich_it2 <- motif_enrich_it2[order(motif_enrich_it2$p.adjust),]
motif_enrich_it2 <- motif_enrich_it2[!duplicated(motif_enrich_it2$motif),]
motif_enrich_it2$contrast <- "NU vs NN"

##
motif_enrich_it <- rbind(motif_enrich_it1, motif_enrich_it2)
motif_enrich_sig <- subset(motif_enrich_it,
                           fold.enrichment > 1 & p.adjust < 0.05,
                           select = c("Term", "motif.name", "fold.enrichment", "pvalue", "p.adjust", "dir", "contrast"))
write.table(motif_enrich_sig, file = glue("{plt_path}motif_enrich_sig.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S7 #####
## load differential test results
diff_auc_ti <- rbind(readRDS("03_output/04_Diff/AUC/TI/CD4T/scAUC_pseudo_pair_ann_level2_final_II_vs_NN.rds"),
                     readRDS("03_output/04_Diff/AUC/TI/CD4T/scAUC_pseudo_pair_ann_level2_final_NU_vs_NN.rds"))
diff_auc_ti <- diff_auc_ti[!is.na(diff_auc_ti$padj),]
# selected Regulons
diff_auc_ti_sel <- diff_auc_ti[grepl("ROR|BACH1|BATF|MAF", diff_auc_ti$Term) &
                                 diff_auc_ti$padj < 0.05,]
write.table(diff_auc_ti_sel, file = glue("{plt_path}diff_sig_auc_CD4_level2_final.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S8 #####
sc_meta_CD4$group_comb  <- paste0(sc_meta_CD4$Condition, ":", sc_meta_CD4$Section)
sc_meta_CD4$group_comb <- factor(sc_meta_CD4$group_comb, 
                                 levels = c("NN:Ileum", "NU:Ileum", "II:Ileum",
                                            "NN:Colon", "NU:Colon", "II:Colon"))
tab_ct <- table(sc_meta_CD4$ann_level4_final, sc_meta_CD4$group_comb)
tab_ct_prop <- round(prop.table(tab_ct, 2)*100, 2) %>% as.data.frame()
tab_ct_comb <- as.data.frame(tab_ct)
tab_ct_comb$comb <- paste0(tab_ct_comb$Freq, " (", tab_ct_prop$Freq, ")")
tab_ct_comb$Var1 <- factor(tab_ct_comb$Var1,
                           levels = levels(sc_meta_CD4$ann_level4_final),
                           labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
tab_ct_rp <- reshape2::dcast(tab_ct_comb, Var1 ~ Var2, value.var = "comb")
write.table(tab_ct_rp, file = glue("{plt_path}tab_ct_rp.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

## Table S9 ##
Treg_levels <- c("CD4_Treg_naive", "CD4_IKZF2low_Treg", "CD4_Treg")
diff_auc <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/diff_auc_treg.rds")
diff_auc$cluster <- gsub("_", " ", diff_auc$cluster) %>%
  factor(., levels = gsub("_", " ", Treg_levels))
diff_auc_sig <- subset(diff_auc, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
diff_auc_sig <- diff_auc_sig[order(diff_auc_sig$cluster,
                                   diff_auc_sig$p_val),
                             c("gene", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
colnames(diff_auc_sig) <- c("Term", "log2FC", "pvalue", "padj", "cluster")
write.table(diff_auc_sig, file = glue("{plt_path}diff_auc_sig_treg.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

## Table S10 ##
ct_markers <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/ct_markers_treg.rds")
ct_markers$cluster <- gsub("_", " ", ct_markers$cluster) %>%
  factor(., levels = gsub("_", " ", Treg_levels))
ct_markers_sig <- subset(ct_markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
ct_markers_sig <- ct_markers_sig[order(ct_markers_sig$cluster,
                                       ct_markers_sig$p_val),
                                 c("gene", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
colnames(ct_markers_sig) <- c("Term", "log2FC", "pvalue", "padj", "cluster")
write.table(ct_markers_sig, file = glue("{plt_path}ct_markers_treg.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

## Table S11 ##
diff_gsva <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/diff_gsva_treg.rds")
diff_gsva$cluster <- gsub("_", " ", diff_gsva$cluster) %>%
  factor(., levels = gsub("_", " ", Treg_levels))
diff_gsva$gene <- gsub("\\-", " ", diff_gsva$gene) %>%
  gsub("KEGG", "KEGG:", .) %>%
  gsub("REACTOME", "REACTOME:", .) %>%
  gsub("GOBP", "GOBP:", .)
diff_gsva_sig <- subset(diff_gsva, p_val_adj < 0.05)
diff_gsva_sig <- diff_gsva_sig[order(diff_gsva_sig$cluster,
                                     diff_gsva_sig$p_val),
                               c("gene", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
colnames(diff_gsva_sig) <- c("Term", "log2FC", "pvalue", "padj", "cluster")
write.table(diff_gsva_sig, file = glue("{plt_path}diff_gsva_treg.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

##### Table S12 #####
use_de <- "DESeq2"
ann_col <- "ann_level4_final"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD4T/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD4T/{use_de}/")
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
                             levels = levels(sc_meta_CD4$ann_level4_final),
                             labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
deg_df_it2$cluster <- factor(deg_df_it2$cluster,
                             levels = levels(sc_meta_CD4$ann_level4_final),
                             labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
deg_df_cr1$cluster <- factor(deg_df_cr1$cluster,
                             levels = levels(sc_meta_CD4$ann_level4_final),
                             labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
deg_df_cr2$cluster <- factor(deg_df_cr2$cluster,
                             levels = levels(sc_meta_CD4$ann_level4_final),
                             labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
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
  gsub("_", " ", .)
write.table(deg_df_sig, file = glue("{plt_path}deg_df_sig_level4.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

## Table S13 ##
use_de <- "DESeq2"
ann_col <- "ann_level4_final"
out_path <- "03_output/04_Diff/DAR/"
dar_path_it <- glue("{out_path}TI/CD4T/{use_de}/")
dar_path_cr <- glue("{out_path}Colon/CD4T/{use_de}/")

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
                             levels = levels(sc_meta_CD4$ann_level4_final),
                             labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
dar_df_it2$cluster <- factor(dar_df_it2$cluster,
                             levels = levels(sc_meta_CD4$ann_level4_final),
                             labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
dar_df_cr1$cluster <- factor(dar_df_cr1$cluster,
                             levels = levels(sc_meta_CD4$ann_level4_final),
                             labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
dar_df_cr2$cluster <- factor(dar_df_cr2$cluster,
                             levels = levels(sc_meta_CD4$ann_level4_final),
                             labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
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
  gsub("_", " ", .)
write.table(dar_df_sig, file = glue("{plt_path}dar_df_sig_level4.txt"), 
            sep = "\t", col.names = T, row.names = F, quote = F)

