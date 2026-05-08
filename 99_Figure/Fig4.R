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
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/"
setwd(project_path)

cd4cd103trm_cols <- c(
  "#FF7B36", "#FE0000",
  "#FDD501", "#DDCB82",
  "#4691CC", "#A6D6FA",
  "#082C8C", "#112E5A", 
  "#EC4B88"
)
pair_list <- list("II_vs_NN_TI" = c("II:TI", "NN:TI"),
                  "NU_vs_NN_TI" = c("NU:TI", "NN:TI"),
                  "II_vs_NN_Colon" = c("II:Colon", "NN:Colon"),
                  "NU_vs_NN_Colon" = c("NU:Colon", "NN:Colon"),
                  "TI_vs_Colon_NN" = c("NN:TI", "NN:Colon"))

source("code/FUNCTION/process/PROCESS_FUN.R")
plt_path <- "04_plot/Fig4/"
seed_use <- 20250528

##### load data #####
data_path <- "03_output/03_clustering/CD4_CD103_TRM_Th17/WNN_RNA_ATAC/"
sc_obj <- readRDS(glue("{data_path}scWNN_obj.rds"))
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)
sc_meta_cd4cd103trm17 <- readRDS(glue("{data_path}sc_meta.rds"))

sc_meta_cd4cd103trm <- readRDS(glue("03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/sc_meta.rds"))

#
sc_obj@meta.data <- sc_meta_cd4cd103trm17[colnames(sc_obj),]
sc_obj$Condition <- factor(sc_obj$Condition, levels = c("NN", "NU", "II"))
sc_obj$Section <- factor(sc_obj$Section, levels = c("Colon", "TI"))
Idents(sc_obj) <-  factor(sc_obj$ann_level5_final,
                          labels = gsub("_", " ", levels(sc_obj$ann_level5_final)))

##### Fig. 4A: volcano plot #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD4T/{use_de}/")
ann_col <- "ann_level4_final"
sel_ct <- "CD4_CD103_TRM_Th17"
source("code/FUNCTION/diff/volca_plot.R")
## 1. Ileum
deg_df_it_cd4x <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[sel_ct]]
deg_df_it_cd4x <- subset(deg_df_it_cd4x, !is.na(deg_df_it_cd4x$padj))
vc_plt_cd4cd103trm17 <- volca_function(deg_df = deg_df_it_cd4x,
                                       term_col = "Term",
                                       fc_col = "log2FC",
                                       p_col = "padj",
                                       n_top = 30,
                                       log_thresh = 0.25, 
                                       top_log_thresh = 0.5,
                                       group_label = c("Down regulation", "Others", "Up regulation"),
                                       valco_col = c4a("classic_blue_red12", 3),
                                       highlight_col = "orange") + 
  ggtitle("CD4 CD103 TRM-Th17") + 
  xlab(bquote(~Log[2]~"(Fold Change) in Ileal II vs NN")) +
  theme(plot.title = element_text(size = 15, face = "bold"))

ggsave(file = glue("{plt_path}/valcano_CD103TRMTh17_II_vs_NN_TI.png"),
       vc_plt_cd4cd103trm17,
       width = 10, height = 6,units = "in", dpi = 300, limitsize = T)
##
deg_path_it <- glue("03_output/04_Diff/DEG/TI/CD4T/DESeq2/")
source("code/FUNCTION/diff/Enrichment.R")
up_gene <- deg_df_it_cd4x$Term[which(deg_df_it_cd4x$padj < 0.05 & deg_df_it_cd4x$log2FC > 0.25)]
down_gene <- deg_df_it_cd4x$Term[which(deg_df_it_cd4x$padj < 0.05 & deg_df_it_cd4x$log2FC < -0.25)]

Enrichment.pipline(gene_list = up_gene,
                   list_name = glue("Up_IIvsNN_{sel_ct}"),
                   num_show = 10,
                   plot = T,
                   GO = T,
                   KEGG = T,
                   Reactome = T,
                   outpath = glue("{deg_path_it}Enrichment"))
Enrichment.pipline(gene_list = down_gene,
                   list_name = glue("Down_IIvsNN_{sel_ct}"),
                   num_show = 10,
                   plot = T,
                   GO = T,
                   KEGG = T,
                   Reactome = T,
                   outpath = glue("{deg_path_it}Enrichment"))

enrich_bar_up <- enrich_barplt(prefix = glue("Up_IIvsNN_{sel_ct}"),
                               df_path = glue("{deg_path_it}Enrichment/"),
                               nn = 5, 
                               min_gene = 10,
                               max_nchar = 50)
enrich_bar_down <- enrich_barplt(prefix = glue("Down_IIvsNN_{sel_ct}"),
                                 df_path = glue("{deg_path_it}Enrichment/"),
                                 nn = 5, 
                                 min_gene = 10,
                                 max_nchar = 50)
ggsave(glue("{plt_path}Enrich_IIvsNN_{sel_ct}_up.png"),
       enrich_bar_up$bar_plt,
       height = 5, width = 7, units = "in", dpi = 300)
ggsave(glue("{plt_path}Enrich_IIvsNN_{sel_ct}_down.png"),
       enrich_bar_down$bar_plt,
       height = 5, width = 7, units = "in", dpi = 300)
##### Fig. 4B: Inflammation_score #####
library(GSVA)
library(GSEABase)
## 
use_de <- "DESeq2"
ann_col <- "ann_level4_final"
ctx <- "CD4_CD103_TRM_Th17"

deg_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{deg_path}TI/CD4T/{use_de}/")
deg_ii_nn_it <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
#
deg_up <- deg_ii_nn_it$Term[which(deg_ii_nn_it$padj < 0.05 &
                                    deg_ii_nn_it$log2FC > 0.25)]
deg_down <- deg_ii_nn_it$Term[which(deg_ii_nn_it$padj < 0.05 &
                                      deg_ii_nn_it$log2FC < -0.25)]
# ##
# deg_path_cr <- glue("{deg_path}Colon/CD4T/{use_de}/")
# deg_ii_nn_cr <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
# #
# deg_up_cr <- deg_ii_nn_cr$Term[which(deg_ii_nn_cr$padj < 0.05 &
#                                     deg_ii_nn_cr$log2FC > 0.25)]
# deg_down_cr <- deg_ii_nn_cr$Term[which(deg_ii_nn_cr$padj < 0.05 &
#                                       deg_ii_nn_cr$log2FC < -0.25)]
## pseudobulk expression data
sample_col <- "Sample_ID_exp"
bulk_meta <- sc_meta_cd4cd103trm17[!duplicated(sc_meta_cd4cd103trm17[[sample_col]]),]
bulk_meta$group_comb  <- paste0(bulk_meta$Condition, ":", bulk_meta$Section)
sample_list <- split(bulk_meta$Sample_exp, f = bulk_meta$group_comb)
rownames(bulk_meta) <- gsub("_", "-", bulk_meta$Sample_ID_exp)
#
clinic_df <- readRDS("03_output/02_clean/clinic_df.rds")
exp_bulk <- AggregateExpression(sc_obj, 
                                assays = "RNA", 
                                return.seurat = F, 
                                group.by = "Sample_ID_exp")[["RNA"]]
exp_bulk_lcpm <- edgeR::cpm(exp_bulk, 
                            normalized.lib.sizes = T,
                            log = T)
## GSVA 
up_set <- GeneSetCollection(GeneSet(deg_up, setName = "deg_up"))
ssgsva_param_up <- ssgseaParam(exprData = exp_bulk_lcpm,
                               geneSets = up_set,
                               normalize = F)
gsva_df_up <- gsva(ssgsva_param_up) %>%
  t %>% as.data.frame()
#
down_set <- GeneSetCollection(GeneSet(deg_down, setName = "deg_down"))
ssgsva_param_down <- ssgseaParam(exprData = exp_bulk_lcpm,
                                 geneSets = down_set,
                                 normalize = F)
gsva_df_down <- gsva(ssgsva_param_down) %>%
  t %>% as.data.frame()
#
bulk_meta$score_cd4cd103trm17_up <- gsva_df_up[rownames(bulk_meta), "deg_up"] %>% 
  scales::rescale(., c(0, 10))
bulk_meta$score_cd4cd103trm17_down <- gsva_df_down[rownames(bulk_meta), "deg_down"] %>% 
  scales::rescale(., c(0, 10))
# add NI_score
ssgsea_df  <- readRDS("ssgsva_NI_score.rds")
bulk_meta$if_score_scale <- ssgsea_df[rownames(bulk_meta), "NIScore"] %>% 
  scales::rescale(., c(0, 10))

## plot
pair_list <- list("II_vs_NN_TI" = c("II:TI", "NN:TI"),
                  "NU_vs_NN_TI" = c("NU:TI", "NN:TI"),
                  "II_vs_NN_Colon" = c("II:Colon", "NN:Colon"),
                  "NU_vs_NN_Colon" = c("NU:Colon", "NN:Colon"),
                  "TI_vs_Colon_NN" = c("NN:TI", "NN:Colon"))
bulk_meta$Section <- factor(bulk_meta$Section, 
                            levels = c("TI", "Colon"),
                            labels = c("Ileum", "Colon"))
bulk_meta$group <- factor(bulk_meta$Condition, 
                          levels = c("NN", "NU", "II"))
## score_cd4cd103trm17_up
bulk_meta$score_cd4cd103trm17 <- bulk_meta$score_cd4cd103trm17_up
ppair_df <- lapply(1:4, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  bulk_metax <- subset(bulk_meta, Sample_exp %in% sample_pair_use)
  bulk_metax$group <- droplevels(bulk_metax$group)
  group_level <- levels(bulk_metax$group)
  
  pair_dfx <- split(bulk_metax, f = bulk_metax$Sample_exp) %>%
    lapply(., function(x){
      propx <- x$score_cd4cd103trm17[match(group_level, x$group)]
    }) %>% Reduce("rbind", .)
  ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
  
  return(ppairx)
  
}) %>% unlist()

padjpair_df <- p.adjust(ppair_df, method = "BH")
names(ppair_df) <- names(padjpair_df) <- names(pair_list)[1:4]

##
padj <- format(padjpair_df, scientific = T, digits = 3) %>%
  gsub("e", "E", .) %>% as.vector()
diff_df <- data.frame(y.position = c(10, 9.8, 10, 9.8),
                      x = c(1.5, 2.5, 1.5, 2.5),
                      xmin = c(2, 2, 2, 2),
                      xmax = c(1, 3, 1, 3),
                      group2 = c("NN", "NN", "NN", "NN"),
                      group1 = c("II", "NU", "II", "NU"),
                      Section = c("Ileum", "Ileum", "Colon", "Colon"),
                      padj_sig = padj)
if_box_up <- ggplot(bulk_meta, aes(x = group, y = score_cd4cd103trm17, color = group)) +
  geom_boxplot(fill = NA, width = 0.5) +
  geom_point(aes(color = group), position = position_jitter(width = 0.1)) +
  geom_line(aes(group = Sample_exp), color = "gray", alpha = 0.5) + 
  scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
  xlim(c("II", "NN", "NU")) + 
  stat_pvalue_manual(diff_df,
                     label = "padj_sig",
                     tip.length = 0.01,
                     label.size = 4,
                     bracket.size = 1) +
  facet_wrap(~ Section, scales = 'free_y', nrow = 1) +
  xlab("") + ylab("Inflammation score") + 
  ggtitle("Up-regulated DEGs") + 
  theme_bw() + 
  theme(legend.position = "none",
        title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold"))
ggsave(glue("{plt_path}/if_box_cd4cd103trm17_up.png"),
       if_box_up,
       height = 6, width = 5, units = "in", dpi = 300)

## score_cd4cd103trm17_down
bulk_meta$score_cd4cd103trm17 <- bulk_meta$score_cd4cd103trm17_down
ppair_df <- lapply(1:4, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  bulk_metax <- subset(bulk_meta, Sample_exp %in% sample_pair_use)
  bulk_metax$group <- droplevels(bulk_metax$group)
  group_level <- levels(bulk_metax$group)
  
  pair_dfx <- split(bulk_metax, f = bulk_metax$Sample_exp) %>%
    lapply(., function(x){
      propx <- x$score_cd4cd103trm17[match(group_level, x$group)]
    }) %>% Reduce("rbind", .)
  ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
  
  return(ppairx)
  
}) %>% unlist()

padjpair_df <- p.adjust(ppair_df, method = "BH")
names(ppair_df) <- names(padjpair_df) <- names(pair_list)[1:4]

##
padj <- format(ppair_df, scientific = T, digits = 3) %>%
  gsub("e", "E", .) %>% as.vector()
diff_df <- data.frame(y.position = c(10, 9.8, 10, 9.8),
                      x = c(1.5, 2.5, 1.5, 2.5),
                      xmin = c(2, 2, 2, 2),
                      xmax = c(1, 3, 1, 3),
                      group2 = c("NN", "NN", "NN", "NN"),
                      group1 = c("II", "NU", "II", "NU"),
                      Section = c("Ileum", "Ileum", "Colon", "Colon"),
                      padj_sig = padj)
if_box_down <- ggplot(bulk_meta, aes(x = group, y = score_cd4cd103trm17, color = group)) +
  geom_boxplot(fill = NA, width = 0.5) +
  geom_point(aes(color = group), position = position_jitter(width = 0.1)) +
  geom_line(aes(group = Sample_exp), color = "gray", alpha = 0.5) + 
  scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
  xlim(c("II", "NN", "NU")) + 
  stat_pvalue_manual(diff_df,
                     label = "padj_sig",
                     tip.length = 0.01,
                     label.size = 4,
                     bracket.size = 1) +
  facet_wrap(~ Section, scales = 'free_y', nrow = 1) +
  xlab("") + ylab("Inflammation score") + 
  ggtitle("Down-regulated DEGs") + 
  theme_bw() + 
  theme(legend.position = "none",
        title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold"))
ggsave(glue("{plt_path}/if_box_cd4cd103trm17_down.png"),
       if_box_down,
       height = 6, width = 5, units = "in", dpi = 300)

##### Fig. 4C: Diff AUC #####
## load SCENIC+ object
sp_outpath <- glue("03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/")
scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
## load differential test results
diff_path_it <- "03_output/04_Diff/AUC/TI/CD4T/"
ann_col <- "ann_level4_final"
sel_ct <- "CD4_CD103_TRM_Th17"
scAUC_obj_sub <- subset(scAUC_obj, Section == "TI" &
                          ann_level4_final == sel_ct)
# selected Regulons
diff_auc_ti <- readRDS(glue("{diff_path_it}scAUC_pseudo_pair_{ann_col}_II_vs_NN.rds"))
diff_auc_ti <- diff_auc_ti[!is.na(diff_auc_ti$padj),]
diff_auc_ti_trm17 <- subset(diff_auc_ti, cluster == sel_ct)
diff_auc_ti_trm17_sig <- subset(diff_auc_ti_trm17, padj < 0.05&
                                  grepl("\\(\\+\\)", Term))
diff_auc_ti_trm17_sig$TF <- str_split_i(diff_auc_ti_trm17_sig$Term, "\\(", 1)
diff_auc_ti_trm17_sig$eRegulon <- str_split_i(diff_auc_ti_trm17_sig$Term, "\\-", 1)
#
sel_tf <- c("BACH1",  "BATF", "MAF", "IKZF2", "STAT4",  "STAT5A", "STAT5B")
diff_auc_ti_sel <- subset(diff_auc_ti_trm17_sig, 
                          TF %in% sel_tf)
sel_er <- unique(diff_auc_ti_sel$Term)

scAUC_obj_sub_sel <- subset(scAUC_obj_sub, 
                        Section == "TI", 
                        features = sel_er)
scAUC_obj_sub_sel <- ScaleData(scAUC_obj_sub_sel)

# pseudobulk AUC matrix
AUC_mat_agg <- AverageExpression(scAUC_obj_sub_sel, 
                                 group.by = "Condition", 
                                 return.seurat = F, 
                                 assays = "AUC", 
                                 layer = "scale.data")[["AUC"]] %>% 
  as.data.frame()
AUC_mat_agg <- AUC_mat_agg[sort(sel_er), c("NN", "NU", "II"), drop = F]
rownames(AUC_mat_agg) <- gsub("\\-\\(", "_\\(", rownames(AUC_mat_agg))
## create symbol mat
symbol_matx <- matrix("", 
                      ncol = ncol(AUC_mat_agg), 
                      nrow = nrow(AUC_mat_agg))
dimnames(symbol_matx) <- dimnames(AUC_mat_agg)
for (nx in 1:nrow(diff_auc_ti_sel)) {
  
  termx <- diff_auc_ti_sel[nx, "Term"] %>%
    gsub("\\-\\(", "_\\(", .)
  posx <- diff_auc_ti_sel[nx, "avg_diff"] > 0
  sig_g <- diff_auc_ti_sel[nx, "contrast"] %>% 
    stringr::str_split(., "\\-", simplify = T)
  if (posx) {
    sig_g <- sig_g[,1:2, drop = T]
  } else {
    sig_g <- sig_g[,2:1, drop = T]
  }
  #
  if (symbol_matx[termx, sig_g[1]] == "-") {
    symbol_matx[termx, sig_g[1]] <- "+/-"
  } else {
    symbol_matx[termx, sig_g[1]] <- "+"
  }
  #
  if (symbol_matx[termx, sig_g[2]] == "+") {
    symbol_matx[termx, sig_g[2]] <- "+/-"
  } else {
    symbol_matx[termx, sig_g[2]] <- "-"
  }
}
symbol_matx[,"NN"] <- ""

## plot heatmap
clust_row <- nrow(AUC_mat_agg) > 1
png(glue("{plt_path}/heatAUC_sel_CD4TRM17.png"),
    height = 6, width = 4, 
    units = "in", res = 300)
pheatmap::pheatmap(AUC_mat_agg,
                   scale = "none", 
                   legend_labels = c("min","max"),
                   fontsize = 12, 
                   main = "Ileum",
                   display_numbers = symbol_matx,
                   fontsize_number = 15, 
                   treeheight_col = 0,
                   cluster_rows = F,
                   treeheight_row = 0,
                   cluster_cols = T, 
                   angle_col = 0) %>% print()
dev.off()

##### Fig. 4D and Supp 4.: Network #####
library(tidyverse)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(clusterProfiler)

use_de <- "DESeq2"
ann_col <- "ann_level4_final"
ctx <- "CD4_CD103_TRM_Th17"

deg_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{deg_path}TI/CD4T/{use_de}/")
## enrichment plot
enrich_gene_up <- readRDS(glue("{deg_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}_up.rds"))
enrich_dot_up <- dotplot(enrich_gene_up,
                         showCategory = 10) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  ggtitle("Enrichment of DEGs up-regulated") + 
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, color = "black"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10, color = "black"))
ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_up.png"),
       enrich_dot_up,
       height = 4, width = 6, units = "in", dpi = 300)
##
enrich_gene_down <- readRDS(glue("{deg_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}_down.rds"))
enrich_dot_down <- dotplot(enrich_gene_down,
                           showCategory = 10) +
  ggtitle("Enrichment of DEGs down-regulated") + 
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, color = "black"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10, color = "black"))
ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_down.png"),
       enrich_dot_down,
       height = 4, width = 6, units = "in", dpi = 300)
##
enrich_peak_up <- readRDS(glue("{dar_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}_up.rds"))
enrich_dot_up <- dotplot(enrich_peak_up,
                         showCategory = 10) +
  ggtitle("Enrichment of DARs up-regulated") + 
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, color = "black"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10, color = "black"))
ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_dar_up.png"),
       enrich_dot_up,
       height = 4, width = 6, units = "in", dpi = 300)
##
enrich_peak_down <- readRDS(glue("{dar_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}_down.rds"))
enrich_dot_down <- dotplot(enrich_peak_down,
                           showCategory = 10) +
  ggtitle("Enrichment of DARs down-regulated") + 
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, color = "black"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10, color = "black"))
ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_dar_down.png"),
       enrich_dot_down,
       height = 4, width = 6, units = "in", dpi = 300)

## network plot
sel_tf1 <- c("BATF", "BACH1", "MAF", 
             "IKZF2", "STAT4", "STAT5A", "STAT5B") %>%
  paste0(., "(+)")
#
sel_tf2 <- intersect(enrich_gene_down@result$ID[enrich_gene_down@result$Count > 5 & 
                                                  enrich_gene_down@result$p.adjust < 0.05], 
                     enrich_peak_down@result$ID[enrich_peak_down@result$Count > 10 & 
                                                  enrich_peak_down@result$p.adjust < 0.05])
sel_tf <- union(sel_tf1, sel_tf2)
## load peak2gene
p2g_cd103trm17 <- readRDS(glue("{diff_path}peak2gene/TI/LinkPeaks_links_CD103_TRM_Th17_DEG_II_vs_NN.rds")) %>%
  as.data.frame()
## load eRegulon annotation of CD4T
e_regulon_filter <- readRDS("03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/e_regulon_filter_rho_add.rds")
e_regulon_filter$eRegulon_simplify <- str_split_i(e_regulon_filter$Gene_signature_name_simplify, "_", 1)
e_regulon_filter$Region <- gsub("\\:", "-", e_regulon_filter$Region)

## load DEG and DAR
use_de <- "DESeq2"
ann_col <- "ann_level4_final"
diff_path <- "03_output/04_Diff/"
deg_path_it <- glue("{diff_path}DEG/TI/CD4T/{use_de}/")
dar_path_it <- glue("{diff_path}DAR/TI/CD4T/{use_de}/")
#
deg_cd103trm17_ti <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4_CD103_TRM_Th17"]]
deg_cd103trm17_ti$Sig <- F
deg_cd103trm17_ti$Sig[which(deg_cd103trm17_ti$padj < 0.05 &
                              abs(deg_cd103trm17_ti$log2FC) > 0.25)] <- T
#
dar_cd103trm17_ti <- readRDS(glue("{dar_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4_CD103_TRM_Th17"]]
dar_cd103trm17_ti$Sig <- F
dar_cd103trm17_ti$Sig[which(dar_cd103trm17_ti$padj < 0.05 &
                              abs(dar_cd103trm17_ti$log2FC) > 0.25)] <- T
##
e_regulon_f_sel <- subset(e_regulon_filter, eRegulon_simplify %in% sel_tf &
                            (Gene %in% deg_cd103trm17_ti$Term[deg_cd103trm17_ti$Sig] | 
                               Region %in% dar_cd103trm17_ti$Term[dar_cd103trm17_ti$Sig]))
#
deg_cd103trm17_ti_sel <- subset(deg_cd103trm17_ti, 
                                Sig & Term %in% e_regulon_f_sel$Gene)
deg_cd103trm17_ti_sel$Group <- ifelse(deg_cd103trm17_ti_sel$log2FC > 0,
                                      "DEG-Up", "DEG-Down")
#
dar_cd103trm17_ti_sel <- subset(dar_cd103trm17_ti, 
                                Term %in% e_regulon_f_sel$Region)

dar_cd103trm17_ti_sel$Group <- ifelse(!dar_cd103trm17_ti_sel$Sig, "DAR-No.sig",
                                      ifelse(dar_cd103trm17_ti_sel$log2FC > 0,
                                             "DAR-Up", "DAR-Down"))
diff_cd103trm17_ti_sel <- rbind(deg_cd103trm17_ti_sel, dar_cd103trm17_ti_sel)

##
nodeDf <- data.frame(id = c(unique(e_regulon_f_sel$eRegulon_simplify), 
                            diff_cd103trm17_ti_sel$Term),
                     label = c(unique(e_regulon_f_sel$eRegulon_simplify), 
                               diff_cd103trm17_ti_sel$Term),
                     Group = c(rep("TF", length(unique(e_regulon_f_sel$eRegulon_simplify))),
                               diff_cd103trm17_ti_sel$Group))
##
nodeDf$FC <- diff_cd103trm17_ti_sel$log2FC[match(nodeDf$id, diff_cd103trm17_ti_sel$Term)]
nodeDf$Group <- factor(nodeDf$Group, levels = c("TF", 
                                                "DAR-Up", "DAR-Down", "DAR-No.sig",
                                                "DEG-Up", "DEG-Down"))
nodeDf$FC[nodeDf$Group == "TF"] <- 2
nodeDf$Type <- case_when(nodeDf$Group == "TF" ~ "TF",
                         grepl("DAR", nodeDf$Group) ~ "Region",
                         .default = "Gene")
nodeDf$labelsize <- case_when(nodeDf$Type == "TF" ~ 4,
                              nodeDf$Type == "Region" ~ 0,
                              .default = 2)
nodeDf$label[nodeDf$Type == "Region"] <- NA

##
linkDf <- rbind(with(e_regulon_f_sel[!duplicated(paste0(e_regulon_f_sel$eRegulon_simplify, 
                                                        ":",
                                                        e_regulon_f_sel$Region)), ],
                     data.frame(from = eRegulon_simplify,
                                to = Region,
                                R2 = 0.2,
                                Direction = "T2R")),
                with(e_regulon_f_sel[!duplicated(paste0(e_regulon_f_sel$Region, 
                                                        ":",
                                                        e_regulon_f_sel$Gene)), ],
                     data.frame(from = Region,
                                to = Gene,
                                R2 = rho_R2G,
                                Direction = ifelse(rho_R2G > 0, "R2G-POS", "R2G-NEG"))))
linkDf <- subset(linkDf, from %in% nodeDf$id &
                   to %in% nodeDf$id)

##
nodeDf = nodeDf %>%
  arrange(Group, abs(FC)) 
##
mygraph = graph_from_data_frame(linkDf, 
                                vertices = nodeDf, 
                                directed = FALSE) 
#
png(glue("{plt_path}/net_tf_CD103TRM17.png"),
    height = 6, width = 12, units = "in", res = 600)
set.seed(seed_use)
ggraph(mygraph, layout = 'stress') +  
  geom_edge_link(aes(edge_color = Direction, alpha = R2, edge_width = R2)) +
  geom_node_point(aes(fill = Group, shape = Type, color = Group, size = abs(FC)), 
                  alpha = 0.9) +
  geom_node_text(aes(label = label), 
                 repel = T, 
                 size = nodeDf$labelsize, 
                 segment.size = 0, segment.color = NA) +
  scale_size_continuous(range = c(1, 4),
                        breaks = c(1:4),
                        labels = c(as.character(1:4))) +
  scale_color_manual(values = c("TF" = "brown", 
                                "DEG-Down" = "#2C69B0", 
                                "DEG-Up" = "#F02720",
                                "DAR-Down" = "skyblue", 
                                "DAR-Up" = "pink",
                                "DAR-No.sig" = "grey80")) +
  scale_edge_width_continuous(range = c(0.2, 0.6)) +
  scale_edge_color_manual(values = c("R2G-NEG" = "purple", 
                                     "R2G-POS" = "orange",
                                     "T2R" = "darkgreen")) +
  theme_void() + 
  theme(legend.position = "right")
dev.off()

##### Fig. 4E: UMAP #####
umap_plt <- DimPlot(sc_obj, 
                    reduction = "wnn.umap",
                    cols = cd4cd103trm_cols,
                    label = F, 
                    raster = F)+
  xlab("wnnUMAP1") + ylab("wnnUMAP2")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#
ggsave(file = glue("{plt_path}UmapPlot_all.png"),
       umap_plt,
       width = 10, height = 6,units = "in", dpi = 300, limitsize = F)

##### Fig. 4F: cNMF #####
sub_idx <- "CD103TRM"
n_k <- 10
nmf_outpath <- "03_output/03_clustering/CD4_CD103_TRM_Th17/cNMF/"
nmf_usage <- fread2(file = glue("{nmf_outpath}/CD4_CD103_TRM_harmony_usage_k{n_k}.txt"),
                    header = T) %>%
  tibble::column_to_rownames(., "V1")
top_gene <- fread2(glue("{nmf_outpath}/CD4_CD103_TRM_harmony_top_genes_k{n_k}.txt"))
sc_obj[["NMF"]] <- CreateAssayObject(t(nmf_usage),   
                                     min.cells = 0,
                                     min.features = 0)
##
ft_nmf <- FeaturePlot(sc_obj,
                      reduction = "wnn.umap",
                      raster = F,
                      features = colnames(nmf_usage),
                      max.cutoff = "q95", min.cutoff = "q5",
                      ncol = 4)
ft_nmf <- lapply(ft_nmf, function(ft_nmf14x){
  
  ft_nmf14x + 
    xlab("wnnUMAP1") + ylab("wnnUMAP2")+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_rect(size = 1),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
})
#
ggsave(glue("{plt_path}feature_nmf.png"), 
       wrap_plots(ft_nmf, ncol = 4),
       height = 9, width = 12, units = "in", dpi = 300)

##
GEP_order <- paste0("GEP",
                    c(1, 5, 7, 10, 9, 3, 6, 8, 4, 2))
sc_obj_bulk <- Seurat::AverageExpression(sc_obj,
                                         # features = GEP_order,
                                         group.by = "ann_level5_final",
                                         assays = "NMF",
                                         return.seurat = T,
                                         layer = "data")

ann_level <- gsub("\\-", " ", colnames(sc_obj_bulk))
ann_col <- data.frame(Celltype = factor(ann_level,
                                        levels = ann_level),
                      row.names = colnames(sc_obj_bulk))
anno_colors <- cd4cd103trm_cols
names(anno_colors) <- ann_level

ann_row <- data.frame(
  `Top genes` = lapply(GEP_order, function(x){
    
    paste0(top_gene[1:10, x] %>% paste(collapse = ", "))
    
  }) %>% unlist,
  row.names = GEP_order
)
write.table(ann_row, glue("{plt_path}heat_score_{n_k}_topgenes.txt"), 
            sep = "\t", col.names = T, row.names = T, quote = F)
#
png(glue("{plt_path}heat_score_{n_k}.png"),
    height = 6, width = 9, units = "in", res = 300)
pheatmap::pheatmap(sc_obj_bulk[["NMF"]]$data[GEP_order,],
                   scale = "row",
                   legend = F,
                   annotation_col = ann_col,
                   annotation_colors = list(Celltype = anno_colors),
                   border_color = NA,
                   fontsize = 10,
                   fontsize_row = 15,
                   show_colnames = F,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()

##### Fig. 4G: ct-specific eRegulon #####
sp_outpath <- "03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/"
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))

scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
scAUC_obj_trm <- subset(scAUC_obj, cells = rownames(sc_meta_cd4cd103trm17))
scAUC_obj_trm@meta.data <- sc_meta_cd4cd103trm17[colnames(scAUC_obj_trm),]
#
Idents(scAUC_obj_trm) <- scAUC_obj_trm$ann_level5_final
er_gene <- rownames(scAUC_obj_trm)[grep("g\\)", rownames(scAUC_obj_trm))]
diff_auc <- FindAllMarkers(scAUC_obj_trm, 
                           features = er_gene,
                           only.pos = T)
saveRDS(diff_auc, file = glue("{data_path}diff_auc.rds"))
#
sig_dif <- diff_auc %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05)
sig_dif <- sig_dif[order(-abs(sig_dif$avg_log2FC), sig_dif$p_val_adj),]
sig_dif <- sig_dif[!duplicated(sig_dif$gene),]

sig_diff <- sig_dif %>%
  group_by(cluster) %>%
  top_n(5, -p_val_adj) %>%
  ungroup()
sig_diff <- sig_diff[order(sig_diff$cluster, sig_diff$p_val_adj, -abs(sig_diff$avg_log2FC)),]
sig_marker <- str_split_i(sig_diff$gene, "-\\(", 1)

##
heatdot_CD4trm <- fread2(glue("{sp_outpath}heatdot_CD103trm17_ann_level5_final.txt"))[,-1]
colnames(heatdot_CD4trm) <- c("celltype", "Gene-based AUC", "eRegulon_name", "Region-based AUC")
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))
e_regulon_name_trans$heatdot_name <- stringr::str_split_i(e_regulon_name_trans$signature, "_\\(", 1)
heatdot_CD4trm$signature_name_simplify <- 
  e_regulon_name_trans$signature_name_simplify[match(heatdot_CD4trm$eRegulon_name,
                                                     e_regulon_name_trans$heatdot_name)] %>%
  str_split_i(., "_", 1)

##
sig_marker <- sig_marker[which(sig_marker %in% heatdot_CD4trm$signature_name_simplify)]
heatdot_CD4trm_sub <- subset(heatdot_CD4trm, signature_name_simplify %in% sig_marker)
heatdot_CD4trm_sub$signature_name_simplify <- factor(heatdot_CD4trm_sub$signature_name_simplify,
                                                     levels = rev(sig_marker))
heatdot_CD4trm_sub$celltype <- gsub("_|\\-", " ", heatdot_CD4trm_sub$celltype) %>%
  factor(.,
         levels = gsub("_", " ", levels(scAUC_obj_trm$ann_level5_final)))
heatdot_CD4trm_sub$celltype <- droplevels(heatdot_CD4trm_sub$celltype)
heat_dot <- ggplot(heatdot_CD4trm_sub, 
                   aes(celltype, signature_name_simplify))+ 
  geom_tile(mapping = aes(fill = `Gene-based AUC`)) + 
  scale_fill_distiller(type = 'div', palette = 'RdYlBu') + 
  geom_point(mapping = aes(size = `Region-based AUC`),
             colour = "black") + 
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.text.x = element_blank(),
        legend.position = "left",
        legend.justification = c(0, 1),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10, color = "black"))
#
df_annot <- data.frame(celltype = levels(heatdot_CD4trm_sub$celltype))
p_annot = ggplot(df_annot, aes(x = celltype, y = 1, fill = celltype)) +
  scale_fill_manual(values = cd4cd103trm_cols) + 
  geom_tile() +
  theme_void() +
  theme(legend.position = "none")

ggsave(glue("{plt_path}heat_dot_CD4trm17_ann_level5_final.png"),
       wrap_plots(list(p_annot, heat_dot), ncol = 1, heights = c(0.1, 8)),
       height = 16, width = 8, units = "in", dpi = 300)

##### Fig. 4H: miloR #####
source("code/FUNCTION/miloR.R")
pair_sample_set <- readRDS("03_output/03_clustering/pair_sample_set_1205.rds")
sc_obj$group_pair <- pair_sample_set$group_pair[match(sc_obj$Sample_ID_exp,
                                                      pair_sample_set$Sample_ID_exp)]
## compare between Conditions
all_pair_set <- c("II vs NN in Ileum", "NU vs NN in Ileum",
                  "II vs NN in Colon", "NU vs NN in Colon")
names(all_pair_set) <- c("TI_II_NN", "TI_NU_NN",
                         "Colon_II_NN", "Colon_NU_NN")
for (pair_setx in names(all_pair_set)) {
  
  sel_cellx <- colnames(sc_obj)[grep(pair_setx, sc_obj@meta.data$group_pair)]
  sc_obj_sub <- subset(sc_obj, 
                       cells = sel_cellx) %>%
    DietSeurat(., 
               graphs = "wsnn", 
               dimreducs = c("wnn.umap", "harmony_lsi", "harmony_SCT"))
  bulk_meta_sub <- sc_obj_sub@meta.data
  bulk_meta_sub <- bulk_meta_sub[!duplicated(bulk_meta_sub$Sample_ID_exp),]
  pair_sample_sub <- bulk_meta_sub$Sample_exp[duplicated(bulk_meta_sub$Sample_exp)]
  sc_obj_sub <- subset(sc_obj_sub, 
                       Sample_exp %in% pair_sample_sub)
  sc_obj_sub$Condition <- droplevels(sc_obj_sub$Condition)
  sc_obj_sub$ann_level5_final <- factor(sc_obj_sub$ann_level5_final,
                                        levels = levels(sc_obj_sub$ann_level5_final),
                                        labels = gsub("_", " ",
                                                      levels(sc_obj_sub$ann_level5_final)))
  #
  res_milo <- milor.seurat(seurat_obj = sc_obj_sub,
                           n_dims = 20,
                           k_nn = 30,
                           group_col = "Condition",
                           cov_col = "Sample_exp",
                           block_col = NULL,
                           sample_col = "Sample_ID_exp",
                           ann_col = "ann_level5_final",
                           use_reduc = "harmony_SCT",
                           plot_reduc = "wnn.umap",
                           seed_use = seed_use,
                           n_core = 10)
  res_milo$plt_milo <- res_milo$plt_milo + 
    ggtitle(all_pair_set[pair_setx]) +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  ##
  plt_milo_comb <- (res_milo$plt_milo | res_milo$plt_milo2) + 
    plot_layout(widths = c(5, 5))
  ggsave(glue("{plt_path}traj_milo_{pair_setx}.png"),
         plt_milo_comb,
         height = 7, width = 16, units = "in", dpi = 300)
}

## compare between Sections
pair_setx <- "TI_Colon_NN"
sel_cellx <- colnames(sc_obj)[grep(pair_setx, sc_obj@meta.data$group_pair)]
sc_obj_sub <- subset(sc_obj, 
                     cells = sel_cellx) %>%
  DietSeurat(., 
             graphs = "wsnn", 
             dimreducs = c("wnn.umap", "harmony_lsi", "harmony_SCT"))
bulk_meta_sub <- sc_obj_sub@meta.data
bulk_meta_sub <- bulk_meta_sub[!duplicated(bulk_meta_sub$Sample_ID_exp),]
pair_sample_sub <- bulk_meta_sub$Sample_exp[duplicated(bulk_meta_sub$Sample_exp)]
sc_obj_sub <- subset(sc_obj_sub, 
                     Sample_exp %in% pair_sample_sub)
sc_obj_sub$Section <- factor(sc_obj_sub$Section, levels = c("TI", "Colon"))
sc_obj_sub$ann_level5_final <- factor(sc_obj_sub$ann_level5_final,
                                      levels = levels(sc_obj_sub$ann_level5_final),
                                      labels = gsub("_", " ",
                                                    levels(sc_obj_sub$ann_level5_final)))
#
res_milo <- milor.seurat(seurat_obj = sc_obj_sub,
                         n_dims = 30,
                         k_nn = 40,
                         group_col = "Section",
                         cov_col = "Sample_exp",
                         block_col = NULL,
                         sample_col = "Sample_ID_exp",
                         ann_col = "ann_level5_final",
                         use_reduc = "harmony_SCT",
                         plot_reduc = "wnn.umap",
                         seed_use = seed_use,
                         n_core = 10)
#
res_milo$plt_milo <- res_milo$plt_milo + 
  ggtitle("NN in Ileum vs Colon") + 
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
##
plt_milo_comb <- (res_milo$plt_milo | res_milo$plt_milo2) + 
  plot_layout(widths = c(6, 4))
ggsave(glue("{plt_path}traj_milo_{pair_setx}.png"),
       plt_milo_comb,
       height = 7, width = 14, units = "in", dpi = 300)


##### Fig. 4I: Diff prop ######
sample_col <- "Sample_ID_exp"
block_col <- "Sample_exp"
ann_col <-  "ann_level5_final"
group_col <- "Condition"
section_col <- "Section"

## 0. format proportion table
bulk_meta <- sc_meta_cd4cd103trm17[!duplicated(sc_meta_cd4cd103trm17[[sample_col]]),]
bulk_meta$group_comb  <- paste0(bulk_meta[[group_col]], ":", bulk_meta[[section_col]])
sample_list <- split(bulk_meta$Sample_exp, f = bulk_meta$group_comb)
#
tab_prop <- table(sc_meta_cd4cd103trm17[[ann_col]], sc_meta_cd4cd103trm17[[sample_col]]) %>% 
  prop.table(2) %>% as.data.frame()
colnames(tab_prop) <- c("cluster", "sample", "Proportion")
tab_prop$Proportion <- tab_prop$Proportion * 100
tab_prop$block <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                            block_col]
tab_prop$group <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                            group_col]
tab_prop$Section <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                              section_col] %>%
  factor(., levels = c("TI", "Colon"),
         labels = c("Ileum", "Colon"))
tab_prop$group <- factor(tab_prop$group, 
                         levels = c("NN", "NU", "II"))

## 1. paired test between Conditions
all_ct <- unique(tab_prop$cluster)
ppair_prop_df1 <- lapply(1:4, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  tab_prop_use <- subset(tab_prop, block %in% sample_pair_use)
  tab_prop_use$group <- droplevels(tab_prop_use$group)
  group_level <- levels(tab_prop_use$group)
  ppairx <- lapply(all_ct, function(ctx){
    #
    tab_prop_usex <- subset(tab_prop_use, cluster == ctx)
    pair_dfx <- split(tab_prop_usex, f = tab_prop_usex$block) %>%
      lapply(., function(x){
        propx <- x$Proportion[match(group_level, x$group)]
      }) %>% Reduce("rbind", .)
    wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
    
  }) %>% unlist()
  names(ppairx) <- c(all_ct)
  return(ppairx)
  
}) %>% Reduce("cbind", .) %>% as.data.frame()

# adjust p value for each set of comparison
ppair_prop_adj_df1 <- apply(ppair_prop_df1, 2, function(x){
  p.adjust(x, method = "BH")
}) %>% as.data.frame()
colnames(ppair_prop_df1) <- colnames(ppair_prop_adj_df1) <- names(pair_list)[1:4]

## 2.1 plot for Ileum
sel_ct <- c("CD4_CD103_TRM_Th17", "CD4_CD103_TRM_Th17_CREM", "CD4_CD103_TRM_Th17_HSP")
lab_facet <- c(1)
box_pair_prop_list1_TI <- lapply(sel_ct, function(ctx){
  
  plt_dfx <- subset(tab_prop, cluster == ctx & Section == "Ileum")
  plt_dfx$Condition <- plt_dfx$group
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df1[ctx,2:1], scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1, max_propx * 1.2),
                         x = c(1.5, 2.5),
                         xmin = c(1, 1),
                         xmax = c(2, 3),
                         group2 = c("NN", "NN"),
                         group1 = c("NU", "II"),
                         padj_sig = padj)
  
  ylabx <- ifelse(ctx %in% all_ct[lab_facet], "Proportion (%) in Ileum", "")
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Condition, y = Proportion, color = Condition)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(gsub("_", " ", ctx)) +
    xlab("") + ylab(ylabx) + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})

ggsave(glue("{plt_path}box_pair_level5_Ileum_sel.png"),
       patchwork::wrap_plots(box_pair_prop_list1_TI, nrow = 1),
       height = 5, width = 10, units = "in", dpi = 300)

## 2.2 plot for Colon
lab_facet <- c(1)
box_pair_list1_Colon <- lapply(all_ct, function(ctx){
  
  plt_dfx <- subset(tab_prop, cluster == ctx & Section == "Colon")
  plt_dfx$Condition <- plt_dfx$group
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df1[ctx, 4:3], scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1, max_propx * 1.2),
                         x = c(1.5, 2.5),
                         xmin = c(1, 1),
                         xmax = c(2, 3),
                         group2 = c("NN", "NN"),
                         group1 = c("NU", "II"),
                         padj_sig = padj)
  
  ylabx <- ifelse(ctx %in% all_ct[lab_facet], "Proportion (%) in Ileum", "")
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Condition, y = Proportion, color = Condition)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(gsub("_", " ", ctx)) +
    xlab("") + ylab("Proportion (%) in Colon") + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})

ggsave(glue("{plt_path}box_pair_level4_Colon.png"),
       patchwork::wrap_plots(box_pair_list1_Colon, nrow = 2),
       height = 7, width = 16, units = "in", dpi = 300)

## 3. paired test between Sections
ppair_prop_df2 <- lapply(5, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  tab_prop_use <- subset(tab_prop, block %in% sample_pair_use)
  tab_prop_use$Section <- droplevels(tab_prop_use$Section)
  group_level <- levels(tab_prop_use$Section)
  ppairx <- lapply(all_ct, function(ctx){
    #
    tab_prop_usex <- subset(tab_prop_use, cluster == ctx)
    pair_dfx <- split(tab_prop_usex, f = tab_prop_usex$block) %>%
      lapply(., function(x){
        propx <- x$Proportion[match(group_level, x$Section)]
      }) %>% Reduce("rbind", .)
    wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
    
  }) %>% unlist()
  names(ppairx) <- c(all_ct)
  return(ppairx)
  
}) %>% Reduce("cbind", .) %>% as.data.frame()

# adjust p value for each set of comparison
ppair_prop_adj_df2 <- apply(ppair_prop_df2, 2, function(x){
  p.adjust(x, method = "BH")
}) %>% as.data.frame()
colnames(ppair_prop_df2) <- colnames(ppair_prop_adj_df2) <- names(pair_list)[5]

## 3.1 plot
box_pair_list2 <- lapply(all_ct, function(ctx){
  
  plt_dfx <- subset(tab_prop, cluster == ctx & group == "NN")
  plt_dfx$Condition <- plt_dfx$Section
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df2[ctx, 1], scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1),
                         x = c(1.5),
                         xmin = c(1),
                         xmax = c(2),
                         group2 = c("Ileum"),
                         group1 = c("Colon"),
                         padj_sig = padj)
  
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Section, y = Proportion, color = Section)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("#0073C2FF", "#E41A1C")) + 
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(ctx) +
    xlab("") + ylab("Proportion (%) in NN") + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})

ggsave(glue("{plt_path}box_pair_level4_NN.png"),
       patchwork::wrap_plots(box_pair_list2, nrow = 2),
       height = 7, width = 16, units = "in", dpi = 300)

##### Supp Fig. 4A: Circle Manha plot #####
library(GenomicRanges)
library(CMplot)

project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

## 0.1 diff path
diff_path <- "03_output/04_Diff/"
deg_path_it <- glue("{diff_path}DEG/TI/CD4T/DESeq2/")
dar_path_it <- glue("{diff_path}DAR/TI/CD4T/DESeq2/")
p2g_path_it <- glue("{diff_path}peak2gene/TI/")

## 0.2 SNP position GRCh 38
summ_df <- fread2("/ix1/wchen/xiangyu/Ref_data/GWAS/Raw/28067908-GCST004132-EFO_0000384.h.tsv.gz")
summ_dff <- subset(summ_df, 
                   hm_chrom %in% c(1:22) &
                     hm_other_allele %in% c("T", "G", "C", "A") &
                     hm_effect_allele %in% c("T", "G", "C", "A") &
                     standard_error > 0,
                   select = c("hm_chrom", "hm_pos", "hm_rsid", "p_value"))
colnames(summ_dff) <- c("chr", "pos", "rsid", "pval_gwas")
gr_gwas <- GRanges(seqnames = paste0("chr", summ_dff$chr), 
                   ranges = IRanges(start = summ_dff$pos, 
                                    width = 1))
## 0.3 Gene coord
gene_coord <- fread2("/ix1/wchen/xiangyu/Ref_data/genome/cellranger/arc_GRCh38/genes/gene_coordinates.tsv")
gene_coord$chr <- gene_coord$seqnames %>% 
  gsub("chr", "", .) %>% as.integer()

## 1. DEG
deg_CD4 <- readRDS(glue("{deg_path_it}pseudo_bulk_ann_level4_final_II_vs_NN.rds"))[["CD4_CD103_TRM_Th17"]]
deg_CD4$chr <- gene_coord$seqnames[match(deg_CD4$Term, gene_coord$gene_name)] %>%
  gsub("chr", "", .) %>% as.integer()
# map to nearest SNP pos
gr_gene <- GRanges(seqnames = gene_coord$seqnames,
                   ranges = IRanges(start = gene_coord$start,
                                    end = gene_coord$end))
idx_gene <- IRanges::nearest(gr_gene, gr_gwas)
gene_coord$pos <- summ_dff$pos[idx_gene]
gene_coord$rsid_nearest <- summ_dff$rsid[idx_gene]
deg_CD4$rsid_nearest <- gene_coord$rsid_nearest[match(deg_CD4$Term, gene_coord$gene_name)]

## 2. MAGMA
magma_df <- fread2("03_output/01_scDRS/01_magma/CD.genes.out")
magma_zdf <- fread2("03_output/01_scDRS/01_magma/zscore_CD.tsv")
magma_df$Term <- magma_zdf$GENE
deg_CD4$p_magma <- magma_df$P[match(deg_CD4$Term, magma_df$Term)]

## 3. DAR
dar_CD4 <- readRDS(glue("{dar_path_it}pseudo_bulk_ann_level4_final_II_vs_NN.rds"))[["CD4_CD103_TRM_Th17"]]
# map to nearest SNP pos
gr_peak <- StringToGRanges(dar_CD4$Term)
idx_peak <- IRanges::nearest(gr_peak, gr_gwas)
dar_CD4$pos <- summ_dff$pos[idx_peak]
dar_CD4$rsid_nearest <- summ_dff$rsid[idx_peak]
#
dar_CD4 <- subset(dar_CD4, !is.na(padj) & !is.na(pos))
dar_CD4 <- dar_CD4[order(dar_CD4$padj),]
dar_CD4 <- dar_CD4[!duplicated(dar_CD4$rsid_nearest),]

## 4. combine DEG, DAR, and MAGMA
summ_dff$Gene <- deg_CD4$Term[match(summ_dff$rsid, deg_CD4$rsid_nearest)]
summ_dff$log2FC_gene <- deg_CD4$log2FC[match(summ_dff$rsid, deg_CD4$rsid_nearest)]
summ_dff$padj_gene <- deg_CD4$padj[match(summ_dff$rsid, deg_CD4$rsid_nearest)]
summ_dff$p_magma <- deg_CD4$p_magma[match(summ_dff$rsid, deg_CD4$rsid_nearest)]

summ_dff$Peak <- dar_CD4$Term[match(summ_dff$rsid, dar_CD4$rsid_nearest)]
summ_dff$log2FC_peak <- dar_CD4$log2FC[match(summ_dff$rsid, dar_CD4$rsid_nearest)]
summ_dff$padj_peak <- dar_CD4$padj[match(summ_dff$rsid, dar_CD4$rsid_nearest)]

# summ_dff <- subset(summ_dff, Peak %in% df_links$peak)

plt_df <- data.frame(SNP = summ_dff$rsid,
                     Chromosome = summ_dff$chr,
                     Position = summ_dff$pos,
                     DAP = summ_dff$padj_peak,
                     DEG = summ_dff$padj_gene)
plt_df <- subset(plt_df, 
                 !is.na(DAP) | !is.na(DEG))

## 5. parameters for plot
# dot color for three modalities
col_mat <- matrix(c("royalblue4", "darksalmon",
                    "grey", "lightblue"), 
                  nrow = 2, byrow = T)

# labels genes: DEGs & MAGMA Sig
top_gene_magma <- deg_CD4$Term[which(deg_CD4$padj < 0.05 & 
                                       abs(deg_CD4$log2FC) > 0.25 &
                                       deg_CD4$p_magma < 1E-5)]
top_gene_magma_pos <- deg_CD4$rsid_nearest[match(top_gene_magma, deg_CD4$Term)]
gene_high_col_magma <- rep("brown", length(top_gene_magma))

# other highlight genes: DEGs & link to DAPs
links_gr <- readRDS(glue("{p2g_path_it}LinkPeaks_links_CD103_TRM_Th17_DEG_II_vs_NN.rds"))
df_links <- as.data.frame(links_gr)
df_links$padj <- p.adjust(df_links$pvalue, method = "BH")
dff_links <- subset(df_links, df_links$padj < 0.05 & df_links$score > 0.1)
dff_links_sig <- subset(dff_links, peak %in% 
                          dar_CD4$Term[which(dar_CD4$padj < 0.05 & abs(dar_CD4$log2FC) > 0.25)] &
                          gene %in% 
                          deg_CD4$Term[which(deg_CD4$padj < 0.05 & abs(deg_CD4$log2FC) > 0.25)])
#
top_gene_link <- setdiff(unique(dff_links_sig$gene), top_gene_magma)
top_gene_link_pos <- deg_CD4$rsid_nearest[match(top_gene_link, deg_CD4$Term)]
gene_high_col_link <- rep("orange", length(top_gene_link))

# highlight peaks: DEGs & link to DAPs
top_peak_link <- unique(dff_links_sig$peak)
top_peak_link_pos <- dar_CD4$rsid_nearest[match(top_peak_link, dar_CD4$Term)]

# show top one peak for each gene
dff_links_sig <- dff_links_sig[order(dff_links_sig$gene, dff_links_sig$padj, decreasing = F),]
dff_links_sig_unique <- dff_links_sig[!duplicated(dff_links_sig$gene),]
top_peak_link[!top_peak_link %in% dff_links_sig_unique$peak] <- NA

## 6. Circle plot
plt_df_md <- plt_df
colnames(plt_df_md)[-c(1:3)] <- c("DAP in Ileal II vs NN",
                               "DEG in Ileal II vs NN")
CMplot(plt_df_md, 
       type = "p",
       plot.type = "m",
       r = 1,
       col = col_mat,
       cir.chr.h = 1.5,
       amplify = F,
       signal.line = 0.2,
       threshold = 0.05,
       threshold.col = "red",
       threshold.lty = 2,   
       highlight = list("DAP in Ileal II vs NN" = top_peak_link_pos,
                        "DEG in Ileal II vs NN" = c(top_gene_magma_pos, top_gene_link_pos)),
       highlight.col = list("DAP in Ileal II vs NN" = "orange",
                            "DEG in Ileal II vs NN" = c(gene_high_col_magma, gene_high_col_link)),
       highlight.text = list("DAP in Ileal II vs NN" = top_peak_link,
                             "DEG in Ileal II vs NN" = c(top_gene_magma, top_gene_link)),
       highlight.cex = 1,
       highlight.text.cex = 1.5,
       highlight.text.col = "black",
       outward = T,
       file = "png",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 16,
       height = 8,
       axis.cex = 1.5,
       lab.cex = 2,
       legend.cex = 2.5,
       multracks = T)
system(glue("mv -f 'Multi-tracks_Manhtn.DAP in Ileal II vs NN_DEG in Ileal II vs NN.png' {plt_path}Manhtn.DAP.DEG.png"))

##### Supp Fig. 4: Heatmap for ct-specific genes #####
DefaultAssay(sc_obj) <- "RNA"
Idents(sc_obj) <- sc_obj$ann_level5_final
ct_markers <- FindAllMarkers(sc_obj, 
                             only.pos = T)
ct_markers <- ct_markers[order(ct_markers$cluster),]
saveRDS(ct_markers, file = glue("{data_path}ct_markers.rds"))
##
ct_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
sc_obj <- NormalizeData(sc_obj)
sc_obj_bulk <- Seurat::AverageExpression(sc_obj,
                                         features = unique(top5$gene),
                                         group.by = "ann_level5_final",
                                         assays = "RNA",
                                         return.seurat = T,
                                         layer = "data")
ann_level <- gsub("\\-", " ", colnames(sc_obj_bulk))
ann_col <- data.frame(Celltype = factor(ann_level,
                                        levels = ann_level),
                      row.names = colnames(sc_obj_bulk))
anno_colors <- cd4cd103trm_cols
names(anno_colors) <- ann_level
##
png(glue("{plt_path}pseudo_heat_ann.png"),
    height = 8, width = 6, units = "in", res = 300)
pheatmap::pheatmap(sc_obj_bulk[["RNA"]]$data[unique(top5$gene),],
                   scale = "row",
                   annotation_col = ann_col,
                   annotation_colors = list(Celltype = anno_colors),
                   border_color = NA,
                   fontsize = 8,
                   show_colnames = F,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()

##### Supp Fig. 4: pie plot #####
ann_col <-  "ann_level5_final"
group_col <- "Condition"
section_col <- "Section"
##
sc_meta_cd4cd103trm17$group_comb <- paste0(sc_meta_cd4cd103trm17[[group_col]], ":", sc_meta_cd4cd103trm17[[section_col]])
tab_prop <- table(sc_meta_cd4cd103trm17[[ann_col]], sc_meta_cd4cd103trm17$group_comb) %>% 
  prop.table(2) %>% as.data.frame()
colnames(tab_prop) <- c("Cell type", "group_comb", "Proportion")
tab_prop$Proportion <- tab_prop$Proportion * 100
tab_prop$Section <- str_split_i(tab_prop$group_comb, ":", 2)
tab_prop$Group <- str_split_i(tab_prop$group_comb, ":", 1)

tab_prop$Section <- factor(tab_prop$Section, 
                         levels = c("Colon", "TI"),
                         labels = c("Colon", "Ileum"))
tab_prop$Group <- factor(tab_prop$Group, 
                         levels = c("NN", "NU", "II"))
##
pie_plt <- ggplot(tab_prop, aes(x = 3, y = Proportion, fill = `Cell type`)) + 
  geom_col(width = 1.5, color = NA) + 
  facet_grid(Group ~ Section, switch = "y") + 
  coord_polar(theta = "y") + 
  xlim(c(0.2, 3.8)) + 
  scale_fill_manual(values = cd4cd103trm_cols) + 
  theme_void()+
  theme(strip.text = element_text(size = 15, face = "bold"), 
        legend.position = "right")
ggsave(file = glue("{plt_path}pie_plt_cd4cd103trm17.png"),
       pie_plt,
       width = 8, height = 6,units = "in", dpi = 300)

##### Fig. 4A: volcano plot #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD4_CD103_TRM_Th17/{use_de}/")
ann_col <- "ann_level5_final"
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/volca_plot.R")
source("code/FUNCTION/diff/Enrichment.R")
##
lapply(c("CD4_CD103_TRM_Th17_HSP", "CD4_CD103_TRM_Th17_CREM"), function(sel_ct1){
  
  # sel_ct1 <- "CD4_CD103_TRM_Th17_HSP"
  deg_df_it_cd4x <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[sel_ct1]]
  deg_df_it_cd4x <- subset(deg_df_it_cd4x, !is.na(deg_df_it_cd4x$padj))
  vc_plt_cd4cd103trm17 <- volca_function(deg_df = deg_df_it_cd4x,
                                         term_col = "Term",
                                         fc_col = "log2FC",
                                         p_col = "padj",
                                         n_top = 30,
                                         log_thresh = 0.25, 
                                         top_log_thresh = 0.5,
                                         group_label = c("Down regulation", "Others", "Up regulation"),
                                         valco_col = c4a("classic_blue_red12", 3),
                                         highlight_col = "orange") + 
    ggtitle(gsub("_", " ", sel_ct1)) + 
    xlab(bquote(~Log[2]~"(Fold Change) in Ileal II vs NN")) +
    theme(plot.title = element_text(size = 15, face = "bold"))
  
  ggsave(file = glue("{plt_path}/valcano_{sel_ct1}_II_vs_NN_TI.png"),
         vc_plt_cd4cd103trm17,
         width = 10, height = 6,units = "in", dpi = 300, limitsize = T)
  ##
  up_gene <- deg_df_it_cd4x$Term[which(deg_df_it_cd4x$padj < 0.05 & deg_df_it_cd4x$log2FC > 0.25)]
  down_gene <- deg_df_it_cd4x$Term[which(deg_df_it_cd4x$padj < 0.05 & deg_df_it_cd4x$log2FC < -0.25)]
  if (length(up_gene) > 10) {
    Enrichment.pipline(gene_list = up_gene,
                       list_name = glue("Up_IIvsNN_{sel_ct1}"),
                       num_show = 10,
                       plot = T,
                       GO = T,
                       KEGG = T,
                       Reactome = T,
                       outpath = glue("{deg_path_it}Enrichment"))
    enrich_bar_up <- enrich_barplt(prefix = glue("Up_IIvsNN_{sel_ct1}"),
                                   df_path = glue("{deg_path_it}Enrichment/"),
                                   nn = 5, 
                                   min_gene = 10,
                                   max_nchar = 50)
    if (nrow(enrich_bar_up$enrich_df) > 0) {
      ggsave(glue("{plt_path}Enrich_IIvsNN_{sel_ct1}_up.png"),
           enrich_bar_up$bar_plt,
           height = 5, width = 7, units = "in", dpi = 300)
    }
  }
  if (length(down_gene) > 10) {
    Enrichment.pipline(gene_list = down_gene,
                       list_name = glue("Down_IIvsNN_{sel_ct1}"),
                       num_show = 10,
                       plot = T,
                       GO = T,
                       KEGG = T,
                       Reactome = T,
                       outpath = glue("{deg_path_it}Enrichment"))
    
    enrich_bar_down <- enrich_barplt(prefix = glue("Down_IIvsNN_{sel_ct1}"),
                                     df_path = glue("{deg_path_it}Enrichment/"),
                                     nn = 5, 
                                     min_gene = 10,
                                     max_nchar = 50)
    if (nrow(enrich_bar_down$enrich_df) > 0) {
      ggsave(glue("{plt_path}Enrich_IIvsNN_{sel_ct1}_down.png"),
             enrich_bar_down$bar_plt,
             height = 5, width = 7, units = "in", dpi = 300)
    }
    
  }
  return(sel_ct1)
})

##### CCC #####
ccc_path_ti <- "03_output/05_Interaction/multinichenetr/TI/CD4T/DESeq2/ann_level5_final/"
ccc_path_cr <- "03_output/05_Interaction/multinichenetr/Colon/CD4T/DESeq2/ann_level5_final/"
system(glue("cp -f {ccc_path_ti}circos_CD4T_to_CD4_CD103_TRM_Th17_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4T_to_CD4_CD103_TRM_Th17_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_CD4T_to_CD4_CD103_TRM_Th17_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4T_to_CD4_CD103_TRM_Th17_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4T_to_CD4_CD103_TRM_Th17_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4T_to_CD4_CD103_TRM_Th17_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4T_to_CD4_CD103_TRM_Th17_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4T_to_CD4_CD103_TRM_Th17_NU_vs_NN_Colon_top50.pdf"))

ccc_path_ti <- "03_output/05_Interaction/multinichenetr/TI/CD4_CD103_TRM_Th17/DESeq2/ann_level5_final/"
ccc_path_cr <- "03_output/05_Interaction/multinichenetr/Colon/CD4_CD103_TRM_Th17/DESeq2/ann_level5_final/"
system(glue("cp -f {ccc_path_ti}circos_CD4_CD103_TRM_Th17_to_CD4_CD103_TRM_Th17_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4_CD103_TRM_Th17_to_CD4_CD103_TRM_Th17_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_CD4_CD103_TRM_Th17_to_CD4_CD103_TRM_Th17_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4_CD103_TRM_Th17_to_CD4_CD103_TRM_Th17_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4_CD103_TRM_Th17_to_CD4_CD103_TRM_Th17_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4_CD103_TRM_Th17_to_CD4_CD103_TRM_Th17_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4_CD103_TRM_Th17_to_CD4_CD103_TRM_Th17_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4_CD103_TRM_Th17_to_CD4_CD103_TRM_Th17_NU_vs_NN_Colon_top50.pdf"))

system(glue("cp -f {ccc_path_ti}circos_CD4_CD103_TRM_Th17_to_CD4T_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4_CD103_TRM_Th17_to_CD4T_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_CD4_CD103_TRM_Th17_to_CD4T_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4_CD103_TRM_Th17_to_CD4T_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4_CD103_TRM_Th17_to_CD4T_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4_CD103_TRM_Th17_to_CD4T_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4_CD103_TRM_Th17_to_CD4T_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4_CD103_TRM_Th17_to_CD4T_NU_vs_NN_Colon_top50.pdf"))


# ##### test #####
# contrax <- "II_vs_NN"
# ct <- "CD4T"
# de_method <- "DESeq2"
# ann_col <- "ann_level4_final"
# 
# deg_path_it <- glue("03_output/04_Diff/DEG/TI/{ct}/{de_method}/")
# dar_path_it <- glue("03_output/04_Diff/DAR/TI/{ct}/{de_method}/")
# deg_path_cr <- glue("03_output/04_Diff/DEG/Colon/{ct}/{de_method}/")
# dar_path_cr <- glue("03_output/04_Diff/DAR/Colon/{ct}/{de_method}/")
# ##
# deg_list <- readRDS(glue("{deg_path_it}pseudo_bulk_{ann_col}_{contrax}.rds"))
# dar_list <- readRDS(glue("{dar_path_it}pseudo_bulk_{ann_col}_{contrax}.rds"))
# 
# e_regulon_filter <- readRDS("03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/e_regulon_filter_rho_add.rds")
# e_regulon_filter$eRegulon_simplify <- str_split_i(e_regulon_filter$Gene_signature_name_simplify, "_", 1)
# enrich_grn <- e_regulon_filter[,c("eRegulon_simplify", "Gene")]
# colnames(enrich_grn) <- c("term", "gene")
# ##
# ctx <- "CD4_CD103_TRM_Th17"
# deg_df_ctx <- deg_list[[ctx]]
# deg_df_ctx <- subset(deg_df_ctx, !is.na(deg_df_ctx$padj))
# deg_df_ctx$Sig <- deg_df_ctx$padj < 0.05 & (deg_df_ctx$log2FC) < -0.25
# #
# library(clusterProfiler)
# enrich_gene <- enricher(
#   deg_df_ctx$Term[deg_df_ctx$Sig],
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   # universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.05,
#   TERM2GENE = enrich_grn,
#   TERM2NAME = NA
# )
# enrich_dot <- dotplot(enrich_gene, 
#                       showCategory = 10) + 
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) + 
#   theme(axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 13, color = "black"),
#         legend.title = element_text(size = 13, face = "bold"),
#         legend.text = element_text(size = 10, color = "black"))
# ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_down.png"),
#        enrich_dot,
#        height = 4, width = 5.5, units = "in", dpi = 300)
# saveRDS(enrich_gene@result, file = glue("{deg_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}.rds"))
# ##
# deg_df_ctx$Sig <- deg_df_ctx$padj < 0.05 & (deg_df_ctx$log2FC) > 0.25
# #
# library(clusterProfiler)
# enrich_gene <- enricher(
#   deg_df_ctx$Term[deg_df_ctx$Sig],
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   # universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.05,
#   TERM2GENE = enrich_grn,
#   TERM2NAME = NA
# )
# enrich_dot <- dotplot(enrich_gene, 
#                       showCategory = 10) + 
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) + 
#   theme(axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 13, color = "black"),
#         legend.title = element_text(size = 13, face = "bold"),
#         legend.text = element_text(size = 10, color = "black"))
# ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_up.png"),
#        enrich_dot,
#        height = 4, width = 5.5, units = "in", dpi = 300)
# saveRDS(enrich_gene@result, file = glue("{deg_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}_up.rds"))
# 
# ##
# ctx <- "CD4_CD103_TRM_Th17"
# enrich_grnr <- e_regulon_filter[,c("eRegulon_simplify", "Region")]
# colnames(enrich_grnr) <- c("term", "gene")
# enrich_grnr$gene <- gsub("\\:", "-", enrich_grnr$gene)
# 
# dar_df_ctx <- dar_list[[ctx]]
# dar_df_ctx <- subset(dar_df_ctx, !is.na(dar_df_ctx$padj))
# dar_df_ctx$Sig <- dar_df_ctx$padj < 0.05 & (dar_df_ctx$log2FC) < -0.25
# #
# enrich_gene <- enricher(
#   dar_df_ctx$Term[dar_df_ctx$Sig],
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   # universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.05,
#   TERM2GENE = enrich_grnr,
#   TERM2NAME = NA
# )
# enrich_dot <- dotplot(enrich_gene, 
#                       showCategory = 10) + 
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) + 
#   theme(axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 13, color = "black"),
#         legend.title = element_text(size = 13, face = "bold"),
#         legend.text = element_text(size = 10, color = "black"))
# ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_dar_down.png"),
#        enrich_dot,
#        height = 4, width = 5.5, units = "in", dpi = 300)
# # saveRDS(enrich_gene@result, file = glue("{deg_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}.rds"))
# ##
# dar_df_ctx$Sig <- dar_df_ctx$padj < 0.05 & (dar_df_ctx$log2FC) > 0.25
# #
# library(clusterProfiler)
# enrich_gene <- enricher(
#   dar_df_ctx$Term[dar_df_ctx$Sig],
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   # universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.05,
#   TERM2GENE = enrich_grnr,
#   TERM2NAME = NA
# )
# enrich_dot <- dotplot(enrich_gene, 
#                       showCategory = 10) + 
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) + 
#   theme(axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 13, color = "black"),
#         legend.title = element_text(size = 13, face = "bold"),
#         legend.text = element_text(size = 10, color = "black"))
# ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_dar_up.png"),
#        enrich_dot,
#        height = 4, width = 5.5, units = "in", dpi = 300)
# # saveRDS(enrich_gene@result, file = glue("{deg_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}_up.rds"))
# 
# ##
# 
# contrax <- "II_vs_NN"
# ct <- "CD4_CD103_TRM_Th17"
# de_method <- "DESeq2"
# ann_col <- "ann_level5_final"
# 
# deg_path_it <- glue("03_output/04_Diff/DEG/TI/{ct}/{de_method}/")
# dar_path_it <- glue("03_output/04_Diff/DAR/TI/{ct}/{de_method}/")
# deg_path_cr <- glue("03_output/04_Diff/DEG/Colon/{ct}/{de_method}/")
# dar_path_cr <- glue("03_output/04_Diff/DAR/Colon/{ct}/{de_method}/")
# ##
# deg_list <- readRDS(glue("{deg_path_it}pseudo_bulk_{ann_col}_{contrax}.rds"))
# dar_list <- readRDS(glue("{dar_path_it}pseudo_bulk_{ann_col}_{contrax}.rds"))
# 
# enrich_grn_region <- e_regulon_filter[,c("eRegulon_simplify", "Region")]
# enrich_grn_region$Region <- gsub("\\:", "-", enrich_grn_region$Region)
# colnames(enrich_grn_region) <- c("term", "gene")
# 
# enrich_grn <- e_regulon_filter[,c("eRegulon_simplify", "Gene")]
# colnames(enrich_grn) <- c("term", "gene")
# 
# #
# ctx <- "CD4_CD103_TRM_Th17_HSP"
# deg_df_ctx <- deg_list[[ctx]]
# deg_df_ctx <- subset(deg_df_ctx, !is.na(deg_df_ctx$padj))
# deg_df_ctx$Sig <- deg_df_ctx$padj < 0.05 & deg_df_ctx$log2FC < -0.25
# #
# library(clusterProfiler)
# enrich_region <- enricher(
#   deg_df_ctx$Term[deg_df_ctx$Sig],
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   # universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.05,
#   TERM2GENE = enrich_grn,
#   TERM2NAME = NA
# )
# enrich_dot <- dotplot(enrich_region, 
#                       showCategory = 10) + 
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) + 
#   theme(axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 13, color = "black"),
#         legend.title = element_text(size = 13, face = "bold"),
#         legend.text = element_text(size = 10, color = "black"))
# ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_down.png"),
#        enrich_dot,
#        height = 4, width = 5.5, units = "in", dpi = 300)
# # saveRDS(enrich_region@result, file = glue("{dar_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}.rds"))
# ##
# #
# ctx <- "CD4_CD103_TRM_Th17_HSP"
# dar_df_ctx <- dar_list[[ctx]]
# dar_df_ctx <- subset(dar_df_ctx, !is.na(dar_df_ctx$padj))
# dar_df_ctx$Sig <- dar_df_ctx$padj < 0.05 & dar_df_ctx$log2FC < -0.25
# #
# library(clusterProfiler)
# enrich_region <- enricher(
#   dar_df_ctx$Term[dar_df_ctx$Sig],
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   # universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.05,
#   TERM2GENE = enrich_grn_region,
#   TERM2NAME = NA
# )
# enrich_dot <- dotplot(enrich_region, 
#                       showCategory = 10) + 
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) + 
#   theme(axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 13, color = "black"),
#         legend.title = element_text(size = 13, face = "bold"),
#         legend.text = element_text(size = 10, color = "black"))
# ggsave(glue("{plt_path}enrich_dot_{ann_col}_{contrax}_{ctx}_down_dar.png"),
#        enrich_dot,
#        height = 4, width = 5.5, units = "in", dpi = 300)
# # saveRDS(enrich_region@result, file = glue("{dar_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}.rds"))
# 
# ##
# contrax <- "I_I_vs_N_N"
# ct <- "CD4_CD103_TRM"
# de_method <- "DESeq2"
# ann_col <- "ann_level5_final"
# 
# deg_path_it <- glue("03_output/04_Diff/DEG/small_intestine/{ct}/{de_method}/")
# dar_path_it <- glue("03_output/04_Diff/DAR/small_intestine/{ct}/{de_method}/")
# deg_path_cr <- glue("03_output/04_Diff/DEG/colorect/{ct}/{de_method}/")
# dar_path_cr <- glue("03_output/04_Diff/DAR/colorect/{ct}/{de_method}/")
# ##
# deg_list <- readRDS(glue("{deg_path_it}pseudo_bulk_{ann_col}_{contrax}.rds"))
# dar_list <- readRDS(glue("{dar_path_it}pseudo_bulk_{ann_col}_{contrax}.rds"))
# 
# e_regulon_filter <- readRDS("03_output/07_SCENIC/test/CD4T/scplus_pipeline/Snakemake/e_regulon_filter_rho_add.rds")
# e_regulon_filter$eRegulon_simplify <- str_split_i(e_regulon_filter$Gene_signature_name_simplify, "_", 1)
# enrich_grn <- e_regulon_filter[,c("eRegulon_simplify", "Gene")]
# colnames(enrich_grn) <- c("term", "gene")
# ##
# ctx <- "CD4_CD103_TRM_Th17_MAIT_like"
# deg_df_ctx <- deg_list[[ctx]]
# deg_df_ctx <- subset(deg_df_ctx, !is.na(deg_df_ctx$padj))
# deg_df_ctx$Sig <- deg_df_ctx$padj < 0.05 & abs(deg_df_ctx$log2FC) > 0.25
# #
# library(clusterProfiler)
# enrich_gene <- enricher(
#   deg_df_ctx$Term[deg_df_ctx$Sig],
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   # universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.05,
#   TERM2GENE = enrich_grn,
#   TERM2NAME = NA
# )
# enrich_dot <- dotplot(enrich_gene, 
#                       showCategory = 10) + 
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) + 
#   theme(axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 13, color = "black"),
#         legend.title = element_text(size = 13, face = "bold"),
#         legend.text = element_text(size = 10, color = "black"))
# ggsave(glue("{deg_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}.png"),
#        enrich_dot,
#        height = 4, width = 5.5, units = "in", dpi = 300)
# saveRDS(enrich_gene@result, file = glue("{deg_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}.rds"))
# ##
# enrich_grn_region <- e_regulon_filter[,c("eRegulon_simplify", "Region")]
# enrich_grn_region$Region <- gsub("\\:", "-", enrich_grn_region$Region)
# colnames(enrich_grn_region) <- c("term", "gene")
# #
# ctx <- "CD4_CD103_TRM_Th17_CREM"
# dar_df_ctx <- dar_list[[ctx]]
# dar_df_ctx <- subset(dar_df_ctx, !is.na(dar_df_ctx$padj))
# dar_df_ctx$Sig <- dar_df_ctx$padj < 0.05 & abs(dar_df_ctx$log2FC) > 0.25
# #
# library(clusterProfiler)
# enrich_region <- enricher(
#   dar_df_ctx$Term[dar_df_ctx$Sig],
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   # universe,
#   minGSSize = 10,
#   maxGSSize = 500,
#   qvalueCutoff = 0.05,
#   TERM2GENE = enrich_grn_region,
#   TERM2NAME = NA
# )
# enrich_dot <- dotplot(enrich_region, 
#                       showCategory = 10) + 
#   scale_fill_distiller(palette = "RdYlBu", direction = 1) + 
#   theme(axis.title = element_text(size = 14, face = "bold"),
#         axis.text = element_text(size = 13, color = "black"),
#         legend.title = element_text(size = 13, face = "bold"),
#         legend.text = element_text(size = 10, color = "black"))
# ggsave(glue("{dar_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}.png"),
#        enrich_dot,
#        height = 4, width = 5.5, units = "in", dpi = 300)
# saveRDS(enrich_region@result, file = glue("{dar_path_it}enrich_dot_{ann_col}_{contrax}_{ctx}.rds"))
# 
# 
# ##
# ctx <- "CD4_CD103_TRM_Th17"
# deg_cd103trm17 <- readRDS("03_output/04_Diff/DEG/small_intestine/CD4T/DESeq2/pseudo_bulk_ann_level4_final_I_I_vs_N_N.rds")[["CD4_CD103_TRM_Th17"]]
# deg_cd103trm17 <- deg_cd103trm17[!is.na(deg_cd103trm17$padj),]
# #
# links_gr <- readRDS("03_output/03_clustering_test/CD4T/test153000/test302015/LinkPeaks_links_cd103trm17.rds")
# df_links <- as.data.frame(links_gr)
# df_links$padj <- p.adjust(df_links$pvalue, method = "BH")
# dff_links <- subset(df_links, df_links$padj < 0.05 & df_links$score > 0.1)
# dff_links_sig <- subset(dff_links, peak %in% 
#                           dar_cd103trm17$Term[which(dar_cd103trm17$padj < 0.05 & abs(dar_cd103trm17$log2FC) > 0.25)])
# target_df2 <- deg_cd103trm17[which(deg_cd103trm17$padj < 0.05 &
#                                      abs(deg_cd103trm17$log2FC) > 0.25 &
#                                      deg_cd103trm17$Term %in% magma_df$Term[magma_df$P < 1E-5]),]
# target_df2$sig <- -log10(target_df2$padj)
# 
# source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/volca_plot.R")
# vc_plt <- volca_function(deg_df = deg_cd103trm17,
#                          term_col = "Term",
#                          fc_col = "log2FC",
#                          p_col = "padj",
#                          # target_gene = unique(dff_links_sig$gene),
#                          n_top = 20,
#                          log_thresh = 0.25, 
#                          top_log_thresh = 0.5,
#                          group_label = c("Down regulation", "Others", "Up regulation"),
#                          valco_col = c4a("classic_blue_red12", 3),
#                          highlight_col = "orange")
# ggsave(file = glue("{deg_path_it}/valcano_{ctx}_I_I_vs_N_N_TI.png"),
#        vc_plt,
#        width = 7, height = 5,units = "in", dpi = 300, limitsize = T)

