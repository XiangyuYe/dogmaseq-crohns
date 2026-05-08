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
ct <- "all"
cov_use = c("Age", "Gender", 
            "use_Aminosalic", "use_Immunomodu", "use_MTX",
            "use_AntiTNF", "use_AntiIL23", "use_AntiIntegrin")
block_col = "Sample_exp"
sample_col = "Sample_ID_exp"
ann_col <- "ann_level2_refine"
group_col <- "Condition"
subgroup_col <- "Section"

##
data_path <- "03_output/03_clustering/"
out_path <- "03_output/04_Diff/DAR/"
sc_meta_all <- readRDS(glue("{data_path}/all/WNN_ADT_RNA/sc_meta_ann_level2_refine.rds"))
sc_meta_all$Age <- scale(sc_meta_all$Age, center = T, scale = T) %>% as.numeric()
#
scATAC_obj <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_scATAC_obj.rds")
scATAC_obj <- subset(scATAC_obj, cells = rownames(sc_meta_all))
DefaultAssay(scATAC_obj)<-"peaks"
##
sc_meta <- sc_meta_all[, c(block_col, sample_col, ann_col, 
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
apply(comb_group, 2, function(comb_groupx){
  g1 <- comb_groupx[1]
  g2 <- comb_groupx[2]
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
  #
  return(out_prefix)
})

# ##### limma #####
# de_method <- "limma"
# for (ngx in c(2, 3, 1)) {
#   
#   g1 <- comb_group[1, ngx]
#   g2 <- comb_group[2, ngx]
#   if (g1 == "I_I" & g2 == "N_U") {
#     cov_use = c("Age", "Gender", "Anti_Inflam", "Immunomodulatory", "Biologics")
#     block_col = NULL
#   } else {
#     cov_use = c("Age", "Gender", "Anti_Inflam", "Immunomodulatory", "Biologics")
#     block_col = "sample_exp"
#   }
#   ## all
#   scATAC_obj_sub <- subset(scATAC_obj, condition %in% c(g1, g2))
#   scATAC_obj_sub$condition <- factor(scATAC_obj_sub$condition, 
#                                     levels = c(g2, g1))
#   out_prefix <- glue("pseudo_bulk_{ann_col}_{g1}_vs_{g2}.rds")
#   ## small intestine
#   deg_list_it <- pseudo.bulk.deg(scATAC_obj = subset(scATAC_obj_sub,
#                                                     section_comb == "TI"),
#                                  ident_col = ann_col,
#                                  group_col = "condition",
#                                  cov_col = cov_use,
#                                  sample_col = "sample_id_exp",
#                                  block_col = block_col,
#                                  min_cell_per_sample = 10,
#                                  filter_ByExpr = filter_ByCPM,
#                                  min_percent = 0.3,
#                                  test_use = "limma",
#                                  n_core = 10)
#   saveRDS(deg_list_it, file = glue("{out_path_it}/{de_method}/{out_prefix}"))
#   ## colorect
#   deg_list_cr <- pseudo.bulk.deg(scATAC_obj = subset(scATAC_obj_sub, 
#                                                     section_comb == "Colon"),
#                                  ident_col = ann_col,
#                                  group_col = "condition",
#                                  cov_col = cov_use,
#                                  sample_col = "sample_id_exp",
#                                  block_col = block_col,
#                                  min_cell_per_sample = 10,
#                                  filter_ByExpr = filter_ByCPM,
#                                  min_percent = 0.3,
#                                  test_use = "limma",
#                                  n_core = 10)
#   saveRDS(deg_list_cr, file = glue("{out_path_cr}/{de_method}/{out_prefix}"))
#   
# }
# 
# ##### edgeR #####
# de_method <- "edgeR"
# for (ngx in 1:ncol(comb_group)) {
#   
#   g1 <- comb_group[1, ngx]
#   g2 <- comb_group[2, ngx]
#     cov_use = c("Age", "Gender", "Anti_Inflam", "Immunomodulatory", "Biologics")
#     block_col = "sample_exp"
#   ## all
#   scATAC_obj_sub <- subset(scATAC_obj, condition %in% c(g1, g2))
#   scATAC_obj_sub$condition <- factor(scATAC_obj_sub$condition, 
#                                     levels = c(g2, g1))
#   out_prefix <- glue("pseudo_bulk_{ann_col}_{g1}_vs_{g2}.rds")
#   ## small intestine
#   deg_list_it <- pseudo.bulk.deg(scATAC_obj = subset(scATAC_obj_sub,
#                                                     section_comb == "TI"),
#                                  ident_col = ann_col,
#                                  group_col = "condition",
#                                  cov_col = cov_use,
#                                  sample_col = "sample_id_exp",
#                                  block_col = block_col,
#                                  min_cell_per_sample = 10,
#                                  filter_ByExpr = filter_ByCPM,
#                                  min_percent = 0.3,
#                                  test_use = "edgeR-QLF",
#                                  n_core = 10)
#   saveRDS(deg_list_it, file = glue("{out_path_it}/{de_method}/{out_prefix}"))
#   ## colorect
#   deg_list_cr <- pseudo.bulk.deg(scATAC_obj = subset(scATAC_obj_sub, 
#                                                     section_comb == "Colon"),
#                                  ident_col = ann_col,
#                                  group_col = "condition",
#                                  cov_col = cov_use,
#                                  sample_col = "sample_id_exp",
#                                  block_col = block_col,
#                                  min_cell_per_sample = 10,
#                                  filter_ByExpr = filter_ByCPM,
#                                  min_percent = 0.3,
#                                  test_use = "edgeR-QLF",
#                                  n_core = 10)
#   saveRDS(deg_list_cr, file = glue("{out_path_cr}/{de_method}/{out_prefix}"))
#   
# }
# ##### MAST-lmm #####
# source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/sc_DEG.R")
# Idents(scATAC_obj) <- scATAC_obj@meta.data$condition
# scATAC_obj <- NormalizeData(scATAC_obj)
# de_method <- "MAST"
# for (ngx in c(2, 3, 1)) {
#   
#   g1 <- comb_group[1, ngx]
#   g2 <- comb_group[2, ngx]
#   if (g1 == "I_I" & g2 == "N_U") {
#     cov_use = c("Age", "Gender", "Anti_Inflam", "Immunomodulatory", "Biologics")
#     block_col = NULL
#   } else {
#     cov_use = c("Age", "Gender", "Anti_Inflam", "Immunomodulatory", "Biologics")
#     block_col = "sample_exp"
#   }
#   scATAC_obj_sub <- subset(scATAC_obj, condition %in% c(g1, g2))
#   scATAC_obj_sub$condition <- factor(scATAC_obj_sub$condition, 
#                                     levels = c(g2, g1))
#   out_prefix <- glue("lmm_{ann_col}_{g1}_vs_{g2}.rds")
#   deg_list_it <- sc.deg(seurat_obj = subset(scATAC_obj_sub, 
#                                             section_comb == "TI"),
#                         use_assay = "RNA",
#                         ident_col = ann_col,
#                         group_col = "condition",
#                         cov_col = cov_use,
#                         sample_col = "sample_id_exp",
#                         block_col = block_col,
#                         min_percent = 0.05,
#                         test_use = "MAST-lmm",
#                         do_para = T)
#   
#   saveRDS(deg_list_it, file = glue("{out_path_it}/{de_method}/{out_prefix}"))
#   deg_list_cr <- sc.deg(seurat_obj = subset(scATAC_obj_sub, 
#                                             section_comb == "Colon"),
#                         use_assay = "RNA",
#                         ident_col = ann_col,
#                         group_col = "condition",
#                         cov_col = cov_use,
#                         sample_col = "sample_id_exp",
#                         block_col = block_col,
#                         min_percent = 0.05,
#                         test_use = "MAST-lmm",
#                         do_para = T)
#   saveRDS(deg_list_cr, file = glue("{out_path_cr}/{de_method}/{out_prefix}"))
# }

# 
# source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/volca_plot.R")
# deg_df_it <- readRDS(glue("{out_path}small_intestine/pseudo_bulk_CD4_I_I_vs_N_N_limma.rds"))[[1]]
# # Term    log2FC     stat       pvalue         padj cluster
# vc_plt <- volca_function(deg_df = deg_df_it,
#                          term_col = "Term",
#                          fc_col = "log2FC",
#                          p_col = "padj",
#                          # target_gene = NULL,
#                          n_top = 20,
#                          log_thresh = 0.5, 
#                          group_label = c("Down regulation", "Others", "Up regulation"),
#                          valco_col = c4a("classic_blue_red12", 3),
#                          highlight_col = "orange")
# tiff(file = glue("{out_path}small_intestine/valcano_CD4_I_I_vs_N_N_limma.tiff"),
#      width = 8, height = 9, units = "in", res = 600, compression = "lzw")
# vc_plt
# dev.off()
# 
# ##
# source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/Enrichment.R")
# enrich_file_path <- glue("{out_path}small_intestine/Enrichment/")
# deg_gene_list <- list("CD4T_up" = deg_df_it$Term[deg_df_it$log2FC > 0.5 & deg_df_it$padj < 0.05],
#                       "CD4T_down" = deg_df_it$Term[deg_df_it$log2FC < -0.5 & deg_df_it$padj < 0.05])
# lapply(seq_along(deg_gene_list), function(x){
#   Enrichment.pipline(gene_list = deg_gene_list[[x]],
#                      list_name = paste0("Enrichment_", names(deg_gene_list)[x]),
#                      num_show = 10,
#                      plot = T,
#                      GO = T,
#                      KEGG = T,
#                      Reactome = T,
#                      outpath = enrich_file_path)
# })
# 
# ##
# deg_df_cr <- readRDS(glue("{out_path}colorect/pseudo_bulk_CD4_I_I_vs_N_N_limma.rds"))[[1]]
# # Term    log2FC     stat       pvalue         padj cluster
# vc_plt <- volca_function(deg_df = deg_df_cr,
#                          term_col = "Term",
#                          fc_col = "log2FC",
#                          p_col = "padj",
#                          # target_gene = NULL,
#                          n_top = 20,
#                          log_thresh = 0.5, 
#                          group_label = c("Down regulation", "Others", "Up regulation"),
#                          valco_col = c4a("classic_blue_red12", 3),
#                          highlight_col = "orange")
# tiff(file = glue("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/04_Diff/DEG/colorect/valcano_CD4_I_I_vs_N_N_limma.tiff"),
#      width = 8, height = 9, units = "in", res = 600, compression = "lzw")
# vc_plt
# dev.off()
# ##
# enrich_file_path <- glue("03_output/04_Diff/DEG/colorect/Enrichment/")
# deg_gene_list <- list("CD4T_up" = deg_df_cr$Term[deg_df_cr$log2FC > 0.5 & deg_df_cr$padj < 0.05],
#                       "CD4T_down" = deg_df_cr$Term[deg_df_cr$log2FC < -0.5 & deg_df_cr$padj < 0.05])
# lapply(seq_along(deg_gene_list), function(x){
#   Enrichment.pipline(gene_list = deg_gene_list[[x]],
#                      list_name = paste0("Enrichment_", names(deg_gene_list)[x]),
#                      num_show = 10,
#                      plot = T,
#                      GO = T,
#                      KEGG = T,
#                      Reactome = T,
#                      outpath = enrich_file_path)
# })
# 
# ##
# library(bigreadr)
# library(dplyr)
# library(stringr)
# library(glue)
# library(tibble)
# library(ggplot2)
# library(ggrepel)
# 
# setwd("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/")
# ct <- "CD4T"
# tissuex <- "TI"
# deg_df_it <- deg_list_it[[1]]
# deg_df_cr <- deg_list_cr[[1]]
# 
# deg_df_it$Expression <- ifelse(deg_df_it$padj >= 0.05 | abs(deg_df_it$log2FC) < 0.5,
#                                "Not Significant",
#                                ifelse(deg_df_it$log2FC > 0, "Up-regulated", "Down-regulated"))
# deg_df_cr$Expression <- ifelse(deg_df_cr$padj >= 0.05 | abs(deg_df_cr$log2FC) < 0.5,
#                                "Not Significant",
#                                ifelse(deg_df_cr$log2FC > 0, "Up-regulated", "Down-regulated"))
# 
# ##
# deg_df_it <- deg_df_it[, c("Term", "log2FC", "padj", "Expression")]
# deg_df_cr <- deg_df_cr[, c("Term", "log2FC", "padj", "Expression")]
# 
# colnames(deg_df_it) <- paste0(colnames(deg_df_it), "_it")
# colnames(deg_df_cr) <- paste0(colnames(deg_df_cr), "_cr")
# 
# deg_comb <- merge(deg_df_it, 
#                   deg_df_cr, 
#                   by.x = "Term_it", 
#                   by.y = "Term_cr")
# deg_comb$Group <- ifelse(deg_comb$Expression_it == "Not Significant" & 
#                            deg_comb$Expression_cr == "Not Significant",
#                          "Neither",
#                          ifelse(deg_comb$Expression_it != "Not Significant" & 
#                                   deg_comb$Expression_cr != "Not Significant",
#                                 "Both",
#                                 ifelse(deg_comb$Expression_it != "Not Significant",
#                                        "Significant in Ileum only",
#                                        "Significant in Colon only"))) %>% 
#   factor(., levels = c("Both", "Significant in Ileum only", "Significant in Colon only", "Neither"))
# 
# # deg_comb <- subset(deg_comb, Group != "Neither")
# top_df <- top_n(deg_comb[deg_comb$Group == "Both",], 50, 
#                 abs(log2FC_it + log2FC_cr))
# 
# max_X <- max(abs(deg_comb$log2FC_cr)) %>% ceiling()
# max_Y <- max(abs(deg_comb$log2FC_it)) %>% ceiling()
# use_color <- c("#F02720", "#E9C39B", "#B5C8E2", "grey")
# plt_inter <- ggplot(deg_comb) + 
#   geom_point(aes(x = log2FC_cr, y = log2FC_it, color = Group),
#              shape = 19, size = 1.5) + 
#   geom_hline(yintercept = 0, linetype = 2, size = 1, color = "grey20")+ 
#   geom_vline(xintercept = 0, linetype = 2, size = 1, color = "grey20")+ 
#   scale_x_continuous(limits = c(-max_X, max_X)) + 
#   scale_y_continuous(limits = c(-max_Y, max_Y)) + 
#   xlab(bquote(""~Log[2]~"(Fold Change) in Colon"))+
#   ylab(bquote(""~Log[2]~"(Fold Change) in Ileum"))+
#   geom_label_repel(data = top_df, 
#                    aes(label = Term_it, x = log2FC_cr, y = log2FC_it), 
#                    segment.size = 0.25, size = 2.5, 
#                    max.iter = 1E7, min.segment.length = 2,
#                    max.overlaps = getOption("ggrepel.max.overlaps", default = 100))+
#   scale_colour_manual(values = use_color) + 
#   theme_bw()+
#   theme(legend.title = element_blank(),
#         legend.position = "top",
#         legend.text = element_text(size = 12, color = "black"),
#         axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 15, face= "bold", color = "black"))
# 
# tiff(glue("03_output/04_Diff/DEG/valcano_comb_CD4.tiff"),
#      height = 8, width = 7, units = "in", compression = "lzw", res = 300)
# plt_inter
# dev.off()
