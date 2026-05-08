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
plt_path = "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/04_plot/Fig3/"
pair_list <- list("II_vs_NN_TI" = c("II:TI", "NN:TI"),
                  "NU_vs_NN_TI" = c("NU:TI", "NN:TI"),
                  "II_vs_NN_Colon" = c("II:Colon", "NN:Colon"),
                  "NU_vs_NN_Colon" = c("NU:Colon", "NN:Colon"),
                  "TI_vs_Colon_NN" = c("NN:TI", "NN:Colon"))
# ct_cols <- c("#00897B", "#06C3F2", "#016FAD", "#004F78", "#00A0E9",
#              "#FFEB3B", "#FFB300", "#F57C33", "#FF5722", "#D32F2F",  "#A14141",  "#5D4037",
#              "#F8AAC0", "#D81B60", "#F06292")
ct_cols <- c("#F8AAC0",  "#D81B60",  "#F06292")
gp_cols <- c("#0073C2FF","#EFC000FF", "#E41A1C")

seed_use <- 20250528
##### load data #####
sc_obj <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/scWNN_obj.rds")
sc_meta_CD4 <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta.rds")
sc_obj@meta.data <- sc_meta_CD4
sc_obj$Condition <- factor(sc_obj$Condition, levels = c("NN", "NU", "II"))
sc_obj$Section <- factor(sc_obj$Section, levels = c("Colon", "TI"))
sc_obj_treg <- subset(sc_obj, ann_level3_final == "CD4_Treg")

## load SCENIC+ output (for 3A and B)
sp_outpath <- "03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/"
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))

##### Fig. 3A: RSS of Tregs #####
## load rss data
rss_CD4 <- fread2(glue("{sp_outpath}rss_CD4_ann_level4_final.txt")) %>%
  tibble::column_to_rownames(., "V1")
colnames(rss_CD4) <- e_regulon_name_trans$signature_name_simplify[match(colnames(rss_CD4), 
                                                                        e_regulon_name_trans$signature)]
rss_CD4$celltype <- rownames(rss_CD4)
## reshape and set top 5 eRegulons for each celltype
rss_CD4_melt <- melt(rss_CD4,
                     id.vars = "celltype",
                     variable.name = "eRegulon",
                     value.name = "RSS")
rss_CD4_melt <- split(rss_CD4_melt, rss_CD4_melt$celltype) %>%
  lapply(., function(x){
    x$Rank <- rank(-x$RSS)
    x$Group <- ifelse(x$Rank <= 10,
                      "Top", "Others")
    return(x)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
rss_CD4_melt$Group <- factor(rss_CD4_melt$Group, levels = c("Others", "Top"))
rss_CD4_melt$celltype <- factor(rss_CD4_melt$celltype,
                                levels = levels(sc_meta_CD4$ann_level4_final),
                                labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))
rss_CD4_melt$Label <- rss_CD4_melt$eRegulon
rss_CD4_melt$Label[rss_CD4_melt$Group == "Others"] <- NA
## subset Treg
rss_treg_melt <- subset(rss_CD4_melt, grepl("Treg", rss_CD4_melt$celltype))
rss_treg_melt$celltype <- droplevels(rss_treg_melt$celltype)
plt_rss <- ggplot(data = rss_treg_melt) + 
  geom_point(aes(x = Rank, y = RSS, color = Group, size = Group)) + 
  scale_size_manual(values = c(1, 2)) +
  scale_color_manual(values = c("#4682B4", "red")) + 
  facet_wrap(~rss_treg_melt$celltype, nrow = 1, scale = "free_y") + 
  geom_text_repel(aes(x = Rank, y = RSS, label = Label),
                  color = "red", size = 3,
                  max.iter = 1E7, 
                  max.overlaps = 2000) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, color = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        strip.background = element_blank())
#
ggsave(glue("{plt_path}rss_Treg_ann_level4_final.png"),
       plt_rss,
       height = 5, width = 15, units = "in", dpi = 300)

####### Fig. 3B: ct-specific eRegulon #########
## load AUC object
scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
scAUC_obj_treg <- subset(scAUC_obj, cells = colnames(sc_obj_treg))
scAUC_obj_treg@meta.data <- sc_obj_treg@meta.data
Idents(scAUC_obj_treg) <- scAUC_obj_treg$ann_level4_final
er_gene <- rownames(scAUC_obj_treg)[grep("g\\)", rownames(scAUC_obj_treg))]
## eRegulon markers
diff_auc <- FindAllMarkers(scAUC_obj_treg, 
                           features = er_gene,
                           only.pos = T)
saveRDS(diff_auc, "03_output/03_clustering/CD4T/WNN_ADT_RNA/diff_auc_treg.rds")
#
sig_diff_auc <- diff_auc %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  top_n(20, -p_val_adj) %>%
  ungroup()
sig_diff_auc <- rbind(sig_diff_auc, diff_auc[grep("RORC|CREM|IKZF2", diff_auc$gene),])
sig_diff_auc <- sig_diff_auc[order(sig_diff_auc$p_val_adj, -abs(sig_diff_auc$avg_log2FC)),]
sig_diff_auc <- sig_diff_auc[!duplicated(sig_diff_auc$gene),]
sig_diff_auc <- sig_diff_auc[order(sig_diff_auc$cluster, sig_diff_auc$p_val_adj, -abs(sig_diff_auc$avg_log2FC)),]
sig_auc_marker <- str_split_i(sig_diff_auc$gene, "-\\(", 1)

## load object for heatdot plot
heatdot_treg <- fread2(glue("{sp_outpath}heatdot_Treg_ann_level4_final.txt"))[,-1]
colnames(heatdot_treg) <- c("celltype", "Gene-based AUC", "eRegulon_name", "Region-based AUC")
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))
e_regulon_name_trans$heatdot_name <- stringr::str_split_i(e_regulon_name_trans$signature, "_\\(", 1)
heatdot_treg$signature_name_simplify <- 
  e_regulon_name_trans$signature_name_simplify[match(heatdot_treg$eRegulon_name,
                                                     e_regulon_name_trans$heatdot_name)] %>%
  str_split_i(., "_", 1)

## intersect and format
sig_auc_marker_use <- sig_auc_marker[which(sig_auc_marker %in% heatdot_treg$signature_name_simplify)]
heatdot_treg_sub <- subset(heatdot_treg, signature_name_simplify %in% sig_auc_marker_use)
heatdot_treg_sub$signature_name_simplify <- factor(heatdot_treg_sub$signature_name_simplify,
                                                   levels = rev(sig_auc_marker_use))
heatdot_treg_sub$celltype <- factor(heatdot_treg_sub$celltype,
                                    levels = levels(scAUC_obj_treg$ann_level4_final),
                                    labels = gsub("_", " ", levels(scAUC_obj_treg$ann_level4_final)))
## plot
heat_dot <- ggplot(heatdot_treg_sub, 
                   aes(celltype, signature_name_simplify))+ 
  geom_tile(mapping = aes(fill = `Gene-based AUC`)) + 
  scale_fill_distiller(type = 'div', palette = 'RdYlBu') + 
  geom_point(mapping = aes(size = `Region-based AUC`),
             colour = "black") + 
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        legend.position = "left",
        legend.justification = c("left", "bottom"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.text.x = element_blank(),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10, color = "black"))
#
df_annot <- data.frame(celltype = levels(scAUC_obj_treg$ann_level4_final))
p_annot = ggplot(df_annot, aes(x = celltype, y = 1, fill = celltype)) +
  scale_fill_manual(values = ct_cols) + 
  geom_tile() +
  theme_void() +
  theme(legend.position = "none")

ggsave(glue("{plt_path}heat_dot_Treg_ann_level4_final.png"),
       wrap_plots(list(p_annot, heat_dot), ncol = 1, heights = c(0.1, 8)),
       height = 11, width = 5, units = "in", dpi = 300)

##### Fig. 3C: bulk heatmap for celltype marker genes #####
## DE analysis across cell types
DefaultAssay(sc_obj_treg) <- "RNA"
sc_obj_treg <- NormalizeData(sc_obj_treg) %>% ScaleData()
Idents(sc_obj_treg) <- sc_obj_treg$ann_level4_final
ct_markers <- FindAllMarkers(sc_obj_treg, 
                             logfc.threshold = 1,
                             only.pos = T)
ct_markers <- ct_markers[order(ct_markers$cluster),]
saveRDS(ct_markers, "03_output/03_clustering/CD4T/WNN_ADT_RNA/ct_markers_treg.rds")
#
top_marker <- ct_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup()
top_marker <- unique(top_marker$gene)

## format pseudo matrix
sc_mat_treg_bulk <- Seurat::AverageExpression(sc_obj_treg,
                                              features = top_marker,
                                              group.by = "ann_level4_final",
                                              assays = "RNA",
                                              return.seurat = F,
                                              layer = "scale.data")[["RNA"]] %>%
  as.data.frame()
#
ann_level <- gsub("\\-", " ", colnames(sc_mat_treg_bulk))
ann_col <- data.frame(Celltype = factor(ann_level,
                                        levels = ann_level),
                      row.names = colnames(sc_mat_treg_bulk))
## plot
anno_colors <- ct_cols
names(anno_colors) <- ann_level
png(glue("{plt_path}pseudo_heatRNA_treg.png"),
    height = 5.5, width = 3.8, units = "in", res = 300)
pheatmap::pheatmap(sc_mat_treg_bulk[top_marker,],
                   annotation_col = ann_col,
                   annotation_colors = list(Celltype = anno_colors),
                   scale = "none",
                   border_color = NA,
                   fontsize = 8,
                   show_colnames = F,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()

##### Fig. 3D: bulk heatmap for GSVA #####
## load GSVA results
gsva_obj <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/gsva_obj.rds")
gsva_obj_treg <- subset(gsva_obj, grepl("Treg", gsva_obj$ann_level4_final))
## differnetial analysis
DefaultAssay(gsva_obj_treg) <- "GSVA"
Idents(gsva_obj_treg) <- gsva_obj_treg$ann_level4_final
diff_gsva <- FindAllMarkers(gsva_obj_treg, 
                            features = rownames(gsva_obj_treg)[grep("^KEGG|^GOBP|^REACTOME", rownames(gsva_obj_treg))],
                            only.pos = T)
saveRDS(diff_gsva, "03_output/03_clustering/CD4T/WNN_ADT_RNA/diff_gsva_treg.rds")
#
top_df <- diff_gsva %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup()
top_df <- top_df[order(top_df$p_val_adj, -abs(top_df$avg_log2FC)),]
top_df <- top_df[!duplicated(top_df$gene),]
top_df <- top_df[order(top_df$cluster, top_df$p_val_adj, -abs(top_df$avg_log2FC)),]

## format pseudo matrix
gsva_mat_treg_bulk <- Seurat::AverageExpression(gsva_obj_treg,
                                                  features = unique(top_df$gene),
                                                  group.by = "ann_level4_final",
                                                  assays = "GSVA",
                                                  return.seurat = F,
                                                  layer = "data")[["GSVA"]] %>%
  as.data.frame()
#
gsva_mat_treg_bulk <- gsva_mat_treg_bulk[top_df$gene, ]
rownames(gsva_mat_treg_bulk) <- gsub("\\-", " ", rownames(gsva_mat_treg_bulk)) %>%
  gsub("KEGG", "KEGG:", .) %>%
  gsub("REACTOME", "REACTOME:", .) %>%
  gsub("GOBP", "GOBP:", .)
ann_level <- gsub("\\-", " ", colnames(gsva_mat_treg_bulk))
ann_col <- data.frame(Celltype = factor(ann_level,
                                        levels = ann_level),
                      row.names = colnames(gsva_mat_treg_bulk))
## plot
anno_colors <- ct_cols
names(anno_colors) <- ann_level
#
png(glue("{plt_path}gsva_heat_treg.png"),
    height = 5, width = 6, units = "in", res = 300)
pheatmap::pheatmap(gsva_mat_treg_bulk,
                   scale = "row",
                   annotation_col = ann_col,
                   annotation_colors = list(Celltype = anno_colors),
                   border_color = NA,
                   fontsize = 8,
                   legend = T,
                   annotation_legend = F,
                   show_colnames = F,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()

gsva_obj_treg$group_ann <- paste(gsva_obj_treg$ann_level4_final,
                                 gsva_obj_treg$Section,
                                 gsva_obj_treg$Condition,
                              sep = ":")
dot_rna_cd4 <- DotPlot(gsva_obj_treg, 
                       features = , 
                       assay = "GSVA", 
                       group.by = "group_ann",
                       cols = c("lightgrey", "blue")) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 12, color = "black"),
        axis.title = element_blank(),
        panel.grid = element_blank())
#
tiff(file = glue("{plt_path}dot_rna_cd4.tiff"),
     width = 12, height = 10, units = "in", res = 600, compression = "lzw")
print(dot_rna_cd4)
dev.off()

# ##### Fig. 3E: Monocle2 #####
# sel_term <- c("GOBP-CAMP-MEDIATED-SIGNALING", "GOBP-CGMP-MEDIATED-SIGNALING")
# sample_col <- "Sample_ID_exp"
# block_col <- "Sample_exp"
# ann_col <-  "ann_level4_final"
# group_col <- "Condition"
# section_col <- "Section"
# 
# ## 0. format proportion table
# gsva_meta <- FetchData(gsva_obj_treg, 
#                        vars = c(sel_term,
#                                 group_col, section_col, ann_col, sample_col, block_col))
# gsva_meta$group_comb  <- paste0(gsva_meta[[group_col]], ":", gsva_meta[[section_col]])
# sample_list <- split(gsva_meta$Sample_exp, f = gsva_meta$group_comb)
# #
# gsva_meta$cluster <- gsva_meta[[ann_col]]
# gsva_meta$Section <- factor(gsva_meta[[section_col]], 
#                             levels = c("TI", "Colon"),
#                             labels = c("Ileum", "Colon"))
# gsva_meta$group <- factor(gsva_meta[[group_col]], 
#                           levels = c("NN", "NU", "II"))
# gsva_meta$block <- gsva_meta[[block_col]]
# 
# sel_termx <- sel_term[1]
# gsva_meta$Score <- gsva_meta[[sel_termx]]
# 
# ## 1. paired test between Conditions
# all_ct <- unique(gsva_meta$cluster)
# ppair_prop_df1 <- lapply(1:4, function(x){
#   
#   sample_listx <- pair_list[[x]]
#   sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
#   gsva_meta_use <- subset(gsva_meta, block %in% sample_pair_use)
#   gsva_meta_use$group <- droplevels(gsva_meta_use$group)
#   group_level <- levels(gsva_meta_use$group)
#   ppairx <- lapply(all_ct, function(ctx){
#     #
#     gsva_meta_usex <- subset(gsva_meta_use, cluster == ctx)
#     pair_dfx <- split(gsva_meta_usex, f = gsva_meta_usex$block) %>%
#       lapply(., function(x){
#         propx <- x$Score[match(group_level, x$group)]
#       }) %>% Reduce("rbind", .)
#     wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
#     
#   }) %>% unlist()
#   # ppairx <- c(p.adjust(ppairx, method = "BH"))
#   names(ppairx) <- c(all_ct)
#   return(ppairx)
#   
# }) %>% Reduce("cbind", .) %>% as.data.frame()
# colnames(ppair_prop_df1) <- names(pair_list)[1:4]
# 
# ## 2.1 plot for Ileum
# box_pair_prop_list1_TI <- lapply(all_ct, function(ctx){
#   
#   plt_dfx <- subset(gsva_meta, cluster == ctx & Section == "Ileum")
#   plt_dfx$Condition <- plt_dfx$group
#   max_propx <-max(plt_dfx$Score)
#   padj <- ifelse(ppair_prop_df1[ctx,2:1] < 0.001, "***", 
#                  ifelse(ppair_prop_df1[ctx, 2:1] < 0.01, "**",
#                         ifelse(ppair_prop_df1[ctx,2:1] < 0.05, "*",
#                                ifelse(ppair_prop_df1[ctx,2:1] < 0.1, ".",
#                                       "ns"))))%>% as.vector()
#   diff_dfx <- data.frame(y.position = c(max_propx * 1.1, max_propx * 1.2),
#                          x = c(1.5, 2.5),
#                          xmin = c(1, 1),
#                          xmax = c(2, 3),
#                          group2 = c("NN", "NN"),
#                          group1 = c("NU", "II"),
#                          padj_sig = padj)
#   
#   box_pairx <- ggplot(data = plt_dfx, 
#                       aes(x = Condition, y = Score, color = Condition)) + 
#     geom_boxplot(width = 0.5) +
#     geom_jitter(width = 0.2) + 
#     scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
#     scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
#     ggtitle(ctx) +
#     xlab("") + ylab("Proportion (%) in Ileum") + 
#     stat_pvalue_manual(diff_dfx,
#                        label = "padj_sig",
#                        tip.length = 0.01,
#                        label.size = 5,
#                        bracket.size = 1) +
#     theme_bw() + 
#     theme(legend.position = "none",
#           plot.title = element_text(hjust = 0.5),
#           title = element_text(size = 10, face = "bold", hjust = 0.5),
#           axis.text = element_text(size = 12, color = "black"),
#           axis.title = element_text(size = 15, face = "bold"))
#   
# })
# 
# box_pair_prop_list1_TI[[1]] <- box_pair_prop_list1_TI[[1]] + 
#   theme(axis.title = element_text(size = 15, face = "bold"))
# ggsave(glue("{plt_path}box_pair_level4_Ileum_sel2.png"),
#        patchwork::wrap_plots(box_pair_prop_list1_TI, nrow = 1),
#        height = 4, width = 14, units = "in", dpi = 300)

##### Fig. 3E: Monocle2 #####
library(monocle)
cluster_col <- "ann_level4_final"
group_col <- "Condition"
source("code/FUNCTION/monocle2/monocle_fun_modify.R")
##
monocle_obj <- readRDS("03_output/06_Traj/Monocle2/CD4_Treg/cds_monocle2.rds")
mono_plt <- monocle.plot(monocle_obj = monocle_obj,
                         ct_cols = ct_cols,
                         cluster_col = "ann_level4_final")
mono_plt$traj_plt1 <- mono_plt$traj_plt1 + 
  theme(legend.position = "bottom")
mono_plt$traj_plt2 <- mono_plt$traj_plt2 + 
  theme(legend.position = "bottom")
ggsave(glue("{plt_path}/traj_monocle2.png"), 
       mono_plt$traj_plt1 | mono_plt$traj_plt2,
       height = 5, width = 8, units = "in", dpi = 300)
##
monocle_obj$Condition <- factor(monocle_obj$Condition, 
                        levels = c("NN", "NU", "II"))
traj_plt1 <- monocle::plot_cell_trajectory(monocle_obj, 
                                           color_by = "State",
                                           cell_size = 1,
                                           show_branch_points = F) +
  theme_bw()
ggsave(glue("{plt_path}/traj_monocle_state.png"), 
       traj_plt1,
       height = 4, width = 5, units = "in", dpi = 300)
####
traj_sel <- 4
cds_sub <- monocle_obj[, pData(monocle_obj)$State %in% c(1, 2, 3, traj_sel)]
deg_pseudo <- readRDS(glue("03_output/06_Traj/Monocle2/CD4_Treg/deg_pseudo_monocle2_sub{traj_sel}.rds"))
top_gene_pseudo <- dplyr::top_n(deg_pseudo, n = 60, wt = -qval) %>%
  dplyr::top_n(., n = 60, wt = vf_vst_counts_variance)
ann_col <- make_add_annotation_col(cds_sub, 
                                   cols = c(cluster_col, group_col))
colnames(ann_col) <- c("Cell type", "Condition")
ann_col[["Cell type"]] <- factor(ann_col[["Cell type"]],
                                 levels = c("CD4_Treg_naive", "CD4_IKZF2low_Treg", "CD4_Treg"),
                                 labels = c("CD4 Treg naive", "CD4 IKZF2low Treg", "CD4 Treg"))
##
my_ann_colors <- list(
  Condition = c(NN = gp_cols[1], 
                NU = gp_cols[2], 
                II = gp_cols[3]),
  `Cell type` = c(`CD4 Treg naive` = ct_cols[1], 
                  `CD4 IKZF2low Treg` = ct_cols[2], 
                  `CD4 Treg` = ct_cols[3])
)
#

png(glue("{plt_path}/ptime_heat_monocle2_sub{traj_sel}.png"),
    height = 6, width = 6, units = "in", res = 300)
plot_pseudotime_heatmap_modify(cds_sub[top_gene_pseudo$gene_short_name,],
                               num_clusters = 2,
                               add_annotation_col = ann_col,
                               annotation_colors = my_ann_colors,
                               show_rownames = T,
                               return_heatmap = F)
dev.off()

##### Fig. 3: CCC #####
ccc_path_ti <- "03_output/05_Interaction/multinichenetr/TI/CD4T/DESeq2/ann_level4_final/"
ccc_path_cr <- "03_output/05_Interaction/multinichenetr/Colon/CD4T/DESeq2/ann_level4_final/"
system(glue("cp -f {ccc_path_ti}circos_Treg_to_CD4T_II_vs_NN_top50.pdf {plt_path}/circos_Treg_to_CD4T_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_Treg_to_CD4T_NU_vs_NN_top50.pdf {plt_path}/circos_Treg_to_CD4T_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_Treg_to_CD4T_II_vs_NN_top50.pdf {plt_path}/circos_Treg_to_CD4T_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_Treg_to_CD4T_NU_vs_NN_top50.pdf {plt_path}/circos_Treg_to_CD4T_NU_vs_NN_Colon_top50.pdf"))

system(glue("cp -f {ccc_path_ti}circos_CD4T_to_Treg_II_vs_NN_top50.pdf {plt_path}/circos_CD4T_to_Treg_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_CD4T_to_Treg_NU_vs_NN_top50.pdf {plt_path}/circos_CD4T_to_Treg_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4T_to_Treg_II_vs_NN_top50.pdf {plt_path}/circos_CD4T_to_Treg_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4T_to_Treg_NU_vs_NN_top50.pdf {plt_path}/circos_CD4T_to_Treg_NU_vs_NN_Colon_top50.pdf"))

# ##### Inflammation_score #####
# library(GSVA)
# library(GSEABase)
# 
# deg_level4_ii_nn_ti <- readRDS("03_output/04_Diff/DEG/small_intestine/CD4T/DESeq2/pseudo_bulk_ann_level4_final_I_I_vs_N_N.rds")
# deg_level4_ii_nn_cr <- readRDS("03_output/04_Diff/DEG/colorect/CD4T/DESeq2/pseudo_bulk_ann_level4_final_I_I_vs_N_N.rds")
# deg_cd4_ti <- deg_level4_ii_nn_ti[["CD4_IKZF2low_Treg"]]
# deg_cd4_cr <- deg_level4_ii_nn_cr[["CD4_Treg"]]
# 
# ##
# # top_deg <- deg_cd4_ti$Term[which(deg_cd4_ti$padj < 1E-5 & 
# # abs(deg_cd4_ti$log2FC) > 1)]
# # ref_g <- deg_cd4_ti$Term[which(deg_cd4_ti$padj > 0.5 & 
# # abs(deg_cd4_ti$log2FC) < 0.1)]
# inter_deg_ti <- deg_cd4_ti$Term[which(deg_cd4_ti$padj < 0.05 &
#                                         (deg_cd4_ti$log2FC) > 0.25)]
# inter_deg_cr <- deg_cd4_cr$Term[which(deg_cd4_cr$padj < 0.05 &
#                                           (deg_cd4_cr$log2FC) > 0.25)]
# inter_deg <- c(inter_deg_up, inter_deg_down)
# 
# ##
# bulk_meta <- sc_obj@meta.data[!duplicated(sc_obj$sample_id_exp),]
# bulk_meta <- sc_obj_sub@meta.data[!duplicated(sc_obj_sub$sample_id_exp),]
# sample_col <- "sample_id_exp"
# bulk_meta$group_comb  <- paste0(bulk_meta$condition, ":", bulk_meta$section_comb)
# sample_list <- split(bulk_meta$sample_exp, f = bulk_meta$group_comb)
# rownames(bulk_meta) <- gsub("_", "-", bulk_meta$sample_id_exp)
# clinic_df <- readRDS("03_output/02_clean/clinic_df.rds")
# ssgsea_df  <- readRDS("ssgsea_df.rds")
# bulk_meta$nancy_index_score <- clinic_df$nancy_index_score[match(bulk_meta$sample_id_exp,
#                                                                  clinic_df$sample_id_exp)]
# bulk_meta$rhi_score_calculated <- clinic_df$rhi_score_calculated[match(bulk_meta$sample_id_exp,
#                                                                        clinic_df$sample_id_exp)]
# exp_bulk1 <- AggregateExpression(subset(sc_obj_treg, 
#                                         ann_level4_final == "CD4_IKZF2low_Treg"), 
#                                  assays = "RNA", 
#                                  return.seurat = F, 
#                                  group.by = "sample_id_exp")[["RNA"]]
# exp_bulk_lcpm1 <- edgeR::cpm(exp_bulk1, 
#                              normalized.lib.sizes = T,
#                              log = T)
# exp_bulk2 <- AggregateExpression(subset(sc_obj_treg, 
#                                         ann_level4_final == "CD4_Treg"), 
#                                  assays = "RNA", 
#                                  return.seurat = F, 
#                                  group.by = "sample_id_exp")[["RNA"]]
# exp_bulk_lcpm2 <- edgeR::cpm(exp_bulk2, 
#                              normalized.lib.sizes = T,
#                              log = T)
# saveRDS(exp_bulk1, file = "03_output/03_clustering_test/CD4T/test153000/test302015/exp_bulk1.rds")
# saveRDS(exp_bulk2, file = "03_output/03_clustering_test/CD4T/test153000/test302015/exp_bulk2.rds")
# ## GSVA 
# inter_ti_set <- GeneSetCollection(GeneSet(inter_deg_ti, setName = "inter_deg_ti"))
# inter_cr_set <- GeneSetCollection(GeneSet(inter_deg_cr, setName = "inter_deg_cr"))
# ssgsva_param1 <- ssgseaParam(exprData = exp_bulk_lcpm1,
#                              geneSets = inter_ti_set,
#                              normalize = F)
# gsva_df1 <- gsva(ssgsva_param1)
# ssgsva_param2 <- ssgseaParam(exprData = exp_bulk_lcpm2,
#                              geneSets = inter_cr_set,
#                              normalize = F)
# gsva_df2 <- gsva(ssgsva_param2)
# #
# gsva_df <- data.frame(inter_deg_ti = as.data.frame(t(gsva_df1))[rownames(bulk_meta),], 
#                       inter_deg_cr =as.data.frame(t(gsva_df2))[rownames(bulk_meta),],
#                       row.names = rownames(bulk_meta))
# bulk_meta$score1 <- gsva_df[rownames(bulk_meta), "inter_deg_ti"] %>% scales::rescale(., c(0, 10))
# bulk_meta$score2 <- gsva_df[rownames(bulk_meta), "inter_deg_cr"] %>% scales::rescale(., c(0, 10))
# 
# 
# ## plot
# pair_list <- list("II_vs_NN_TI" = c("I_I:TI", "N_N:TI"),
#                   "NU_vs_NN_TI" = c("N_U:TI", "N_N:TI"),
#                   "II_vs_NN_Colon" = c("I_I:Colon", "N_N:Colon"),
#                   "NU_vs_NN_Colon" = c("N_U:Colon", "N_N:Colon"),
#                   "TI_vs_Colon_NN" = c("N_N:TI", "N_N:Colon"))
# bulk_meta$Section <- factor(bulk_meta$section_comb, levels = c("TI", "Colon"),
#                             labels = c("Ileum", "Colon"))
# bulk_meta$group <- factor(bulk_meta$condition, 
#                           levels = c("N_N", "N_U", "I_I"), 
#                           labels = c("NN", "NU", "II"))
# ##
# bulk_meta1 <- subset(bulk_meta, Section == "Ileum")
# bulk_meta1$score0 <- bulk_meta1$score1
# ppair_df <- lapply(1:2, function(x){
#   
#   sample_listx <- pair_list[[x]]
#   sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
#   bulk_metax <- subset(bulk_meta1, sample_exp %in% sample_pair_use)
#   bulk_metax$group <- droplevels(bulk_metax$group)
#   group_level <- levels(bulk_metax$group)
#   
#   pair_dfx <- split(bulk_metax, f = bulk_metax$sample_exp) %>%
#     lapply(., function(x){
#       propx <- x$score0[match(group_level, x$group)]
#     }) %>% Reduce("rbind", .)
#   ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
#   
#   return(ppairx)
#   
# }) %>% unlist()
# 
# names(ppair_df) <- names(pair_list)[1:2]
# ppair_df_adj <- p.adjust(ppair_df, method = "BH")
# ##
# padj <- ifelse(ppair_df_adj < 0.001, "***", 
#                ifelse(ppair_df_adj < 0.01, "**",
#                       ifelse(ppair_df_adj < 0.05, "*",
#                              ifelse(ppair_df < 0.05, ".",
#                                     "ns"))))%>% as.vector()
# diff_df <- data.frame(y.position = c(10.5, 10),
#                       x = c(1.5, 2.5),
#                       xmin = c(2, 2),
#                       xmax = c(1, 3),
#                       group2 = c("NN", "NN"),
#                       group1 = c("II", "NU"),
#                       Section = c("Ileum", "Ileum"),
#                       padj_sig = padj)
# if_box <- ggplot(bulk_meta1, 
#                  aes(x = group, y = score0, color = group)) +
#   geom_boxplot(fill = NA, width = 0.5) +
#   geom_point(aes(color = group), position = position_jitter(width = 0.1)) +
#   geom_line(aes(group = sample_exp), color = "gray", alpha = 0.5) + 
#   scale_color_manual(values = c("#0073C2FF","#EFC000FF", "#E41A1C")) +
#   ylim(c(0, 11)) + 
#   # scale_y_continuous(limits = c(0, 11))
#   xlim(c("II", "NN", "NU")) + 
#   stat_pvalue_manual(diff_df,
#                      label = "padj_sig",
#                      tip.length = 0.01,
#                      label.size = 5,
#                      bracket.size = 1) +
#   # facet_wrap(~ Section, scales = 'fixed', nrow = 1) +
#   xlab("") + ylab("Inflammation score (Up-regulated in ileal CD4+IKZF2low Tregs)") + 
#   theme_bw() + 
#   theme(legend.position = "none",
#         title = element_text(size = 13, face = "bold"),
#         strip.text = element_text(size = 13, face = "bold"),
#         strip.background = element_blank(),
#         axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 15, face = "bold"))
# ggsave("03_output/03_clustering_test/CD4T/test153000/test302015/if_box_treg1.png",
#        if_box,
#        height = 7, width = 4, units = "in", dpi = 300)
# ##
# bulk_meta2 <- subset(bulk_meta, Section == "Colon")
# bulk_meta2$score0 <- bulk_meta2$score2
# ##
# ppair_df <- lapply(3:4, function(x){
#   
#   sample_listx <- pair_list[[x]]
#   sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
#   bulk_metax <- subset(bulk_meta2, sample_exp %in% sample_pair_use)
#   bulk_metax$group <- droplevels(bulk_metax$group)
#   group_level <- levels(bulk_metax$group)
#   
#   pair_dfx <- split(bulk_metax, f = bulk_metax$sample_exp) %>%
#     lapply(., function(x){
#       propx <- x$score0[match(group_level, x$group)]
#     }) %>% Reduce("rbind", .)
#   ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
#   
#   return(ppairx)
#   
# }) %>% unlist()
# 
# names(ppair_df) <- names(pair_list)[3:4]
# ppair_df_adj <- p.adjust(ppair_df, method = "BH")
# ##
# padj <- ifelse(ppair_df_adj < 0.001, "***", 
#                ifelse(ppair_df_adj < 0.01, "**",
#                       ifelse(ppair_df_adj < 0.05, "*",
#                              ifelse(ppair_df < 0.05, ".",
#                                     "ns"))))%>% as.vector()
# diff_df <- data.frame(y.position = c(10.5, 10),
#                       x = c(1.5, 2.5),
#                       xmin = c(2, 2),
#                       xmax = c(1, 3),
#                       group2 = c("NN", "NN"),
#                       group1 = c("II", "NU"),
#                       Section = c("Colon", "Colon"),
#                       padj_sig = padj)
# if_box <- ggplot(bulk_meta2, aes(x = group, y = score0, color = group)) +
#   geom_boxplot(fill = NA, width = 0.5) +
#   geom_point(aes(color = group), position = position_jitter(width = 0.1)) +
#   geom_line(aes(group = sample_exp), color = "gray", alpha = 0.5) + 
#   scale_color_manual(values = c("#0073C2FF","#EFC000FF", "#E41A1C")) +
#   ylim(c(0, 11)) + 
#   # scale_y_continuous(limits = c(0, 11))
#   xlim(c("II", "NN", "NU")) + 
#   stat_pvalue_manual(diff_df,
#                      label = "padj_sig",
#                      tip.length = 0.01,
#                      label.size = 5,
#                      bracket.size = 1) +
#   # facet_wrap(~ Section, scales = 'fixed', nrow = 1) +
#   xlab("") + ylab("Inflammation score (Up-regulated in colonic CD4+Tregs)") + 
#   theme_bw() + 
#   theme(legend.position = "none",
#         title = element_text(size = 13, face = "bold"),
#         strip.text = element_text(size = 13, face = "bold"),
#         strip.background = element_blank(),
#         axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 15, face = "bold"))
# ggsave("03_output/03_clustering_test/CD4T/test153000/test302015/if_box_treg2.png",
#        if_box,
#        height = 8, width = 4, units = "in", dpi = 300)
# ###
# exp_bulk1 <- AggregateExpression(subset(sc_obj_treg, 
#                                         ann_level4_final == "CD4_IKZF2low_Treg"), 
#                                  assays = "RNA", 
#                                  return.seurat = F, 
#                                  group.by = "sample_id_exp")[["RNA"]]
# exp_bulk_lcpm1 <- edgeR::cpm(exp_bulk1, 
#                              normalized.lib.sizes = T,
#                              log = T) %>% t %>% as.data.frame()
# bulk_meta1 <- subset(bulk_meta, Section == "Ileum")
# bulk_meta1$score0 <- exp_bulk_lcpm1[match(gsub("_", "-", bulk_meta1$sample_id_exp), rownames(exp_bulk_lcpm1)),
#                                     "THADA"]
# ppair_df <- lapply(1:2, function(x){
#   
#   sample_listx <- pair_list[[x]]
#   sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
#   bulk_metax <- subset(bulk_meta1, sample_exp %in% sample_pair_use)
#   bulk_metax$group <- droplevels(bulk_metax$group)
#   group_level <- levels(bulk_metax$group)
#   
#   pair_dfx <- split(bulk_metax, f = bulk_metax$sample_exp) %>%
#     lapply(., function(x){
#       propx <- x$score0[match(group_level, x$group)]
#     }) %>% Reduce("rbind", .)
#   ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
#   
#   return(ppairx)
#   
# }) %>% unlist()
# 
# names(ppair_df) <- names(pair_list)[1:2]
# ppair_df_adj <- p.adjust(ppair_df, method = "BH")
# ##
# padj <- ifelse(ppair_df_adj < 0.001, "***", 
#                ifelse(ppair_df_adj < 0.01, "**",
#                       ifelse(ppair_df_adj < 0.05, "*",
#                              ifelse(ppair_df < 0.05, ".",
#                                     "ns"))))%>% as.vector()
# diff_df <- data.frame(y.position = c(13, 12.5),
#                       x = c(1.5, 2.5),
#                       xmin = c(2, 2),
#                       xmax = c(1, 3),
#                       group2 = c("NN", "NN"),
#                       group1 = c("II", "NU"),
#                       padj_sig = c("*", "ns"))
# if_box <- ggplot(bulk_meta1, 
#                  aes(x = group, y = score0, color = group)) +
#   geom_boxplot(fill = NA, width = 0.5) +
#   geom_point(aes(color = group), position = position_jitter(width = 0.1)) +
#   geom_line(aes(group = sample_exp), color = "gray", alpha = 0.5) + 
#   scale_color_manual(values = c("#0073C2FF","#EFC000FF", "#E41A1C")) +
#   # ylim(c(0, 11)) + 
#   # scale_y_continuous(limits = c(0, 11))
#   xlim(c("II", "NN", "NU")) + 
#   stat_pvalue_manual(diff_df,
#                      label = "padj_sig",
#                      tip.length = 0.01,
#                      label.size = 5,
#                      bracket.size = 1) +
#   xlab("") + ylab("Expression of RORA") + 
#   theme_bw() + 
#   theme(legend.position = "none",
#         title = element_text(size = 13, face = "bold"),
#         strip.text = element_text(size = 13, face = "bold"),
#         strip.background = element_blank(),
#         axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 15, face = "bold"))
# ggsave("03_output/03_clustering_test/CD4T/test153000/test302015/gene_box_treg_THADA.png",
#        if_box,
#        height = 7, width = 4, units = "in", dpi = 300)
# 
# plots <- VlnPlot(sc_obj_treg, 
#                  features = c("THADA", "CSNK1A1", "MALT1", "SATB1-AS1"), 
#                  split.by = "section_group", 
#                  group.by = "ann_level4_final",
#                  pt.size = 0, combine = FALSE)
# 
# ggsave("03_output/03_clustering_test/CD4T/test153000/test302015/gene_vln_treg_test.png",
#        wrap_plots(plots = plots, ncol = 1),
#        height = 15, width = 10, units = "in", dpi = 300)

##### Supplementary Fig. 3: volcano plots #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD4T/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD4T/{use_de}/")
ann_col <- "ann_level4_final"
ct_treg <- levels(sc_obj_treg@meta.data[[ann_col]])
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/volca_plot.R")
## 1. Ileum
deg_df_it_list <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))
sel_gene <- c("RORA", "CCR6", "THADA", "SATB1", "SATB1-AS1", "PDE3B", "PDE4A", "PDE4D", "PDE7B")
lapply(ct_treg, function(ctx){
  
  deg_df_itx <- deg_df_it_list[[ctx]]
  if(nrow(deg_df_itx) > 1){
    deg_df_itx <- subset(deg_df_itx, !is.na(deg_df_itx$padj))
    vc_pltx <- volca_function(deg_df = deg_df_itx,
                              term_col = "Term",
                              fc_col = "log2FC",
                              p_col = "padj",
                              target_gene = sel_gene,
                              n_top = 30,
                              log_thresh = 0.25, 
                              top_log_thresh = 0.25,
                              group_label = c("Down regulation", "Others", "Up regulation"),
                              valco_col = c4a("classic_blue_red12", 3),
                              highlight_col = "orange") + 
      xlab(bquote(~Log[2]~"(Fold Change) in Ileal II vs NN")) +
      ggtitle(gsub("_", " ", ctx)) + 
      theme(plot.title = element_text(size = 15, face = "bold"),
            legend.position = "right")
    
    ggsave(file = glue("{plt_path}/valcano_{ctx}_II_vs_NN_TI.png"),
           vc_pltx,
           width = 9, height = 5,units = "in", dpi = 300, limitsize = T)
    
  }
  return(ctx)
  
})
## 2. Colon
deg_df_cr_list <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))
lapply(ct_treg, function(ctx){
  
  deg_df_crx <- deg_df_cr_list[[ctx]]
  if(nrow(deg_df_crx) > 1){
    deg_df_crx <- subset(deg_df_crx, !is.na(deg_df_crx$padj))
    vc_pltx <- volca_function(deg_df = deg_df_crx,
                              term_col = "Term",
                              fc_col = "log2FC",
                              p_col = "padj",
                              # target_gene = NULL,
                              n_top = 30,
                              log_thresh = 0.25, 
                              top_log_thresh = 0.25,
                              group_label = c("Down regulation", "Others", "Up regulation"),
                              valco_col = c4a("classic_blue_red12", 3),
                              highlight_col = "orange") + 
      xlab(bquote(~Log[2]~"(Fold Change) in Colonic II vs NN")) +
      ggtitle(gsub("_", " ", ctx)) + 
      theme(plot.title = element_text(size = 15, face = "bold"),
            legend.position = "right")
    ggsave(file = glue("{plt_path}/valcano_{ctx}_II_vs_NN_CR.png"),
           vc_pltx,
           width = 9, height = 5,units = "in", dpi = 300, limitsize = T)
    
  }
  return(ctx)
  
})

