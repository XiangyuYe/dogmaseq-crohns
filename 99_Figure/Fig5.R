## model load r/4.5.0
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
library(glue)
library(stringr)
library(reshape2)
library(cols4all)

###data input and select parameters
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/"
setwd(project_path)
#
source("code/FUNCTION/process/PROCESS_FUN.R")
cd8_cols <- c("#03ADF0", "#437DBF", "#7BAFDE", 
              "#E7E419", "#E0C16E", "#FF4500", "#800000", 
              "#EDCC19", "#FF9500", "#F66C2D", "#F0C0C0")
data_path <- "03_output/03_clustering/CD8T/WNN_ADT_RNA/"
plt_path <- "04_plot/Fig5/"
#
pair_list <- list("II_vs_NN_TI" = c("II:Ileum", "NN:Ileum"),
                  "NU_vs_NN_TI" = c("NU:Ileum", "NN:Ileum"),
                  "II_vs_NN_Colon" = c("II:Colon", "NN:Colon"),
                  "NU_vs_NN_Colon" = c("NU:Colon", "NN:Colon"),
                  "TI_vs_Colon_NN" = c("NN:Ileum", "NN:Colon"))
seed_use <- 20250528

##
sc_obj <- readRDS(glue("{data_path}scWNN_obj_drcl.rds"))
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
sc_obj[["RNA"]] <- subset(scRNA_obj, cells = colnames(sc_obj))[["RNA"]]
scATAC_obj <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_addmotif.rds")
sc_obj[["peaks"]] <- subset(scATAC_obj, cells = colnames(sc_obj))[["peaks"]]

#
sc_meta_cd8 <- readRDS(glue("{data_path}sc_meta.rds"))
sc_meta_cd8$Condition <- factor(sc_meta_cd8$Condition, levels = c("NN", "NU", "II"))
sc_meta_cd8$Section <- factor(sc_meta_cd8$Section, levels = c("Colon", "TI"), 
                         labels = c("Colon", "Ileum"))
sc_obj@meta.data <- sc_meta_cd8[colnames(sc_obj),]
Idents(sc_obj) <-  factor(sc_obj$ann_level4_final,
                          levels = levels(sc_obj$ann_level4_final),
                          labels = gsub("_", " ", levels(sc_obj$ann_level4_final)))
##### Fig. 5A: volcano plot #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/all/{use_de}/")
ann_col <- "ann_level2_refine"
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/volca_plot.R")
## 1. Ileum
deg_df_it_cd8 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD8T"]]
deg_df_it_cd8 <- subset(deg_df_it_cd8, !is.na(deg_df_it_cd8$padj))
vc_plt_cd8 <- volca_function(deg_df = deg_df_it_cd8,
                          term_col = "Term",
                          fc_col = "log2FC",
                          p_col = "padj",
                          # target_gene = NULL,
                          n_top = 30,
                          log_thresh = 0.25, 
                          top_log_thresh = 0.5,
                          group_label = c("Down regulation", "Others", "Up regulation"),
                          valco_col = c4a("classic_blue_red12", 3),
                          highlight_col = "orange") + 
  xlab(bquote(~Log[2]~"(Fold Change) in Ileal II vs NN"))

ggsave(file = glue("{plt_path}/valcano_CD8_II_vs_NN_TI.png"),
       vc_plt_cd8,
       width = 10, height = 6,units = "in", dpi = 300, limitsize = T)

##### Fig. 5B: inflammation score #####
library(GSVA)
library(GSEABase)

## DEG list Ileal II vs NN
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/all/{use_de}/")
ann_col <- "ann_level2_refine"
deg_df_it_cd8 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD8T"]]
deg_up <- deg_df_it_cd8$Term[which(deg_df_it_cd8$padj < 0.05 &
                                     deg_df_it_cd8$log2FC > 0.25)]
deg_down <- deg_df_it_cd8$Term[which(deg_df_it_cd8$padj < 0.05 &
                                       deg_df_it_cd8$log2FC < -0.25)]
## pseudobulk expression data
sample_col <- "Sample_ID_exp"
bulk_meta <- sc_meta_cd8[!duplicated(sc_meta_cd8[[sample_col]]),]
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
bulk_meta$score_CD8_up <- gsva_df_up[rownames(bulk_meta), "deg_up"] %>% 
  scales::rescale(., c(0, 10))
bulk_meta$score_CD8_down <- gsva_df_down[rownames(bulk_meta), "deg_down"] %>% 
  scales::rescale(., c(0, 10))
# add NI_score
ssgsea_df  <- readRDS("ssgsva_NI_score.rds")
bulk_meta$if_score_scale <- ssgsea_df[rownames(bulk_meta), "NIScore"] %>% 
  scales::rescale(., c(0, 10))

## plot
bulk_meta$Section <- factor(bulk_meta$Section, 
                            levels = c("Ileum", "Colon"))
bulk_meta$group <- factor(bulk_meta$Condition, 
                          levels = c("NN", "NU", "II"))
## score_CD8_up
bulk_meta$score_CD8 <- bulk_meta$score_CD8_up
ppair_df <- lapply(1:4, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  bulk_metax <- subset(bulk_meta, Sample_exp %in% sample_pair_use)
  bulk_metax$group <- droplevels(bulk_metax$group)
  group_level <- levels(bulk_metax$group)
  
  pair_dfx <- split(bulk_metax, f = bulk_metax$Sample_exp) %>%
    lapply(., function(x){
      propx <- x$score_CD8[match(group_level, x$group)]
    }) %>% Reduce("rbind", .)
  ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
  
  return(ppairx)
  
}) %>% unlist()

padjpair_df <- p.adjust(ppair_df, method = "BH")
names(ppair_df) <- names(padjpair_df) <- names(pair_list)[1:4]

##
padj <- format(padjpair_df, scientific = T, digits = 3) %>%
  gsub("e", "E", .) %>% as.vector()
diff_df <- data.frame(y.position = c(10, 9, 10, 9),
                      x = c(1.5, 2.5, 1.5, 2.5),
                      xmin = c(2, 2, 2, 2),
                      xmax = c(1, 3, 1, 3),
                      group2 = c("NN", "NN", "NN", "NN"),
                      group1 = c("II", "NU", "II", "NU"),
                      Section = c("Ileum", "Ileum", "Colon", "Colon"),
                      padj_sig = padj)
if_box_up <- ggplot(bulk_meta, aes(x = group, y = score_CD8, color = group)) +
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
  facet_wrap(~ Section, scales = 'fixed', nrow = 1) +
  xlab("") + ylab("Inflammation score") + 
  ggtitle("Up-regulated DEGs") + 
  theme_bw() + 
  theme(legend.position = "none",
        title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold"))
ggsave(glue("{plt_path}/if_box_CD8_up.png"),
       if_box_up,
       height = 7, width = 7, units = "in", dpi = 300)

## score_CD8_down
bulk_meta$score_CD8 <- bulk_meta$score_CD8_down
ppair_df <- lapply(1:4, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  bulk_metax <- subset(bulk_meta, Sample_exp %in% sample_pair_use)
  bulk_metax$group <- droplevels(bulk_metax$group)
  group_level <- levels(bulk_metax$group)
  
  pair_dfx <- split(bulk_metax, f = bulk_metax$Sample_exp) %>%
    lapply(., function(x){
      propx <- x$score_CD8[match(group_level, x$group)]
    }) %>% Reduce("rbind", .)
  ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
  
  return(ppairx)
  
}) %>% unlist()

padjpair_df <- p.adjust(ppair_df, method = "BH")
names(ppair_df) <- names(padjpair_df) <- names(pair_list)[1:4]

##
padj <- format(padjpair_df, scientific = T, digits = 3) %>%
  gsub("e", "E", .) %>% as.vector()
diff_df <- data.frame(y.position = c(10, 9, 10, 9),
                      x = c(1.5, 2.5, 1.5, 2.5),
                      xmin = c(2, 2, 2, 2),
                      xmax = c(1, 3, 1, 3),
                      group2 = c("NN", "NN", "NN", "NN"),
                      group1 = c("II", "NU", "II", "NU"),
                      Section = c("Ileum", "Ileum", "Colon", "Colon"),
                      padj_sig = padj)
if_box_down <- ggplot(bulk_meta, aes(x = group, y = score_CD8, color = group)) +
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
  facet_wrap(~ Section, scales = 'fixed', nrow = 1) +
  xlab("") + ylab("Inflammation score") + 
  ggtitle("Down-regulated DEGs") + 
  theme_bw() + 
  theme(legend.position = "none",
        title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold"))
ggsave(glue("{plt_path}/if_box_CD8_down.png"),
       if_box_down,
       height = 7, width = 7, units = "in", dpi = 300)

##### Fig. 5C: Enrichment plot #####
source("code/FUNCTION/diff/Enrichment.R")
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/all/{use_de}/")
ann_col <- "ann_level2_refine"
deg_df_it_cd8 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD8T"]]
deg_up <- deg_df_it_cd8$Term[which(deg_df_it_cd8$padj < 0.05 &
                                     deg_df_it_cd8$log2FC > 0.25)]
deg_down <- deg_df_it_cd8$Term[which(deg_df_it_cd8$padj < 0.05 &
                                       deg_df_it_cd8$log2FC < -0.25)]
#
Enrichment.pipline(gene_list = deg_up,
                   list_name = "Up_IIvsNN_CD8T",
                   num_show = 10,
                   plot = T,
                   GO = T,
                   KEGG = T,
                   Reactome = T,
                   outpath = glue("{deg_path_it}Enrichment"))
Enrichment.pipline(gene_list = deg_down,
                   list_name = "Down_IIvsNN_CD8T",
                   num_show = 10,
                   plot = T,
                   GO = T,
                   KEGG = T,
                   Reactome = T,
                   outpath = glue("{deg_path_it}Enrichment"))
##
enrich_up <- enrich_barplt(prefix = "Up_IIvsNN_CD8T",
                           df_path = glue("{deg_path_it}Enrichment/"),
                           nn = 5,
                           min_gene = 5,
                           max_nchar = 50)
n_up <- nrow(enrich_up$enrich_df)
enrich_bar_up <- enrich_up[["bar_plt"]] + 
  ggtitle("Up-regulated") + 
  theme(plot.title = element_text(size = 15, face = "bold"))
enrich_down <- enrich_barplt(prefix = "Down_IIvsNN_CD8T",
                                 df_path = glue("{deg_path_it}Enrichment/"),
                                 nn = 5,
                                 min_gene = 5,
                                 max_nchar = 50)
n_down <- nrow(enrich_down$enrich_df)
enrich_bar_down <- enrich_down[["bar_plt"]] + 
  ggtitle("Down-regulated") + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ggsave(glue("{plt_path}Enrichment_bar.png"),
       wrap_plots(list(enrich_bar_up, enrich_bar_down), 
                  heights = c(n_up * 0.1 + 5, n_down * 0.1 + 5)),
       height = 8, width = 7, units = "in",  dpi = 300)

##### Fig. 5D: UMAP #####
umap_plt <- DimPlot(sc_obj, 
                    cols = cd8_cols,
                    label = F, 
                    # pt.size = 1E-3,
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
ggsave(file = paste0("03_output/03_clustering_test/CD8T/test153000/test302015/UmapPlot_all.png"),
       umap_plt,
       width = 10, height = 7,units = "in", dpi = 300, limitsize = F)

##### Fig. 5E: markers plot #####
adt_cd8 <- c("TCR-Vd2", "CD142", "CD45RA","CD45RO", 
             "CD62L", "CD69", "CD103", "CD39")
rna_cd8 <- c("CCR7", "SELL", "GPR183",
            "CCL5", "LYST", "GZMB", "GZMK", "KLRB1", "IL26", "IL23R", "CCR6", 
            "HOMER1", "IKZF1", "LY6E", "EGR1",
             "HSPA1A", "HSPA1B", "IL2", "EGR2", "IKZF2")
atac_cd8 <- c("TBX21", "IL17A", "RORC", "IKZF2")

sel_motif <- c("MA0690.1", "MA1151.1", "MA0071.1", "MA0072.1", 
               "MA0740.1", "MA0079.4", 
               "MA1141.1", "MA0476.1")
names_motif <- c("TBX21", "RORC", "RORA.1", "RORA.2", 
                 "KLF14", "SP1", 
                 "FOS::JUND", "FOS")

Idents(sc_obj) <- sc_obj$ann_level4_final
dot_adt_cd8 <- vln.plt(seurat_obj = sc_obj,
                       clus_col = "ann_level4_final",
                       assay_use = "ADT",
                       feature_list = adt_cd8,
                       color_use = cd8_cols)
#
dot_rna_cd8 <- dot.minmax(seurat_obj = sc_obj,
                          assay = "RNA",
                          features = rna_cd8,
                          col_min = 0,
                          col_max = 1,
                          group = "ann_level4_final") + 
  scale_fill_distiller(type = 'div', palette = 'RdYlBu')+
  scale_y_discrete(limits = levels(sc_obj$ann_level4_final) %>% rev) + 
  theme(legend.position = "none",
        axis.text.y =  element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
#
dot_mt_cd8 <- dot.minmax(seurat_obj = sc_obj,
                          assay = "chromvar",
                          features = sel_motif,
                          col_min = 0,
                          col_max = 1,
                          group = "ann_level4_final") + 
  scale_fill_distiller(type = 'div', palette = 'RdYlBu')+
  scale_x_discrete(limits = sel_motif, labels = names_motif) + 
  scale_y_discrete(limits = levels(sc_obj$ann_level4_final) %>% rev) + 
  theme(legend.position = "none",
        axis.text.y =  element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
dot_plt <- (dot_adt_cd8 | dot_rna_cd8 | dot_mt_cd8) + 
  plot_layout(widths = c(length(adt_cd8) + 2, 
                         length(rna_cd8),
                         length(sel_motif)))

ggsave(file = glue("{plt_path}dot_plt_CD8.png"),
       dot_plt,
       height = 6, width = 15,units = "in", dpi = 300, limitsize = F)
##
cv_plt <- peak.set.plt(sc_obj = sc_obj,
                        assay_use = "peaks",
                        mk_list = atac_cd8,
                        extend_kb_up = 3000,
                        extend_kb_down = 0,
                        color_use = cd8_cols)
cv_plt[[1]] <- cv_plt[[1]] + theme(strip.text.y.left = element_blank())
dot_plt2 <- (dot_adt_cd8 | dot_rna_cd8 | cv_plt) + 
  plot_layout(widths = c(length(adt_cd8) + 2, 
                         length(rna_cd8),
                         length(atac_cd8) * 2))
#
ggsave(file = glue("{plt_path}dot_plt_CD8_cv.png"),
       dot_plt2,
       height = 6, width = 15,units = "in", dpi = 300, limitsize = F)

##### Fig. 5F: Heatmap for ct-specific genes #####
DefaultAssay(sc_obj) <- "RNA"
Idents(sc_obj) <- sc_obj$ann_level4_final
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
                                         group.by = "ann_level4_final",
                                         assays = "RNA",
                                         return.seurat = T,
                                         layer = "data")
ann_level <- gsub("\\-", " ", colnames(sc_obj_bulk))
ann_col <- data.frame(Celltype = factor(ann_level,
                                        levels = ann_level),
                      row.names = colnames(sc_obj_bulk))
anno_colors <- cd8_cols
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

####### Fig. 5G: ct-specific eRegulon #########
sp_outpath <- "03_output/07_SCENIC/CD8T/scplus_pipeline/Snakemake/"
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))

scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
scAUC_obj_trm <- subset(scAUC_obj, cells = rownames(sc_meta_cd8)[sc_meta_cd8$ann_level3_final == "CD8_TRM"])
scAUC_obj_trm@meta.data <- sc_meta_cd8[colnames(scAUC_obj_trm),]
#
Idents(scAUC_obj_trm) <- scAUC_obj_trm$ann_level4_final
er_gene <- rownames(scAUC_obj_trm)[grep("g\\)", rownames(scAUC_obj_trm))]
diff_auc <- FindAllMarkers(scAUC_obj_trm, 
                           features = er_gene,
                           only.pos = T)
saveRDS(diff_auc, file = glue("{data_path}diff_auc.rds"))
sig_diff <- diff_auc %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  top_n(5, -p_val_adj) %>%
  ungroup()
sig_diff <- sig_diff[order(sig_diff$p_val_adj, -abs(sig_diff$avg_log2FC)),]
sig_diff <- sig_diff[!duplicated(sig_diff$gene),]
sig_diff <- sig_diff[order(sig_diff$cluster, sig_diff$p_val_adj, -abs(sig_diff$avg_log2FC)),]
sig_marker <- str_split_i(sig_diff$gene, "-\\(", 1)

##
heatdot_CD8trm <- fread2(glue("{sp_outpath}heatdot_trm_ann_level4_final.txt"))[,-1]
colnames(heatdot_CD8trm) <- c("celltype", "Gene-based AUC", "eRegulon_name", "Region-based AUC")
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))
e_regulon_name_trans$heatdot_name <- stringr::str_split_i(e_regulon_name_trans$signature, "_\\(", 1)
heatdot_CD8trm$signature_name_simplify <- 
  e_regulon_name_trans$signature_name_simplify[match(heatdot_CD8trm$eRegulon_name,
                                                     e_regulon_name_trans$heatdot_name)] %>%
  str_split_i(., "_", 1)

##
sig_marker <- sig_marker[which(sig_marker %in% heatdot_CD8trm$signature_name_simplify)]
heatdot_CD8trm_sub <- subset(heatdot_CD8trm, signature_name_simplify %in% sig_marker)
heatdot_CD8trm_sub$signature_name_simplify <- factor(heatdot_CD8trm_sub$signature_name_simplify,
                                                     levels = rev(sig_marker))
heatdot_CD8trm_sub$celltype <- factor(heatdot_CD8trm_sub$celltype,
                                      levels = levels(scAUC_obj_trm$ann_level4_final),
                                      labels = gsub("_", " ", levels(scAUC_obj_trm$ann_level4_final)))
heatdot_CD8trm_sub$celltype <- droplevels(heatdot_CD8trm_sub$celltype)
heat_dot <- ggplot(heatdot_CD8trm_sub, 
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
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10, color = "black"))
#
df_annot <- data.frame(celltype = levels(heatdot_CD8trm_sub$celltype))
p_annot = ggplot(df_annot, aes(x = celltype, y = 1, fill = celltype)) +
  scale_fill_manual(values = cd8_cols[-c(1:3)]) + 
  geom_tile() +
  theme_void() +
  theme(legend.position = "none")

ggsave(glue("{plt_path}heat_dot_CD8trm_ann_level4_final.png"),
       wrap_plots(list(p_annot, heat_dot), ncol = 1, heights = c(0.1, 8)),
       height = 15, width = 8, units = "in", dpi = 300)

##### Fig. 6H: pie plot #####
ann_col <-  "ann_level4_final"
group_col <- "Condition"
section_col <- "Section"
##
sc_meta_cd8$group_comb <- paste0(sc_meta_cd8[[group_col]], ":", sc_meta_cd8[[section_col]])
tab_prop <- table(sc_meta_cd8[[ann_col]], sc_meta_cd8$group_comb) %>% 
  prop.table(2) %>% as.data.frame()
colnames(tab_prop) <- c("Cell type", "group_comb", "Proportion")
tab_prop$Proportion <- tab_prop$Proportion * 100
tab_prop$Section <- str_split_i(tab_prop$group_comb, ":", 2)
tab_prop$Group <- str_split_i(tab_prop$group_comb, ":", 1)

tab_prop$Group <- factor(tab_prop$Group, 
                         levels = c("NN", "NU", "II"))
##
pie_plt <- ggplot(tab_prop, aes(x = 3, y = Proportion, fill = `Cell type`)) + 
  geom_col(width = 1.5, color = NA) + 
  facet_grid(Group ~ Section, switch = "y") + 
  coord_polar(theta = "y") + 
  xlim(c(0.2, 3.8)) + 
  scale_fill_manual(values = cd8_cols) + 
  theme_void()+
  theme(strip.text = element_text(size = 15, face = "bold"), 
        legend.title = element_text(size = 15, face = "bold"), 
        legend.text = element_text(size = 12))
ggsave(file = glue("{plt_path}pie_plt_CD8.png"),
       pie_plt,
       width = 6, height = 6,units = "in", dpi = 300)

##### Fig. 5I: miloR #####
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
  sc_obj_sub$ann_level4_final <- factor(sc_obj_sub$ann_level4_final,
                                         levels = levels(sc_obj_sub$ann_level4_final),
                                         labels = gsub("_", " ",
                                                       levels(sc_obj_sub$ann_level4_final)))
  #
  res_milo <- milor.seurat(seurat_obj = sc_obj_sub,
                           n_dims = 30,
                           k_nn = 50,
                           group_col = "Condition",
                           cov_col = "Sample_exp",
                           block_col = NULL,
                           sample_col = "Sample_ID_exp",
                           ann_col = "ann_level4_final",
                           use_reduc = "harmony_SCT",
                           plot_reduc = "wnn.umap",
                           seed_use = seed_use,
                           n_core = 10)
  res_milo$plt_milo <- res_milo$plt_milo + 
    ggtitle(all_pair_set[pair_setx]) +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  ##
  plt_milo_comb <- (res_milo$plt_milo | res_milo$plt_milo2) + 
    plot_layout(widths = c(6, 4))
  ggsave(glue("{plt_path}traj_milo_{pair_setx}.png"),
         plt_milo_comb,
         height = 7, width = 14, units = "in", dpi = 300)
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
sc_obj_sub$ann_level4_final <- factor(sc_obj_sub$ann_level4_final,
                                       levels = levels(sc_obj_sub$ann_level4_final),
                                       labels = gsub("_", " ",
                                                     levels(sc_obj_sub$ann_level4_final)))
#
res_milo <- milor.seurat(seurat_obj = sc_obj_sub,
                         n_dims = 30,
                         k_nn = 50,
                         group_col = "Section",
                         cov_col = "Sample_exp",
                         block_col = NULL,
                         sample_col = "Sample_ID_exp",
                         ann_col = "ann_level4_final",
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

##### Fig. 5J and Extended Fig. 5: Diff prop ######
sample_col <- "Sample_ID_exp"
block_col <- "Sample_exp"
ann_col <-  "ann_level4_final"
group_col <- "Condition"
section_col <- "Section"

## 0. format proportion table
bulk_meta <- sc_meta_cd8[!duplicated(sc_meta_cd8[[sample_col]]),]
bulk_meta$group_comb  <- paste0(bulk_meta[[group_col]], ":", bulk_meta[[section_col]])
sample_list <- split(bulk_meta$Sample_exp, f = bulk_meta$group_comb)
#
tab_prop <- table(sc_meta_cd8[[ann_col]], sc_meta_cd8[[sample_col]]) %>% 
  prop.table(2) %>% as.data.frame()
colnames(tab_prop) <- c("cluster", "sample", "Proportion")
tab_prop$Proportion <- tab_prop$Proportion * 100
tab_prop$block <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                            block_col]
tab_prop$group <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                            group_col]
tab_prop$Section <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                              section_col] %>%
  factor(., levels = c("Ileum", "Colon"))
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
sel_ct <- c("CD8_TCM", "CD8_TEM_Tc1", "CD8_TRM_HSP")
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
  
  ylabx <- ifelse(ctx %in% sel_ct[lab_facet], "Proportion (%) in Ileum", "")
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

ggsave(glue("{plt_path}box_pair_level4_Ileum_sel.png"),
       patchwork::wrap_plots(box_pair_prop_list1_TI, nrow = 1),
       height = 5, width = 7, units = "in", dpi = 300)

## 2.2 plot for Colon
lab_facet <- c(1)
box_pair_list1_Colon <- lapply(sel_ct, function(ctx){
  
  plt_dfx <- subset(tab_prop, cluster == ctx & Section == "Colon")
  plt_dfx$Condition <- plt_dfx$group
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df1[ctx,4:3], scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1, max_propx * 1.2),
                         x = c(1.5, 2.5),
                         xmin = c(1, 1),
                         xmax = c(2, 3),
                         group2 = c("NN", "NN"),
                         group1 = c("NU", "II"),
                         padj_sig = padj)
  
  ylabx <- ifelse(ctx %in% all_ct[lab_facet], "Proportion (%) in Colon", "")
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

ggsave(glue("{plt_path}box_pair_level4_Colon_sel.png"),
       patchwork::wrap_plots(box_pair_list1_Colon, nrow = 1),
       height = 5, width = 7, units = "in", dpi = 300)

## 3. paired test between Sections
all_ct <- unique(tab_prop$cluster)
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
lab_facet <- c(1, 7)
box_pair_list2 <- lapply(all_ct, function(ctx){
  
  plt_dfx <- subset(tab_prop, cluster == ctx & group == "NN")
  plt_dfx$Condition <- plt_dfx$Section
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df2[ctx,1], scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1),
                         x = c(1.5),
                         xmin = c(1),
                         xmax = c(2),
                         group2 = c("Ileum"),
                         group1 = c("Colon"),
                         padj_sig = padj)
  
  ylabx <- ifelse(ctx %in% all_ct[lab_facet], "Proportion (%) in NN", "")
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Section, y = Proportion, color = Section)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("#0073C2FF", "#E41A1C")) + 
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

ggsave(glue("{plt_path}box_pair_level4_NN.png"),
       patchwork::wrap_plots(box_pair_list2, nrow = 2),
       height = 7, width = 12, units = "in", dpi = 300)

##### Fig. 5K: scDRS #####
scDRS_out_path <- glue("03_output/01_scDRS/04_downstream/cov/")
##
scDRS_comb_level4 <- fread2(glue("{scDRS_out_path}CD.scdrs_group.ann_level4_final"))
scDRS_comb_level4_CD8 <- subset(scDRS_comb_level4, group %in% levels(sc_meta_cd8$ann_level4_final))
scDRS_comb_level4_CD8$section <- "Combined"
scDRS_comb_level4_CD8$celltype <- factor(scDRS_comb_level4_CD8$group, 
                                         levels = levels(sc_meta_cd8$ann_level4_final) %>% rev,
                                         labels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final) %>% rev))
##
scDRS_level4 <- fread2(glue("{scDRS_out_path}CD.scdrs_group.sample_ct_level4"))
scDRS_level4$section <- str_split_i(scDRS_level4$group, ":", 1)
scDRS_level4$celltype <- str_split_i(scDRS_level4$group, ":", 2)
scDRS_level4_CD8 <- subset(scDRS_level4, celltype %in% levels(sc_meta_cd8$ann_level4_final))
scDRS_level4_CD8$celltype <- factor(scDRS_level4_CD8$celltype, 
                                    levels = levels(sc_meta_cd8$ann_level4_final) %>% rev,
                                    labels = gsub("_", " ", levels(sc_meta_cd8$ann_level4_final) %>% rev))
## 
scDRS_CD8_plt <- rbind(scDRS_comb_level4_CD8, scDRS_level4_CD8)
scDRS_CD8_plt$section <- factor(scDRS_CD8_plt$section, 
                                levels = c("Combined", "TI", "RC"),
                                labels = c("Combined", "Ileum", "Colon"))
scDRS_CD8_plt <- split(scDRS_CD8_plt, f = scDRS_CD8_plt$section) %>%
  lapply(., function(scDRS_CD8_pltx){
    scDRS_CD8_pltx$FDR <- p.adjust(scDRS_CD8_pltx$assoc_mcp, method = "BH")
    scDRS_CD8_pltx$hetero_padj <- p.adjust(scDRS_CD8_pltx$hetero_mcp, method = "BH")
    return(scDRS_CD8_pltx)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
scDRS_CD8_plt$assoc_sig <- scDRS_CD8_plt$FDR < 0.05
scDRS_CD8_plt$hetero_sig <- scDRS_CD8_plt$hetero_padj < 0.05


scdrs_plt <- ggplot() + 
  geom_tile(data = scDRS_CD8_plt, 
            aes(x = section, y = celltype, fill = -log10(FDR))) + 
  scale_fill_viridis_c() + 
  geom_tile(data = subset(scDRS_CD8_plt, assoc_sig), 
            aes(x = section, y = celltype),
            color = "black", fill = NA, size = 1, show.legend = T) +
  geom_text(data = subset(scDRS_CD8_plt, hetero_sig),
            aes(x = section, y = celltype, label = "\u2716"),
            color = "black", size = 6, show.legend = T) + 
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10))
#
ggsave(file = glue("{plt_path}CD_sample_ct_CD8T.png"),
       scdrs_plt,
       width = 6, height = 8,units = "in", dpi = 300, limitsize = F)

##### Fig. 5L: volcano plot for Tc17 #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD8T/{use_de}/")
ann_col <- "ann_level4_final"
sel_ct <- "CD8_TRM_Tc17_IL26"
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/volca_plot.R")
## 1. Ileum
deg_df_it_cd8x <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[sel_ct]]
deg_df_it_cd8x <- subset(deg_df_it_cd8x, !is.na(deg_df_it_cd8x$padj))
vc_plt_cd8x <- volca_function(deg_df = deg_df_it_cd8x,
                              term_col = "Term",
                              fc_col = "log2FC",
                              p_col = "padj",
                              # target_gene = NULL,
                              n_top = 30,
                              log_thresh = 0.25, 
                              top_log_thresh = 0.5,
                              group_label = c("Down regulation", "Others", "Up regulation"),
                              valco_col = c4a("classic_blue_red12", 3),
                              highlight_col = "orange") + 
  ggtitle(gsub("_", " ", sel_ct)) + 
  xlab(bquote(~Log[2]~"(Fold Change) in Ileal II vs NN")) + 
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.position = "right")

ggsave(file = glue("{plt_path}/valcano_CD8_Tc17_II_vs_NN_TI.png"),
       vc_plt_cd8x,
       width = 12, height = 6,units = "in", dpi = 300, limitsize = T)

##### Fig. 5M: Enrichment plot Tc17 #####
source("code/FUNCTION/diff/Enrichment.R")
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD8T/{use_de}/")
ann_col <- "ann_level4_final"
sel_ct <- "CD8_TRM_Tc17_IL26"
deg_df_it_cd8x <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[sel_ct]]
deg_upx <- deg_df_it_cd8x$Term[which(deg_df_it_cd8x$padj < 0.05 &
                                      deg_df_it_cd8x$log2FC > 0.25)]
deg_downx <- deg_df_it_cd8x$Term[which(deg_df_it_cd8x$padj < 0.05 &
                                         deg_df_it_cd8x$log2FC < -0.25)]
#
Enrichment.pipline(gene_list = deg_upx,
                   list_name = glue("Up_IIvsNN_{sel_ct}"),
                   num_show = 10,
                   plot = T,
                   GO = T,
                   KEGG = T,
                   Reactome = T,
                   outpath = glue("{deg_path_it}Enrichment"))
Enrichment.pipline(gene_list = deg_downx,
                   list_name = glue("Down_IIvsNN{sel_ct}"),
                   num_show = 10,
                   plot = T,
                   GO = T,
                   KEGG = T,
                   Reactome = T,
                   outpath = glue("{deg_path_it}Enrichment"))
##
enrich_upx <- enrich_barplt(prefix = glue("Up_IIvsNN_{sel_ct}"),
                           df_path = glue("{deg_path_it}Enrichment/"),
                           nn = 5,
                           min_gene = 5,
                           max_nchar = 50)
n_upx <- nrow(enrich_upx$enrich_df)
enrich_bar_upx <- enrich_upx[["bar_plt"]] + 
  ggtitle("Up-regulated") + 
  theme(plot.title = element_text(size = 15, face = "bold"))
enrich_downx <- enrich_barplt(prefix = glue("Down_IIvsNN{sel_ct}"),
                             df_path = glue("{deg_path_it}Enrichment/"),
                             nn = 5,
                             min_gene = 5,
                             max_nchar = 50)
n_downx <- nrow(enrich_downx$enrich_df)
enrich_bar_downx <- enrich_downx[["bar_plt"]] + 
  ggtitle("Down-regulated") + 
  theme(plot.title = element_text(size = 15, face = "bold"))

ggsave(glue("{plt_path}Enrichment_bar_{sel_ct}.png"),
       wrap_plots(list(enrich_bar_upx, enrich_bar_downx), 
                  heights = c(n_upx * 0.1 + 5, n_downx * 0.1 + 5)),
       height = 10, width = 7, units = "in",  dpi = 300)


##### Extended Fig. 8: cNMF #####
sub_idx <- "CD8_TRM"
n_k <- 9
nmf_outpath <- glue("03_output/03_clustering/{sub_idx}/cNMF/")
nmf_usage <- fread2(file = glue("{nmf_outpath}/{sub_idx}_harmony_usage_k{n_k}.txt"),
                      header = T) %>%
  tibble::column_to_rownames(., "V1")
top_gene <- fread2(glue("{nmf_outpath}/{sub_idx}_harmony_top_genes_k{n_k}.txt"))
##
sc_obj_sub <- readRDS("03_output/03_clustering/CD8T/WNN_ADT_RNA/scWNN_obj_sub1.rds")
sc_obj_sub[["NMF"]] <- CreateAssayObject(t(nmf_usage),   
                                       min.cells = 0,
                                       min.features = 0)
##
DefaultAssay(sc_obj_sub) <- "NMF"
ft_NMF <- FeaturePlot(sc_obj_sub,
                       reduction = "wnn.umap",
                       raster = F,
                       features = colnames(nmf_usage),
                       max.cutoff = "q95", 
                       min.cutoff = "q5",
                       ncol = 2)
ft_NMF <- lapply(ft_NMF, function(ft_NMFx){
  
  ft_NMFx + 
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
Idents(sc_obj_sub) <- sc_obj_sub$ann_level4_final
umap_sub <- DimPlot(sc_obj_sub, 
                    reduction = "wnn.umap",
                    cols = cd8_cols[-c(1:3)],
                    label = T, 
                    label.size = 3,
                    raster = F)+
  xlab("wnnUMAP1") + ylab("wnnUMAP2")+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = "bold", size = 15),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#
ggsave(glue("{plt_path}feature_NMF_{n_k}.png"), 
       wrap_plots(c(umap_sub, ft_NMF), ncol = 2),
       height = 20, width = 8, units = "in", dpi = 300)

##
GEP_order <- paste0("GEP",
                    c(6, 4, 9, 7, 3, 5, 1, 8, 2))
sc_obj_bulk <- Seurat::AverageExpression(sc_obj_sub,
                                         # features = GEP_order,
                                         group.by = "ann_level4_final",
                                         assays = "NMF",
                                         return.seurat = T,
                                         layer = "data")

ann_level <- gsub("\\-", " ", colnames(sc_obj_bulk))
ann_col <- data.frame(Celltype = factor(ann_level,
                                        levels = ann_level),
                      row.names = colnames(sc_obj_bulk))
anno_colors <- cd8_cols[-c(1:3)]
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
    height = 6, width = 8, units = "in", res = 300)
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

##### Extended Fig. 5: Number of DEGs #####
#### 1. DEGs
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD8T/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD8T/{use_de}/")
ann_col <- "ann_level4_final"
source("code/FUNCTION/diff/volca_plot.R")
## 1.1 DEGs bewteen II and NN in Ileal CD8T
deg_df_ti <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
deg_df_ti$cluster <- factor(deg_df_ti$cluster, 
                            levels = levels(sc_meta_cd8[[ann_col]]),
                            labels = gsub("_", " ", levels(sc_meta_cd8[[ann_col]])))
plt_deg_ti <- diff.num.plot(diff_df = deg_df_ti,
                            cluster_col = "cluster",
                            p_adj_col = "padj",
                            fc_col = "log2FC",
                            p_adj_thresh = 0.05,
                            logfc_thresh = 0.25,
                            col_use = c("Down" = "#2C69B0", "Up" = "#F02720")) + 
  xlab("Number of DEGs in Ileal II vs NN")
#
ggsave(file = glue("{plt_path}/ndeg_{ann_col}_II_vs_NN_ti.png"),
       plt_deg_ti,
       height = 6, width = 6.5,units = "in", dpi = 300, limitsize = T)
## 1.2 DEGs bewteen II and NN in Colonic CD8T
deg_df_cr <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
deg_df_cr$cluster <- factor(deg_df_cr$cluster, 
                            levels = levels(sc_meta_cd8[[ann_col]]),
                            labels = gsub("_", " ", levels(sc_meta_cd8[[ann_col]])))
plt_deg_cr <- diff.num.plot(diff_df = deg_df_cr,
                            cluster_col = "cluster",
                            p_adj_col = "padj",
                            fc_col = "log2FC",
                            p_adj_thresh = 0.05,
                            logfc_thresh = 0.25,
                            col_use = c("Down" = "#2C69B0", "Up" = "#F02720")) + 
  xlab("Number of DEGs in Colonic II vs NN")
#
ggsave(file = glue("{plt_path}/ndeg_{ann_col}_II_vs_NN_cr.png"),
       plt_deg_cr,
       height = 6, width = 6.5,units = "in", dpi = 300, limitsize = T)

#### 2. DARs
## 2.1 DARs bewteen II and NN in Ileal CD8T
out_path <- "03_output/04_Diff/DAR/"
dar_path_it <- glue("{out_path}TI/CD8T/{use_de}/")
dar_path_cr <- glue("{out_path}Colon/CD8T/{use_de}/")
ann_col <- "ann_level4_final"
##
dar_df_ti <- readRDS(glue("{dar_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
dar_df_ti$cluster <- factor(dar_df_ti$cluster, 
                            levels = levels(sc_meta_cd8[[ann_col]]),
                            labels = gsub("_", " ", levels(sc_meta_cd8[[ann_col]])))
plt_dar_ti <- diff.num.plot(diff_df = dar_df_ti,
                            cluster_col = "cluster",
                            p_adj_col = "padj",
                            fc_col = "log2FC",
                            p_adj_thresh = 0.05,
                            logfc_thresh = 0.25,
                            col_use = c("Down" = "#2C69B0", "Up" = "#F02720")) + 
  xlab("Number of DARs in Colonic II vs NN")
#
ggsave(file = glue("{plt_path}/ndar_{ann_col}_II_vs_NN_ti.png"),
       plt_dar_ti,
       height = 8, width = 8,units = "in", dpi = 300, limitsize = T)

## 2.2 DARs bewteen II and NN in Colonic CD8T
dar_df_cr <- readRDS(glue("{dar_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
dar_df_cr$cluster <- factor(dar_df_cr$cluster, 
                            levels = levels(sc_meta_cd8[[ann_col]]),
                            labels = gsub("_", " ", levels(sc_meta_cd8[[ann_col]])))
print(nrow(dar_df_cr))

##### Fig. 3: CCC #####
ccc_path_ti <- "03_output/05_Interaction/multinichenetr/TI/CD8T/DESeq2/ann_level4_final/"
ccc_path_cr <- "03_output/05_Interaction/multinichenetr/Colon/CD8T/DESeq2/ann_level4_final/"
system(glue("cp -f {ccc_path_ti}circos_CD8T_to_CD4T_II_vs_NN_top50.pdf {plt_path}/circos_CD8T_to_CD4T_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_CD8T_to_CD4T_NU_vs_NN_top50.pdf {plt_path}/circos_CD8T_to_CD4T_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD8T_to_CD4T_II_vs_NN_top50.pdf {plt_path}/circos_CD8T_to_CD4T_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD8T_to_CD4T_NU_vs_NN_top50.pdf {plt_path}/circos_CD8T_to_CD4T_NU_vs_NN_Colon_top50.pdf"))

system(glue("cp -f {ccc_path_ti}circos_CD8T_to_CD8T_II_vs_NN_top50.pdf {plt_path}/circos_CD8T_to_CD8T_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_CD8T_to_CD8T_NU_vs_NN_top50.pdf {plt_path}/circos_CD8T_to_CD8T_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD8T_to_CD8T_II_vs_NN_top50.pdf {plt_path}/circos_CD8T_to_CD8T_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD8T_to_CD8T_NU_vs_NN_top50.pdf {plt_path}/circos_CD8T_to_CD8T_NU_vs_NN_Colon_top50.pdf"))

