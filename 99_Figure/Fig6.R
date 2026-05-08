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

# set parameters
source("code/FUNCTION/process/PROCESS_FUN.R")
pair_list <- list("II_vs_NN_TI" = c("II:Ileum", "NN:Ileum"),
                  "NU_vs_NN_TI" = c("NU:Ileum", "NN:Ileum"),
                  "II_vs_NN_Colon" = c("II:Colon", "NN:Colon"),
                  "NU_vs_NN_Colon" = c("NU:Colon", "NN:Colon"),
                  "TI_vs_Colon_NN" = c("NN:Ileum", "NN:Colon"))
b_cols <- c("#1965B0", "#84D4FA", 
            "#FF4500", "#F1932D", "#F7F056",
            "#72190E", "#DC050C" ,
            "#4EB265", "#94BA19")
seed_use <- 20250528

# set path
plt_path = "04_plot/Fig6/"
data_path = "03_output/03_clustering/B/WNN_ADT_RNA/"
##### load data #####
sc_obj <- readRDS("03_output/03_clustering/B/WNN_ADT_RNA/scWNN_obj.rds")
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)
#
sc_meta_b <- readRDS("03_output/03_clustering/B/WNN_ADT_RNA/sc_meta.rds")
sc_meta_b$Condition <- factor(sc_meta_b$Condition, levels = c("NN", "NU", "II"))
sc_meta_b$Section <- factor(sc_meta_b$Section, levels = c("Colon", "TI"), 
                            labels = c("Colon", "Ileum"))
sc_obj@meta.data <- sc_meta_b
Idents(sc_obj) <- sc_obj$ann_level4_final <- sc_meta_b$ann_level4_final %>% 
  factor(.,
         levels = levels(sc_meta_b$ann_level4_final),
         labels = gsub("_", " ", levels(sc_meta_b$ann_level4_final)))


##### Fig. 6C: UMAP #####
umap_plt <- DimPlot(sc_obj, 
                    cols = b_cols,
                    label = F, 
                    reduction = "wnn.umap",
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
ggsave(file = glue("{plt_path}UmapPlot_B.png"),
       umap_plt,
       width = 10, height = 7,units = "in", dpi = 300, limitsize = F)

##### Fig. 6D: Markers #####
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/process_ATAC.R")
adt_b <- c("CD19","IgD", "CD4", "CD11c", "CD20", "CD27", "CD86", "HLA-DR")
rna_b <- c("IGHD", "FCER2", "EGR1", "AIM2", "TNFRSF13B",
           "TOX", "FCRL3", "FCRL4", "ITGAX",
           'TCL1A', 'BCL6', "RGS13", 
           "CD38", "MZB1", "IGHA1", "IGHA2", "IGHG1", "IGHG2")
Idents(sc_obj) <- sc_obj$ann_level4_final
dot_adt_b <- vln.plt(seurat_obj = sc_obj,
                       clus_col = "ann_level4_final",
                       assay_use = "ADT",
                       feature_list = adt_b,
                       color_use = b_cols) + 
  ggtitle("ADT markers") + 
  theme(plot.title = element_text(size = 11, face = "bold"))
#
dot_rna_b <- dot.minmax(seurat_obj = sc_obj,
                          assay = "RNA",
                          features = rna_b,
                          col_min = 0,
                          col_max = 1,
                          group = "ann_level4_final") + 
  scale_fill_distiller(type = 'div', palette = 'RdYlBu')+
  scale_y_discrete(limits = levels(sc_obj$ann_level4_final) %>% rev) + 
  ggtitle("RNA markers") + 
  theme(legend.position = "none",
        axis.text.y =  element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 11, face = "bold"))
#
dot_plt <- (dot_adt_b | dot_rna_b) + 
  plot_layout(widths = c(length(adt_b) + 2, 
                         length(rna_b)))
ggsave(file = glue("{plt_path}dot_plt_B.png"),
       dot_plt,
       height = 5, width = 12,units = "in", dpi = 300, limitsize = F)

##### Fig. 6E: Heatmap for ct-specific genes #####
DefaultAssay(sc_obj) <- "RNA"
# DEGs between B and Plasma
Idents(sc_obj) <- sc_obj$ann_level2_refine
deg_level2<- FindAllMarkers(sc_obj, only.pos = T)

# DEGs among B or Plasma subtypes
deg_level4_by2 <- SplitObject(sc_obj, split.by = "ann_level2_refine") %>%
  lapply(., function(sc_objx){
    
    Idents(sc_objx) <- sc_objx$ann_level4_final
    if (length(unique(sc_objx$ann_level4_final)) > 1) {
      deg_level4_by2x <- FindAllMarkers(sc_objx, only.pos = T)
    } else {
      deg_level4_by2x <- NULL
    }
    return(deg_level4_by2x)
  })
#
ct_markers <- rbind(deg_level2, Reduce("rbind", deg_level4_by2) %>% as.data.frame())
saveRDS(ct_markers, file = glue("{data_path}ct_markers.rds"))

## top DEGs
ct_markers_level2 <- deg_level2 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()
ct_markers_level4_by2 <- 
  lapply(deg_level4_by2, function(deg_level4_by2x){
    deg_level4_by2x %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 1) %>%
      slice_head(n = 5) %>%
      ungroup()
  })
top5<- rbind(ct_markers_level2[ct_markers_level2$cluster == "B",],
             ct_markers_level4_by2[["B"]],
             ct_markers_level2[ct_markers_level2$cluster == "Plasma",],
             ct_markers_level4_by2[["Plasma"]])
## Pseudo bulk matrix
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
anno_colors <- b_cols
names(anno_colors) <- ann_level
## plot
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

####### Fig. 6G: ct-specific eRegulon #########
## load SCENIC+ object
sp_outpath <- "03_output/07_SCENIC/B/scplus_pipeline/Snakemake/"
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))
scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
## differential analysis
Idents(scAUC_obj) <- scAUC_obj$ann_level4_final
er_gene <- rownames(scAUC_obj)[grep("g\\)", rownames(scAUC_obj))]
diff_auc <- FindAllMarkers(scAUC_obj, 
                           features = er_gene,
                           only.pos = T)
## top markers
sig_diff <- diff_auc %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  top_n(5, -p_val_adj) %>%
  ungroup()
sig_diff <- sig_diff[order(sig_diff$p_val_adj, -abs(sig_diff$avg_log2FC)),]
sig_diff <- sig_diff[!duplicated(sig_diff$gene),]
sig_diff <- sig_diff[order(sig_diff$cluster, sig_diff$p_val_adj, -abs(sig_diff$avg_log2FC)),]
sig_marker <- str_split_i(sig_diff$gene, "-\\(", 1)

## heatdot matrix
sp_outpath <- glue("03_output/07_SCENIC/B/scplus_pipeline/Snakemake/")
heatdot_B <- fread2(glue("{sp_outpath}heatdot_B_ann_level4_final.txt"))[,-1]
#                            "scRNA_counts:ann_level4_final"
colnames(heatdot_B) <- c("celltype", "Gene-based AUC", "eRegulon_name", "Region-based AUC")
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))
e_regulon_name_trans$heatdot_name <- stringr::str_split_i(e_regulon_name_trans$signature, "_\\(", 1)
heatdot_B$signature_name_simplify <- 
  e_regulon_name_trans$signature_name_simplify[match(heatdot_B$eRegulon_name,
                                                     e_regulon_name_trans$heatdot_name)] %>%
  str_split_i(., "_", 1)

## plot
sig_marker <- sig_marker[which(sig_marker %in% heatdot_B$signature_name_simplify)]
heatdot_B_sub <- subset(heatdot_B, signature_name_simplify %in% sig_marker)
heatdot_B_sub$signature_name_simplify <- factor(heatdot_B_sub$signature_name_simplify,
                                                levels = rev(sig_marker))
heatdot_B_sub$celltype <- factor(heatdot_B_sub$celltype,
                                 levels = levels(scAUC_obj$ann_level4_final),
                                 labels = gsub("_", " ", levels(scAUC_obj$ann_level4_final)))
heatdot_B_sub$celltype <- droplevels(heatdot_B_sub$celltype)
heat_dot <- ggplot(heatdot_B_sub, 
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
df_annot <- data.frame(celltype = levels(heatdot_B_sub$celltype))
p_annot = ggplot(df_annot, aes(x = celltype, y = 1, fill = celltype)) +
  scale_fill_manual(values = b_cols) + 
  geom_tile() +
  theme_void() +
  theme(legend.position = "none")

ggsave(glue("{plt_path}heat_dot_B_ann_level4_final.png"),
       wrap_plots(list(p_annot, heat_dot), ncol = 1, heights = c(0.1, 8)),
       height = 9, width = 6, units = "in", dpi = 300)

##### Fig. 6G: pie plot #####
ann_col <-  "ann_level4_final"
group_col <- "Condition"
section_col <- "Section"
##
sc_meta_b$group_comb <- paste0(sc_meta_b[[group_col]], ":", sc_meta_b[[section_col]])
tab_prop <- table(sc_meta_b[[ann_col]], sc_meta_b$group_comb) %>% 
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
  scale_fill_manual(values = b_cols) + 
  theme_void()+
  theme(strip.text = element_text(size = 15, face = "bold"), 
        legend.title = element_text(size = 15, face = "bold"), 
        legend.text = element_text(size = 12))
ggsave(file = glue("{plt_path}pie_plt_B.png"),
       pie_plt,
       width = 6, height = 6,units = "in", dpi = 300)

##### Fig. 6A: volcano plot #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_cr <- glue("{out_path}TI/B/{use_de}/")
ann_col <- "ann_level4_final"
sel_ct <- "B_memory_TOX"
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/volca_plot.R")
## 1. Ileum
deg_df_cr_pb <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[sel_ct]]
deg_df_cr_pb <- subset(dd, !is.na(dd$padj))
vc_plt_pb <- volca_function(deg_df = deg_df_cr_pb,
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
  theme(plot.title = element_text(size = 15, face = "bold"))

ggsave(file = glue("{plt_path}/valcano_{sel_ct}_II_vs_NN_CR.png"),
       vc_plt_pb,
       width = 10, height = 6,units = "in", dpi = 300, limitsize = T)

##### Fig. 6H: scDRS #####
scDRS_out_path <- glue("03_output/01_scDRS/04_downstream/cov/")
##
scDRS_comb_level4 <- fread2(glue("{scDRS_out_path}CD.scdrs_group.ann_level4_final"))
scDRS_comb_level4_B <- subset(scDRS_comb_level4, group %in% levels(sc_meta_b$ann_level4_final))
scDRS_comb_level4_B$section <- "Combined"
scDRS_comb_level4_B$celltype <- factor(scDRS_comb_level4_B$group, 
                                       levels = levels(sc_meta_b$ann_level4_final) %>% rev,
                                       labels = gsub("_", " ", levels(sc_meta_b$ann_level4_final) %>% rev))
##
scDRS_level4 <- fread2(glue("{scDRS_out_path}CD.scdrs_group.sample_ct_level4"))
scDRS_level4$section <- str_split_i(scDRS_level4$group, ":", 1)
scDRS_level4$celltype <- str_split_i(scDRS_level4$group, ":", 2)
scDRS_level4_B <- subset(scDRS_level4, celltype %in% levels(sc_meta_b$ann_level4_final))
scDRS_level4_B$celltype <- factor(scDRS_level4_B$celltype, 
                                  levels = levels(sc_meta_b$ann_level4_final) %>% rev,
                                  labels = gsub("_", " ", levels(sc_meta_b$ann_level4_final) %>% rev))
## 
scDRS_B_plt <- rbind(scDRS_comb_level4_B, scDRS_level4_B)
scDRS_B_plt$section <- factor(scDRS_B_plt$section, 
                                levels = c("Combined", "TI", "RC"),
                                labels = c("Combined", "Ileum", "Colon"))
scDRS_B_plt <- split(scDRS_B_plt, f = scDRS_B_plt$section) %>%
  lapply(., function(scDRS_B_pltx){
    scDRS_B_pltx$FDR <- p.adjust(scDRS_B_pltx$assoc_mcp, method = "BH")
    scDRS_B_pltx$hetero_padj <- p.adjust(scDRS_B_pltx$hetero_mcp, method = "BH")
    return(scDRS_B_pltx)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
scDRS_B_plt$assoc_sig <- scDRS_B_plt$FDR < 0.05
scDRS_B_plt$hetero_sig <- scDRS_B_plt$hetero_padj < 0.05

## plot
scdrs_plt <- ggplot() + 
  geom_tile(data = scDRS_B_plt, 
            aes(x = section, y = celltype, fill = -log10(FDR))) + 
  scale_fill_viridis_c() + 
  geom_tile(data = subset(scDRS_B_plt, assoc_sig), 
            aes(x = section, y = celltype),
            color = "black", fill = NA, size = 1, show.legend = T) +
  geom_text(data = subset(scDRS_B_plt, hetero_sig),
            aes(x = section, y = celltype, label = "\u2716"),
            color = "black", size = 6, show.legend = T) + 
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 11))
#
ggsave(file = glue("{plt_path}CD_sample_ct_B.png"),
       scdrs_plt,
       width = 6, height = 7,units = "in", dpi = 300, limitsize = F)

##### Extended Fig. 9: miloR #####
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

##### Extended Fig. 9: Diff prop ######
sample_col <- "Sample_ID_exp"
block_col <- "Sample_exp"
ann_col <-  "ann_level4_final"
group_col <- "Condition"
section_col <- "Section"

## 0. format proportion table
bulk_meta <- sc_meta_b[!duplicated(sc_meta_b[[sample_col]]),]
bulk_meta$group_comb  <- paste0(bulk_meta[[group_col]], ":", bulk_meta[[section_col]])
sample_list <- split(bulk_meta$Sample_exp, f = bulk_meta$group_comb)
#
tab_prop <- table(sc_meta_b[[ann_col]], sc_meta_b[[sample_col]]) %>% 
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
  # ppairx <- c(p.adjust(ppairx, method = "BH"))
  names(ppairx) <- c(all_ct)
  return(ppairx)
  
}) %>% Reduce("cbind", .) %>% as.data.frame()

# adjust p value for each set of comparison
ppair_prop_adj_df1 <- apply(ppair_prop_df1, 2, function(x){
  p.adjust(x, method = "BH")
}) %>% as.data.frame()
colnames(ppair_prop_df1) <- colnames(ppair_prop_adj_df1) <- names(pair_list)[1:4]

## 2.1 plot for Ileum
box_pair_prop_list1_TI <- lapply(all_ct, function(ctx){
  
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
  
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Condition, y = Proportion, color = Condition)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(ctx) +
    xlab("") + ylab("Proportion (%) in Ileum") + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 5,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})

ggsave(glue("{plt_path}box_pair_level4_Ileum.png"),
       patchwork::wrap_plots(box_pair_prop_list1_TI, nrow = 3),
       height = 12, width = 9, units = "in", dpi = 300)

## 2.2 plot for Colon
box_pair_list1_Colon <- lapply(all_ct, function(ctx){
  
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
  
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Condition, y = Proportion, color = Condition)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(ctx) +
    xlab("") + ylab("Proportion (%) in Colon") + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 5,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})

ggsave(glue("{plt_path}box_pair_level4_Colon.png"),
       patchwork::wrap_plots(box_pair_list1_Colon, nrow = 3),
       height = 12, width = 9, units = "in", dpi = 300)
##### CCC #####
ccc_path_ti <- "03_output/05_Interaction/multinichenetr/TI/B/DESeq2/ann_level4_final/"
ccc_path_cr <- "03_output/05_Interaction/multinichenetr/Colon/B/DESeq2/ann_level4_final/"
system(glue("cp -f {ccc_path_ti}circos_B_to_CD4T_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_B_to_CD4T_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_B_to_CD4T_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_B_to_CD4T_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_B_to_CD4T_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_B_to_CD4T_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_B_to_CD4T_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_B_to_CD4T_NU_vs_NN_Colon_top50.pdf"))
#
system(glue("cp -f {ccc_path_ti}circos_B_to_CD8T_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_B_to_CD8T_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_B_to_CD8T_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_B_to_CD8T_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_B_to_CD8T_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_B_to_CD8T_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_B_to_CD8T_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_B_to_CD8T_NU_vs_NN_Colon_top50.pdf"))
## CD4T to B
ccc_path_ti <- "03_output/05_Interaction/multinichenetr/TI/CD4T/DESeq2/ann_level4_final/"
ccc_path_cr <- "03_output/05_Interaction/multinichenetr/Colon/CD4T/DESeq2/ann_level4_final/"

system(glue("cp -f {ccc_path_ti}circos_CD4T_to_B_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4T_to_B_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_CD4T_to_B_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4T_to_B_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4T_to_B_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4T_to_B_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD4T_to_B_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD4T_to_B_NU_vs_NN_Colon_top50.pdf"))
## CD8T to B
ccc_path_ti <- "03_output/05_Interaction/multinichenetr/TI/CD8T/DESeq2/ann_level4_final/"
ccc_path_cr <- "03_output/05_Interaction/multinichenetr/Colon/CD8T/DESeq2/ann_level4_final/"

system(glue("cp -f {ccc_path_ti}circos_CD8T_to_B_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD8T_to_B_II_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_ti}circos_CD8T_to_B_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD8T_to_B_NU_vs_NN_TI_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD8T_to_B_II_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD8T_to_B_II_vs_NN_Colon_top50.pdf"))
system(glue("cp -f {ccc_path_cr}circos_CD8T_to_B_NU_vs_NN_top50.pdf ",
            "{plt_path}/circos_CD8T_to_B_NU_vs_NN_Colon_top50.pdf"))
