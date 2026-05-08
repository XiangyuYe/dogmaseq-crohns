## module load r/4.5.0
library(Seurat)
library(Signac)
library(bigreadr)
library(dplyr)
library(stringr)
library(glue)
library(harmony)
library(tibble)
library(ggplot2)
library(reshape2)
## set work dir
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
code_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/code/FUNCTION/"
source(glue("{code_path}process/PROCESS_FUN.R"))
##
feature_list <- list(
  "Naive_features_list" = c("CCR7","SELL","LEF1","TCF7"),
  "Memery_features_list" = c("GPR183","S100A4"),
  "Effect_features_list" = c("GZMK","GZMA", "GZMB", "GNLY","NKG7"),
  "Exhaust_features_list" = c("HAVCR2", "KLRK1", "LAG3"),
  "Th1_features_list" = c("CCL5", "CCR5", "IFNG"),
  "Th17_features_list" = c("RORC", "CCR6", "rna_IL17A", "IL17F"),
  "HSP_features_list" = c("HSPA1A", "HSPA1B", "HSPH1"),
  "Tfh_features_list" = c("IL21", "BCL6", "CXCR5", "ICOS", "CXCL13"),
  "Treg_features_list" = c("FOXP3", "IKZF2", "IKZF1"),
  "IL_features_list" = c("IL4", "IL5", "IL9", "IL13", "IL2", "IL2RA", "IL2RB", "IL12RB2"),
  "TF_features_list" = c("CREM", "EGR1", "EGR2", "EGR3", "RUNX1", "RUNX2", "RUNX3",
                         "STAT1", "STAT3", "STAT4", "STAT5A", "STAT5B")
)
atac_mk <- c("IFNG", "TBX21", "CCL5", "IL17A", "RORC", "RORA", "IKZF2", "FOXP3")
sel_motif <- c("MA0690.1", "MA1151.1", "MA0071.1", "MA0072.1")
names_motif <- c("TBX21", "RORC", "RORA.1", "RORA.2")
seed_use <- 20250528

#### load data #####
data_path <- "03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/"
out_path <- "03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/"
sc_obj <- readRDS(glue("{data_path}scWNN_obj.rds"))
##
sc_meta_ref <- readRDS("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/03_clustering_test/CD4_CD103_TRM/test302015/sc_meta.rds")
sc_obj$ann_ref <- sc_meta_ref[colnames(sc_obj) %>% 
                                gsub("_DOGMAseq\\-", "", .) %>%
                                gsub("Duerr_", "", .), 
                              "ann_level5_final"]
##### Reclustering in Th17 to level 5 #####
sc_obj_sub <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_final == "CD4_CD103_TRM_Th17"])
DefaultAssay(sc_obj_sub) <- "RNA"
sc_obj_sub <- NormalizeData(sc_obj_sub)
#
sc_obj_sub <- scsub.renorm(sc_obj = sc_obj_sub,
                            do_ADT = F,
                            do_ATAC = T,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")
sc_obj_sub <- dr.cl.wnn(sc_obj = sc_obj_sub,
                         redc_list = list("harmony_SCT", "harmony_lsi"),
                         dim_list = list(1:15, 1:25),
                         k_nn = 20,
                         prune_SNN = 1/20,
                         n_iter = 300,
                         res = 0.8,
                         run_umap = T,
                         n_neig = 30L,
                         n_epochs = 500,
                         neg_rate = 30L,
                         sprd = 0.5,
                         min_dist = 0.4,
                         seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_th17_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        cols = cols4all::c4a("rainbow", 14),
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_th17.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_th17.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "rna_CD79A", 
                         "rna_CD4", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
                         "rna_GZMK", "rna_GNLY", "rna_KLRD1", "rna_KLRG1",
                         "rna_IL7R", "rna_IL26", "rna_IL2", "rna_IKZF2", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1",
                         "chromvar_MA0477.2", "chromvar_MA0476.1", "chromvar_MA0740.1"), 
            raster = F,
            ncol = 6, 
            min.cutoff = "q5", 
            max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_th17.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_th17.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub$ann_level5_temp <- paste0("ct", sc_obj_sub$seurat_clusters)

##### save temp data #####
sc_meta_sub <- cbind(sc_obj_sub@meta.data,
                     Embeddings(sc_obj_sub, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta_sub, file = glue("{out_path}sc_meta_th17_temp.rds"))
fwrite2(sc_meta_sub, file = glue("{out_path}sc_meta_th17_temp.txt"), row.names = T)
##
saveRDS(sc_meta_sub@reductions, file = glue("{out_path}reducWNN_test_th17.rds"))
saveRDS(sc_meta_sub@commands, file = glue("{out_path}cmdWNN_test_th17.rds"))
##
sc_meta_sub_test <- sc_meta_sub
DefaultAssay(sc_meta_sub_test) <- "ADT"
sc_meta_sub_test[["ATAC"]] <- NULL
sc_meta_sub_test[["RNA"]] <- NULL
sc_meta_sub_test[["SCT"]] <- NULL
saveRDS(sc_meta_sub_test, file = glue("{out_path}scWNN_obj_th17_drcl.rds"))
saveRDS(sc_obj_sub, file = glue("{out_path}scWNN_obj_th17.rds"))

##### bulk heatmap #####
exclude_gene <- readRDS("03_output/02_clean/exclude_gene.rds")
DefaultAssay(sc_obj_sub) <- "RNA"
Idents(sc_obj_sub) <- sc_obj_sub$seurat_clusters
ct_markers <- FindAllMarkers(sc_obj_sub, 
                             features = setdiff(rownames(sc_obj_sub), exclude_gene), 
                             only.pos = T)
##
ct_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
##
sc_obj_bulk <- Seurat::AverageExpression(sc_obj_sub,
                                         features = unique(top5$gene),
                                         group.by = "seurat_clusters",
                                         assays = "RNA",
                                         return.seurat = T,
                                         layer = "data")
tiff(glue("{out_path}pseudo_heat_th17.tiff"),
     height = 7, width = 5, units = "in", res = 300, compression = "lzw")
pheatmap::pheatmap(sc_obj_bulk[["RNA"]]$data[unique(top5$gene), ],
                   scale = "row",
                   border_color = NA,
                   fontsize = 8,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()

##### final anotation #####
## ct-specific eRegulon
sp_outpath <- "03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/"
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))

scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
scAUC_obj_trm <- subset(scAUC_obj, cells = colnames(sc_obj))
scAUC_obj_trm@meta.data <- sc_obj@meta.data

Idents(scAUC_obj_trm) <- scAUC_obj_trm$ann_level5_temp
er_gene <- rownames(scAUC_obj_trm)[grep("g\\)", rownames(scAUC_obj_trm))]
diff_auc <- FindAllMarkers(scAUC_obj_trm, 
                           features = er_gene,
                           only.pos = T)
saveRDS(diff_auc, file = glue("{data_path}diff_auc_treg.rds"))
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
sp_outpath <- glue("03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/")
heatdot_CD103trm <- fread2(glue("{sp_outpath}heatdot_CD103trm_ann_level5_temp.txt"))[,-1]
#                            "scRNA_counts:ann_level5_temp"
colnames(heatdot_CD103trm) <- c("celltype", "Gene-based AUC", "eRegulon_name", "Region-based AUC")
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))
e_regulon_name_trans$heatdot_name <- stringr::str_split_i(e_regulon_name_trans$signature, "_\\(", 1)
heatdot_CD103trm$signature_name_simplify <- 
  e_regulon_name_trans$signature_name_simplify[match(heatdot_CD103trm$eRegulon_name,
                                                     e_regulon_name_trans$heatdot_name)] %>%
  str_split_i(., "_", 1)

##
sig_marker <- sig_marker[which(sig_marker %in% heatdot_CD103trm$signature_name_simplify)]
heatdot_CD103trm_sub <- subset(heatdot_CD103trm, signature_name_simplify %in% sig_marker)
heatdot_CD103trm_sub$signature_name_simplify <- factor(heatdot_CD103trm_sub$signature_name_simplify,
                                                       levels = rev(sig_marker))
heatdot_CD103trm_sub$celltype <- factor(heatdot_CD103trm_sub$celltype,
                                        levels = unique(sc_obj$ann_level5_temp) %>% sort,
                                        labels = gsub("_", " ", levels(scAUC_obj_trm$ann_level5_temp)))
heat_dot <- ggplot(heatdot_CD103trm_sub, 
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
df_annot <- data.frame(celltype = unique(sc_obj$ann_level5_temp) %>% sort)
p_annot = ggplot(df_annot, aes(x = celltype, y = 1, fill = celltype)) +
  scale_fill_manual(values = cd4cd103trm_cols) + 
  geom_tile() +
  theme_void() +
  theme(legend.position = "top")

ggsave(glue("{sp_outpath}heat_dot_CD4CD103TRM_ann_level5_temp.png"),
       wrap_plots(list(p_annot, heat_dot), ncol = 1, heights = c(1, 8)),
       height = 17, width = 9, units = "in", dpi = 300)

##
sc_obj_sub$ann_level5_final <- recode(sc_obj_sub$ann_level5_temp,
                                       "ct0" = "CD4_CD103_TRM_Th17_RUNX1",
                                       "ct1" = "CD4_CD103_TRM_Th17_HSP",
                                       "ct2" = "CD4_CD103_TRM_Th17_NFKB",
                                       "ct3" = "CD4_CD103_TRM_Th17",
                                       "ct4" = "CD4_CD103_TRM_Th17_CREM",
                                       "ct5" = "CD4_CD103_TRM_Th17_IKZF1",
                                       "ct6" = "CD4_CD103_TRM_Th17_EGR2",
                                       "ct7" = "CD4_CD103_TRM_Th17_EGR1",
                                       "ct8" = "CD4_CD103_TRM_Th17_IG")
sc_obj_sub$ann_level5_final <- factor(sc_obj_sub$ann_level5_final,
                                      levels = c("CD4_CD103_TRM_Th17",
                                                 "CD4_CD103_TRM_Th17_CREM",
                                                 "CD4_CD103_TRM_Th17_EGR1",
                                                 "CD4_CD103_TRM_Th17_EGR2",
                                                 "CD4_CD103_TRM_Th17_HSP",
                                                 "CD4_CD103_TRM_Th17_IG",
                                                 "CD4_CD103_TRM_Th17_IKZF1",
                                                 "CD4_CD103_TRM_Th17_NFKB",
                                                 "CD4_CD103_TRM_Th17_RUNX1"))

##### combine with all CD103 TRM #####
sc_obj$ann_level5_final <- sc_obj$ann_level4_final
sc_obj$ann_level5_final[match(colnames(sc_obj_sub), colnames(sc_obj))] <- sc_obj_sub$ann_level5_final %>% as.character()
sc_obj$ann_level5_final <- factor(sc_obj$ann_level5_final,
                                  levels = c("CD4_CD103_TRM_Th1", 
                                             levels(sc_obj_sub$ann_level5_final),
                                             "CD4_CD103_TRM_Th17.1"))
## 
tiff(file = glue("{out_path}UMAP_level5_final.tiff"),
     width = 9, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        group.by = "ann_level5_final",
        reduction = "wnn.umap", 
        label = T,
        raster = F)
dev.off()

##### save final data #####
##
sc_meta_sub <- cbind(sc_obj_sub@meta.data,
                     Embeddings(sc_obj_sub, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta_sub, file = glue("{out_path}sc_meta_th17.rds"))
fwrite2(sc_meta_sub, file = glue("{out_path}sc_meta_th17.txt"), row.names = T)
##
saveRDS(sc_obj_sub@reductions, file = glue("{out_path}reducWNN_test_th17.rds"))
saveRDS(sc_obj_sub@commands, file = glue("{out_path}cmdWNN_test_th17.rds"))
##
sc_obj_sub_test <- sc_obj_sub
DefaultAssay(sc_obj_sub_test) <- "chromvar"
sc_obj_sub_test[["peaks"]] <- NULL
sc_obj_sub_test[["RNA"]] <- NULL
sc_obj_sub_test[["SCT"]] <- NULL
saveRDS(sc_obj_sub_test, file = glue("{out_path}scWNN_obj_th17_drcl.rds"))
saveRDS(sc_obj_sub, file = glue("{out_path}scWNN_obj_th17.rds"))

##
sc_meta <- cbind(sc_obj@meta.data,
                 Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta.rds"))
fwrite2(sc_meta, file = glue("{out_path}sc_meta.txt"), row.names = T)
#
saveRDS(sc_obj@reductions, file = glue("{out_path}reducWNN_test.rds"))
saveRDS(sc_obj@commands, file = glue("{out_path}cmdWNN_test.rds"))
#
sc_obj_test <- sc_obj
DefaultAssay(sc_obj_test) <- "chromvar"
sc_obj_test[["peaks"]] <- NULL
sc_obj_test[["RNA"]] <- NULL
sc_obj_test[["SCT"]] <- NULL
saveRDS(sc_obj_test, file = glue("{out_path}scWNN_obj_drcl.rds"))
saveRDS(sc_obj, file = glue("{out_path}scWNN_obj.rds"))
