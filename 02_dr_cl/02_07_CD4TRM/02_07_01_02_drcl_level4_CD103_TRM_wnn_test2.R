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
  "T_features_list" = c("CD3D", "CD4", "CD8A", "CD8B","TRDC", "S1PR1", "CX3CR1"),
  "Prolif_features_list" = c("MKI67","TYMS","PCNA", "TOP2A"),
  "MAIT_features_list" = c("SLC4A10", "KIT"),
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

##### load ADT, RNA, and ATAC data #####
data_path <- "03_output/03_clustering/CD4_CD103_TRM/"
out_path <- "03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/test2/"
scRNA_obj <- readRDS(glue("{data_path}scRNA_obj_test.rds"))
sc_obj <- scRNA_obj
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)
sc_obj$wnnUMAP_1 <- NULL
sc_obj$wnnUMAP_2 <- NULL
#
cutoff_q <- "q70"
scATAC_obj <- readRDS(glue("{data_path}scATAC_obj_test_{cutoff_q}.rds"))
sc_obj[["peaks"]] <- scATAC_obj[["peaks"]]
sc_obj@reductions$harmony_lsi <- scATAC_obj@reductions$harmony_lsi
sc_obj@reductions$umap_lsi <- scATAC_obj@reductions$umap_lsi
#
sc_chromvar <- readRDS(glue("03_output/03_clustering/CD4T/scchromvar_assay_{cutoff_q}.rds")) %>%
  CreateSeuratObject(assay = "chromvar",
                     min.cells = 0,
                     min.features = 0)
sc_obj[["chromvar"]] <- subset(sc_chromvar, cells = colnames(sc_obj))[["chromvar"]]
#
sc_act <- readRDS("03_output/03_clustering/scATAC_act_mat_immune.rds") %>%
  CreateSeuratObject(assay = "activity",
                     min.cells = 0,
                     min.features = 0)
sc_obj[["activity"]] <- subset(sc_act, cells = colnames(sc_obj))[["activity"]]
DefaultAssay(sc_obj) <- "activity"
sc_obj <- NormalizeData(sc_obj)

##
sc_meta_ref <- readRDS("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/03_clustering_test/CD4_CD103_TRM/test302015/sc_meta.rds")
sc_obj$ann_ref <- sc_meta_ref[colnames(sc_obj) %>% 
                                gsub("_DOGMAseq\\-", "", .) %>%
                                gsub("Duerr_", "", .), 
                              "ann_level5_final"]

##### DRCL #####
sc_obj <- dr.cl.wnn(sc_obj = sc_obj,
                    redc_list = list("harmony_SCT", "harmony_lsi"),
                    dim_list = list(1:15, 1:25),
                    k_nn = 20,
                    prune_SNN = 1/20,
                    n_iter = 300,
                    res = 1,
                    cl_method = 1,
                    run_umap = T,
                    n_neig = 30L,
                    n_epochs = 500,
                    neg_rate = 30L,
                    sprd = 0.7,
                    min_dist = 0.4,
                    seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_test.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F, 
        label = T) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_test_ann_ref.tiff"),
     width = 9, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()

## 
tiff(glue("{out_path}ft.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj, 
            reduction = "wnn.umap",
            features = c("rna_EGR1", "rna_EGR2", "rna_IKZF1", "rna_IL2", "rna_IL17A",
                         "rna_CCL5", "rna_CCR6", "rna_KIT", "rna_CREM",
                         "chromvar_MA0690.1", "chromvar_MA1151.1", "activity_TBX21", "activity_IFNG", "activity_CCL5", "activity_RORC"), 
            raster = F,
            ncol = 5, 
            min.cutoff = "q1", 
            max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_test.tiff"), 
     height = 10, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", 
        col.min = -2, 
        col.max = 2, 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_test.tiff"), 
     height = 10, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()
##
dot_act <- DotPlot(sc_obj, 
                   assay = "activity", 
                   col.min = -2,
                   col.max = 2,
                   features = atac_mk,
                   group.by = "seurat_clusters",
                   cols = c("lightgrey", "brown")) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_act_test.tiff"), 
     height = 8, width = 6, units = "in", res = 300, compression = "lzw")
print(dot_act)
dev.off()

## temp anotation
sc_obj$ann_level4_temp1 <- "CD4_CD103_TRM_Th17"
sc_obj$ann_level4_temp1[sc_obj$wsnn_res.1 %in% c(4)] <- "CD4_CD103_TRM_Th1"
sc_obj$ann_level4_temp1[sc_obj$wsnn_res.1 %in% c(9, 11)] <- "CD4_CD103_TRM_Th17.1"
##
tiff(file = glue("{out_path}UMAP_temp1_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1",
        raster = F,
        label = T) + NoLegend()
dev.off()
##### round 1 #####
##### round 1_1 #####
sc_obj_sub1_1 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp1 == "CD4_CD103_TRM_Th17"])
DefaultAssay(sc_obj_sub1_1) <- "RNA"
sc_obj_sub1_1 <- NormalizeData(sc_obj_sub1_1)
#
sc_obj_sub1_1 <- scsub.renorm(sc_obj = sc_obj_sub1_1,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub1_1 <- dr.cl.wnn(sc_obj = sc_obj_sub1_1,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 100,
                           neg_rate = 5L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_1.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_1, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub1_1.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_1, 
            reduction = "wnn.umap",
            features = c("rna_EGR1", "rna_EGR2", "rna_IKZF1", "rna_IL2", "rna_IL17A",
                         "rna_CCL5", "rna_CCR6", "rna_KIT", "rna_CREM",
                         "chromvar_MA0690.1", "chromvar_MA1151.1", "activity_TBX21", "activity_IFNG", "activity_CCL5", "activity_RORC"), 
            raster = F,
            ncol = 5, 
            min.cutoff = "q1", 
            max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_1.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_1, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub1_1, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub1_1.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub1_1$ann_level4_temp1_1 <- "CD4_CD103_TRM_Th17"
sc_obj_sub1_1$ann_level4_temp1_1[sc_obj_sub1_1$wsnn_res.2 %in% c(1, 2, 7, 10, 12, 13, 19)] <- "Mixed"
sc_obj$ann_level4_temp1_1 <- sc_obj$ann_level4_temp1
sc_obj$ann_level4_temp1_1[match(colnames(sc_obj_sub1_1), colnames(sc_obj))] <- sc_obj_sub1_1$ann_level4_temp1_1

##
tiff(file = glue("{out_path}UMAP_sub1_1_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_1, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_1",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_1_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_1",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 1_2 #####
sc_obj_sub1_2 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp1_1 == "Mixed"])
DefaultAssay(sc_obj_sub1_2) <- "RNA"
sc_obj_sub1_2 <- NormalizeData(sc_obj_sub1_2)
#
sc_obj_sub1_2 <- scsub.renorm(sc_obj = sc_obj_sub1_2,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub1_2 <- dr.cl.wnn(sc_obj = sc_obj_sub1_2,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_2.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_2, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub1_2.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_2, 
            reduction = "wnn.umap",
            features = c("rna_EGR1", "rna_EGR2", "rna_IKZF1", "rna_IL2", "rna_IL17A",
                         "rna_CCL5", "rna_CCR6", "rna_KIT", "rna_CREM",
                         "chromvar_MA0690.1", "chromvar_MA1151.1", "activity_TBX21", "activity_IFNG", "activity_CCL5", "activity_RORC"), 
            raster = F,
            ncol = 5, 
            min.cutoff = "q1", 
            max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_2.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_2, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub1_2, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub1_2.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub1_2$ann_level4_temp1_2 <- "CD4_CD103_TRM_Th17"
sc_obj_sub1_2$ann_level4_temp1_2[sc_obj_sub1_2$wsnn_res.2 %in% c(10, 15)] <- "CD4_CD103_TRM_Th17.1"
sc_obj_sub1_2$ann_level4_temp1_2[sc_obj_sub1_2$wsnn_res.2 %in% c(3:7)] <- "Mixed"
sc_obj$ann_level4_temp1_2 <- sc_obj$ann_level4_temp1_1
sc_obj$ann_level4_temp1_2[match(colnames(sc_obj_sub1_2), colnames(sc_obj))] <- sc_obj_sub1_2$ann_level4_temp1_2

##
tiff(file = glue("{out_path}UMAP_sub1_2_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_2, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_2_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_2",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 1_3 #####
sc_obj_sub1_3 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp1_2 == "Mixed"])
DefaultAssay(sc_obj_sub1_3) <- "RNA"
sc_obj_sub1_3 <- NormalizeData(sc_obj_sub1_3)
#
sc_obj_sub1_3 <- scsub.renorm(sc_obj = sc_obj_sub1_3,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub1_3 <- dr.cl.wnn(sc_obj = sc_obj_sub1_3,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_3.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_3, 
        reduction = "wnn.umap",
        label = T)
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub1_3_ann_ref.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_3, 
        reduction = "wnn.umap", 
        group.by = "ann_ref",
        label = T)
dev.off()
##
tiff(glue("{out_path}ft_sub1_3.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_3, 
            reduction = "wnn.umap",
            features = c("rna_EGR1", "rna_EGR2", "rna_IKZF1", "rna_IL2", "rna_IL17A",
                         "rna_CCL5", "rna_CCR6", "rna_KIT", "rna_CREM",
                         "chromvar_MA0690.1", "chromvar_MA1151.1", "activity_TBX21", "activity_IFNG", "activity_CCL5", "activity_RORC"), 
            raster = F,
            ncol = 5, 
            min.cutoff = "q1", 
            max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_3.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_3, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub1_3, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub1_3.tiff"), 
     height = 8, width = 4, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub1_3$ann_level4_temp1_3 <- "Mixed"
sc_obj_sub1_3$ann_level4_temp1_3[sc_obj_sub1_3$wsnn_res.2 %in% c(2, 7, 10, 13, 16, 18, 20)] <- "CD4_CD103_TRM_Th17"
sc_obj_sub1_3$ann_level4_temp1_3[sc_obj_sub1_3$wsnn_res.2 %in% c(19)] <- "CD4_CD103_TRM_Th17.1"
sc_obj$ann_level4_temp1_3 <- sc_obj$ann_level4_temp1_2
sc_obj$ann_level4_temp1_3[match(colnames(sc_obj_sub1_3), colnames(sc_obj))] <- sc_obj_sub1_3$ann_level4_temp1_3

##
tiff(file = glue("{out_path}UMAP_sub1_3_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_3, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_3",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_3_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_3",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 1_4 #####
sc_obj_sub1_4 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp1_3 == "Mixed"])
DefaultAssay(sc_obj_sub1_4) <- "RNA"
sc_obj_sub1_4 <- NormalizeData(sc_obj_sub1_4)
#
sc_obj_sub1_4 <- scsub.renorm(sc_obj = sc_obj_sub1_4,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub1_4 <- dr.cl.wnn(sc_obj = sc_obj_sub1_4,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:10, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_4_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub1_4, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub1_4.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_4, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub1_4.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_4, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub1_4.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_4, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub1_4, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub1_4.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub1_4$ann_level4_temp1_4 <- "Mixed"
sc_obj_sub1_4$ann_level4_temp1_4[sc_obj_sub1_4$wsnn_res.2 %in% c(2, 6, 8, 12, 14, 16)] <- "CD4_CD103_TRM_Th17"
sc_obj_sub1_4$ann_level4_temp1_4[sc_obj_sub1_4$wsnn_res.2 %in% c(15)] <- "CD4_CD103_TRM_Th1"
sc_obj$ann_level4_temp1_4 <- sc_obj$ann_level4_temp1_3
sc_obj$ann_level4_temp1_4[match(colnames(sc_obj_sub1_4), colnames(sc_obj))] <- sc_obj_sub1_4$ann_level4_temp1_4

##
tiff(file = glue("{out_path}UMAP_sub1_4_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_4, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_4",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_4_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_4",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 1_5 #####
sc_obj_sub1_5 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp1_4 == "Mixed"])
DefaultAssay(sc_obj_sub1_5) <- "RNA"
sc_obj_sub1_5 <- NormalizeData(sc_obj_sub1_5)
#
sc_obj_sub1_5 <- scsub.renorm(sc_obj = sc_obj_sub1_5,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub1_5 <- dr.cl.wnn(sc_obj = sc_obj_sub1_5,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:10, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_5.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_5, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub1_5.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_5, 
            reduction = "wnn.umap",
            features = c("rna_EGR1", "rna_EGR2", "rna_IKZF1", "rna_IL2", "rna_IL17A",
                         "rna_CCL5", "rna_CCR6", "rna_KIT", "rna_CREM",
                         "chromvar_MA0690.1", "chromvar_MA1151.1", "activity_TBX21", "activity_IFNG", "activity_CCL5", "activity_RORC"), 
            raster = F,
            ncol = 5, 
            min.cutoff = "q1", 
            max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_5.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_5, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub1_5, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub1_5.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub1_5$ann_level4_temp1_5 <- "CD4_CD103_TRM_Th17"
sc_obj_sub1_5$ann_level4_temp1_5[sc_obj_sub1_5$wsnn_res.2 %in% c(10, 14, 15)] <- "CD4_CD103_TRM_Th17"
sc_obj$ann_level4_temp1_5 <- sc_obj$ann_level4_temp1_4
sc_obj$ann_level4_temp1_5[match(colnames(sc_obj_sub1_5), colnames(sc_obj))] <- sc_obj_sub1_5$ann_level4_temp1_5

##
tiff(file = glue("{out_path}UMAP_sub1_5_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_5, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_5",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_5_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1_5",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2 #####
##### round 2_1 #####
sel_ct <- c("CD4_CD103_TRM_Th17.1")
sc_obj_sub2_1 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp1_5 %in% sel_ct])
DefaultAssay(sc_obj_sub2_1) <- "RNA"
sc_obj_sub2_1 <- NormalizeData(sc_obj_sub2_1)
#
sc_obj_sub2_1 <- scsub.renorm(sc_obj = sc_obj_sub2_1,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_1 <- dr.cl.wnn(sc_obj = sc_obj_sub2_1,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_1.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_1, 
        reduction = "wnn.umap",
        label = T)
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub2_1_ann_ref.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_1, 
        reduction = "wnn.umap", 
        group.by = "ann_ref",
        label = T)
dev.off()
##
tiff(glue("{out_path}ft_sub2_1.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_1, 
            reduction = "wnn.umap",
            features = c("rna_EGR1", "rna_EGR2", "rna_IKZF1", "rna_IL2", "rna_IL17A",
                         "rna_CCL5", "rna_CCR6", "rna_KIT", "rna_CREM",
                         "chromvar_MA0690.1", "chromvar_MA1151.1", "activity_TBX21", "activity_IFNG", "activity_CCL5", "activity_RORC"), 
            raster = F,
            ncol = 5, 
            min.cutoff = "q1", 
            max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_1.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_1, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub2_1, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub2_1.tiff"), 
     height = 8, width = 4, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()
##
dot_act <- DotPlot(sc_obj_sub2_1, 
                   assay = "activity", 
                   col.min = -2,
                   col.max = 2,
                   features = atac_mk,
                   group.by = "seurat_clusters",
                   cols = c("lightgrey", "brown")) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_act_sub2_1.tiff"), 
     height = 8, width = 6, units = "in", res = 300, compression = "lzw")
print(dot_act)
dev.off()

## temp anotation
sc_obj_sub2_1$ann_level4_temp2_1 <- sc_obj_sub2_1$ann_level4_temp1_5
# sc_obj_sub2_1$ann_level4_temp2_1[sc_obj_sub2_1$wsnn_res.2 %in% c(0, 1, 2, 4, 5, 9, 10, 13:15, 17, 22, 24, 25)] <- "CD4_CD103_TRM_Th17.1"
sc_obj$ann_level4_temp2_1 <- sc_obj$ann_level4_temp1_5
sc_obj$ann_level4_temp2_1[match(colnames(sc_obj_sub2_1), colnames(sc_obj))] <- sc_obj_sub2_1$ann_level4_temp2_1

##
tiff(file = glue("{out_path}UMAP_sub2_1_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_1, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_1",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_1_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_1",
        raster = F,
        label = T) + NoLegend()
dev.off()


##### round 2_2 #####
sel_ct <- c("CD4_CD103_TRM_Th1", "CD4_CD103_TRM_Th17")
sc_obj_sub2_2 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp2_1 %in% sel_ct])
DefaultAssay(sc_obj_sub2_2) <- "RNA"
sc_obj_sub2_2 <- NormalizeData(sc_obj_sub2_2)
#
sc_obj_sub2_2 <- scsub.renorm(sc_obj = sc_obj_sub2_2,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_2 <- dr.cl.wnn(sc_obj = sc_obj_sub2_2,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_2_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub2_2, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub2_2.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_2, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub2_2.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_2, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub2_2.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_2, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub2_2, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub2_2.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub2_2$ann_level4_temp2_2 <- sc_obj_sub2_2$ann_level4_temp2_1
# sc_obj_sub2_2$ann_level4_temp2_2[sc_obj_sub2_2$wsnn_res.2 %in% c(5, 14)] <- "CD4_CD103_TRM_Th1"
sc_obj_sub2_2$ann_level4_temp2_2[sc_obj_sub2_2$wsnn_res.2 %in% c(4, 13, 19)] <- "CD4_CD103_TRM_Th17.1"
sc_obj$ann_level4_temp2_2 <- sc_obj$ann_level4_temp2_1
sc_obj$ann_level4_temp2_2[match(colnames(sc_obj_sub2_2), colnames(sc_obj))] <- sc_obj_sub2_2$ann_level4_temp2_2

##
tiff(file = glue("{out_path}UMAP_sub2_2_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_2, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_2_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_2",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2_3 #####
sel_ct <- c("CD4_CD103_TRM_Th17", "CD4_CD103_TRM_Th17.1")
sc_obj_sub2_3 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp2_2 %in% sel_ct])
DefaultAssay(sc_obj_sub2_3) <- "RNA"
sc_obj_sub2_3 <- NormalizeData(sc_obj_sub2_3)
#
sc_obj_sub2_3 <- scsub.renorm(sc_obj = sc_obj_sub2_3,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_3 <- dr.cl.wnn(sc_obj = sc_obj_sub2_3,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:15, 1:25),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_3_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub2_3, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub2_3.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_3, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub2_3.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_3, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub2_3.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_3, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub2_3, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub2_3.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub2_3$ann_level4_temp2_3 <- "CD4_CD103_TRM_Th17"
sc_obj_sub2_3$ann_level4_temp2_3[sc_obj_sub2_3$wsnn_res.2 %in% c(10, 21)] <- "CD4_CD103_TRM_Th17.1"
sc_obj_sub2_3$ann_level4_temp2_3[sc_obj_sub2_3$wsnn_res.2 %in% c(11, 15, 16)] <- "Mixed"
sc_obj$ann_level4_temp2_3 <- sc_obj$ann_level4_temp2_2
sc_obj$ann_level4_temp2_3[match(colnames(sc_obj_sub2_3), colnames(sc_obj))] <- sc_obj_sub2_3$ann_level4_temp2_3

##
tiff(file = glue("{out_path}UMAP_sub2_3_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_3, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_3",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_3_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_3",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2_4 #####
sel_ct <- c("Mixed")
sc_obj_sub2_4 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp2_3 %in% sel_ct])
DefaultAssay(sc_obj_sub2_4) <- "RNA"
sc_obj_sub2_4 <- NormalizeData(sc_obj_sub2_4)
#
sc_obj_sub2_4 <- scsub.renorm(sc_obj = sc_obj_sub2_4,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_4 <- dr.cl.wnn(sc_obj = sc_obj_sub2_4,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 1,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_4_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub2_4, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub2_4.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_4, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub2_4.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_4, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub2_4.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_4, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub2_4, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub2_4.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub2_4$ann_level4_temp2_4 <- "Mixed"
sc_obj_sub2_4$ann_level4_temp2_4[sc_obj_sub2_4$wsnn_res.1 %in% c(6, 7, 10)] <- "CD4_CD103_TRM_Th17.1"
sc_obj_sub2_4$ann_level4_temp2_4[sc_obj_sub2_4$wsnn_res.1 %in% c(0, 4, 12)] <- "CD4_CD103_TRM_Th17"
sc_obj$ann_level4_temp2_4 <- sc_obj$ann_level4_temp2_3
sc_obj$ann_level4_temp2_4[match(colnames(sc_obj_sub2_4), colnames(sc_obj))] <- sc_obj_sub2_4$ann_level4_temp2_4

##
tiff(file = glue("{out_path}UMAP_sub2_4_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_4, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_4",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_4_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_4",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2_5 #####
sel_ct <- c("Mixed")
sc_obj_sub2_5 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp2_4 %in% sel_ct])
DefaultAssay(sc_obj_sub2_5) <- "RNA"
sc_obj_sub2_5 <- NormalizeData(sc_obj_sub2_5)
#
sc_obj_sub2_5 <- scsub.renorm(sc_obj = sc_obj_sub2_5,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_5 <- dr.cl.wnn(sc_obj = sc_obj_sub2_5,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:10, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 1,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_5_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub2_5, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub2_5.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_5, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub2_5.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_5, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub2_5.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_5, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub2_5, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub2_5.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub2_5$ann_level4_temp2_5 <- "CD4_CD103_TRM_Th17"
sc_obj_sub2_5$ann_level4_temp2_5[sc_obj_sub2_5$wsnn_res.1 %in% c(3:5, 11)] <- "CD4_CD103_TRM_Th17.1"
sc_obj$ann_level4_temp2_5 <- sc_obj$ann_level4_temp2_4
sc_obj$ann_level4_temp2_5[match(colnames(sc_obj_sub2_5), colnames(sc_obj))] <- sc_obj_sub2_5$ann_level4_temp2_5

##
tiff(file = glue("{out_path}UMAP_sub2_5_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_5, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_5",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_5_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2_5",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 3 #####
##### round 3_1 #####
sel_ct <- c("CD4_CD103_TRM_Th1")
sc_obj_sub3_1 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp2_3 %in% sel_ct])
DefaultAssay(sc_obj_sub3_1) <- "RNA"
sc_obj_sub3_1 <- NormalizeData(sc_obj_sub3_1)
#
sc_obj_sub3_1 <- scsub.renorm(sc_obj = sc_obj_sub3_1,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub3_1 <- dr.cl.wnn(sc_obj = sc_obj_sub3_1,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_1.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_1, 
        reduction = "wnn.umap",
        label = T)
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3_1_ann_ref.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_1, 
        reduction = "wnn.umap", 
        group.by = "ann_ref",
        label = T)
dev.off()
##
tiff(glue("{out_path}ft_sub3_1.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_1, 
            reduction = "wnn.umap",
            features = c("rna_EGR1", "rna_EGR2", "rna_IKZF1", "rna_IL2", "rna_IL17A",
                         "rna_CCL5", "rna_CCR6", "rna_KIT", "rna_CREM",
                         "chromvar_MA0690.1", "chromvar_MA1151.1", "activity_TBX21", "activity_IFNG", "activity_CCL5", "activity_RORC"), 
            raster = F,
            ncol = 5, 
            min.cutoff = "q1", 
            max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub3_1.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_1, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub3_1, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub3_1.tiff"), 
     height = 8, width = 4, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()
##
dot_act <- DotPlot(sc_obj_sub3_1, 
                   assay = "activity", 
                   col.min = -2,
                   col.max = 2,
                   features = atac_mk,
                   group.by = "seurat_clusters",
                   cols = c("lightgrey", "brown")) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_act_sub3_1.tiff"), 
     height = 8, width = 6, units = "in", res = 300, compression = "lzw")
print(dot_act)
dev.off()

## temp anotation
sc_obj_sub3_1$ann_level4_temp3_1 <- "CD4_CD103_TRM_Th1"
sc_obj_sub3_1$ann_level4_temp3_1[sc_obj_sub3_1$wsnn_res.2 %in% c(2, 5, 11, 13, 14, 16 ,22)] <- "CD4_CD103_TRM_Th17.1"
sc_obj$ann_level4_temp3_1 <- sc_obj$ann_level4_temp2_5
sc_obj$ann_level4_temp3_1[match(colnames(sc_obj_sub3_1), colnames(sc_obj))] <- sc_obj_sub3_1$ann_level4_temp3_1

# ## see in all-cell scale
# sc_obj_sub_test <- sc_obj
# sc_obj_sub_test$ann_test <- sc_obj_sub_test$ann_level4_temp2_5
# sc_obj_sub_test$ann_test[match(colnames(sc_obj_sub3_1),
#                                colnames(sc_obj_sub_test))] <- paste0("a", sc_obj_sub3_1$seurat_clusters)
# dot_motif1 <- DotPlot(sc_obj_sub_test, 
#                       assay = "chromvar", 
#                       col.min = -2,
#                       col.max = 2,
#                       features = sel_motif,
#                       group.by = "ann_test",
#                       cols = c("lightgrey", "brown")) + 
#   scale_x_discrete(labels = names_motif) +
#   theme(axis.text.x = element_text(angle = 90))
# tiff(glue("{out_path}dot_raw_motif_sub3_1_test.tiff"), 
#      height = 10, width = 8, units = "in", res = 300, compression = "lzw")
# print(dot_motif1)
# dev.off()

##
tiff(file = glue("{out_path}UMAP_sub3_1_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_1, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_1",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp3_1_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_1",
        raster = F,
        label = T) + NoLegend()
dev.off()


##### round 3_2 #####
sel_ct <- c("CD4_CD103_TRM_Th17.1")
sc_obj_sub3_2 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp3_1 %in% sel_ct])
DefaultAssay(sc_obj_sub3_2) <- "RNA"
sc_obj_sub3_2 <- NormalizeData(sc_obj_sub3_2)
#
sc_obj_sub3_2 <- scsub.renorm(sc_obj = sc_obj_sub3_2,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub3_2 <- dr.cl.wnn(sc_obj = sc_obj_sub3_2,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_2_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub3_2, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3_2.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_2, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub3_2.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_2, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub3_2.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_2, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub3_2, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub3_2.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub3_2$ann_level4_temp3_2 <- "CD4_CD103_TRM_Th17.1"
sc_obj_sub3_2$ann_level4_temp3_2[sc_obj_sub3_2$wsnn_res.2 %in% c(4, 6, 15, 16)] <- "CD4_CD103_TRM_Th17"
sc_obj$ann_level4_temp3_2 <- sc_obj$ann_level4_temp3_1
sc_obj$ann_level4_temp3_2[match(colnames(sc_obj_sub3_2), colnames(sc_obj))] <- sc_obj_sub3_2$ann_level4_temp3_2

# ## see in all-cell scale
# sc_obj_sub3_2_test <- sc_obj
# sc_obj_sub3_2_test$ann_test <- sc_obj_sub3_2_test$ann_level4_temp3_1
# sc_obj_sub3_2_test$ann_test[match(colnames(sc_obj_sub3_2),
#                                colnames(sc_obj_sub3_2_test))] <- paste0("a", sc_obj_sub3_2$seurat_clusters)
# dot_motif1 <- DotPlot(sc_obj_sub3_2_test,
#                       assay = "chromvar",
#                       col.min = -2,
#                       col.max = 2,
#                       features = sel_motif,
#                       group.by = "ann_test",
#                       cols = c("lightgrey", "brown")) +
#   scale_x_discrete(labels = names_motif) +
#   theme(axis.text.x = element_text(angle = 90))
# tiff(glue("{out_path}dot_raw_motif_sub3_2_test.tiff"),
#      height = 10, width = 8, units = "in", res = 300, compression = "lzw")
# print(dot_motif1)
# dev.off()

##
tiff(file = glue("{out_path}UMAP_sub3_2_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_2, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp3_2_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_2",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 3_3 #####
sel_ct <- c("CD4_CD103_TRM_Th17")
sc_obj_sub3_3 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp3_2 %in% sel_ct])
DefaultAssay(sc_obj_sub3_3) <- "RNA"
sc_obj_sub3_3 <- NormalizeData(sc_obj_sub3_3)
#
sc_obj_sub3_3 <- scsub.renorm(sc_obj = sc_obj_sub3_3,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub3_3 <- dr.cl.wnn(sc_obj = sc_obj_sub3_3,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_3_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub3_3, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3_3.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_3, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub3_3.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_3, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub3_3.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_3, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub3_3, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub3_3.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

# ## see in all-cell scale
# sc_obj_sub_test <- sc_obj
# sc_obj_sub_test$ann_test <- sc_obj_sub_test$ann_level4_temp3_2
# sc_obj_sub_test$ann_test[match(colnames(sc_obj_sub3_3),
#                                colnames(sc_obj_sub_test))] <- paste0("a", sc_obj_sub3_3$seurat_clusters)
# dot_motif1 <- DotPlot(sc_obj_sub_test,
#                       assay = "chromvar",
#                       col.min = -2,
#                       col.max = 2,
#                       features = sel_motif,
#                       group.by = "ann_test",
#                       cols = c("lightgrey", "brown")) +
#   scale_x_discrete(labels = names_motif) +
#   theme(axis.text.x = element_text(angle = 90))
# tiff(glue("{out_path}dot_raw_motif_sub3_3_test.tiff"),
#      height = 10, width = 8, units = "in", res = 300, compression = "lzw")
# print(dot_motif1)
# dev.off()

## temp anotation
sc_obj_sub3_3$ann_level4_temp3_3 <- "CD4_CD103_TRM_Th17"
sc_obj_sub3_3$ann_level4_temp3_3[sc_obj_sub3_3$wsnn_res.2 %in% c(8, 11)] <- "CD4_CD103_TRM_Th17.1"
sc_obj$ann_level4_temp3_3 <- sc_obj$ann_level4_temp3_2
sc_obj$ann_level4_temp3_3[match(colnames(sc_obj_sub3_3), colnames(sc_obj))] <- sc_obj_sub3_3$ann_level4_temp3_3

##
tiff(file = glue("{out_path}UMAP_sub3_3_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_3, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_3",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp3_3_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_3",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 3_4 #####
sel_ct <- c("CD4_CD103_TRM_Th17.1")
sc_obj_sub3_4 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp3_3 %in% sel_ct])
DefaultAssay(sc_obj_sub3_4) <- "RNA"
sc_obj_sub3_4 <- NormalizeData(sc_obj_sub3_4)
#
sc_obj_sub3_4 <- scsub.renorm(sc_obj = sc_obj_sub3_4,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub3_4 <- dr.cl.wnn(sc_obj = sc_obj_sub3_4,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_4_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub3_4, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3_4.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_4, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub3_4.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_4, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub3_4.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_4, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub3_4, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub3_4.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## temp anotation
sc_obj_sub3_4$ann_level4_temp3_4 <- "CD4_CD103_TRM_Th17.1"
sc_obj_sub3_4$ann_level4_temp3_4[sc_obj_sub3_4$wsnn_res.2 %in% c(0, 16)] <- "CD4_CD103_TRM_Th17"
sc_obj$ann_level4_temp3_4 <- sc_obj$ann_level4_temp3_3
sc_obj$ann_level4_temp3_4[match(colnames(sc_obj_sub3_4), colnames(sc_obj))] <- sc_obj_sub3_4$ann_level4_temp3_4

## see in all-cell scale
sc_obj_sub_test <- sc_obj
sc_obj_sub_test$ann_test <- sc_obj_sub_test$ann_level4_temp3_3
sc_obj_sub_test$ann_test[match(colnames(sc_obj_sub3_4),
                               colnames(sc_obj_sub_test))] <- paste0("a", sc_obj_sub3_4$seurat_clusters)
dot_motif1 <- DotPlot(sc_obj_sub_test,
                      assay = "chromvar",
                      col.min = -2,
                      col.max = 2,
                      features = sel_motif,
                      group.by = "ann_test",
                      cols = c("lightgrey", "brown")) +
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub3_4_test.tiff"),
     height = 10, width = 8, units = "in", res = 300, compression = "lzw")
print(dot_motif1)
dev.off()

##
tiff(file = glue("{out_path}UMAP_sub3_4_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_4, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_4",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp3_4_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_4",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 3_5 #####
sel_ct <- c("CD4_CD103_TRM_Th17")
sc_obj_sub3_5 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp3_4 %in% sel_ct])
DefaultAssay(sc_obj_sub3_5) <- "RNA"
sc_obj_sub3_5 <- NormalizeData(sc_obj_sub3_5)
#
sc_obj_sub3_5 <- scsub.renorm(sc_obj = sc_obj_sub3_5,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub3_5 <- dr.cl.wnn(sc_obj = sc_obj_sub3_5,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:15, 1:25),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_5_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub3_5, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3_5.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_5, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub3_5.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_5, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub3_5.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_5, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub3_5, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub3_5.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

# ## see in all-cell scale
# sc_obj_sub_test <- sc_obj
# sc_obj_sub_test$ann_test <- sc_obj_sub_test$ann_level4_temp3_4
# sc_obj_sub_test$ann_test[match(colnames(sc_obj_sub3_5),
#                                colnames(sc_obj_sub_test))] <- paste0("a", sc_obj_sub3_5$seurat_clusters)
# dot_motif1 <- DotPlot(sc_obj_sub_test,
#                       assay = "chromvar",
#                       col.min = -2,
#                       col.max = 2,
#                       features = sel_motif,
#                       group.by = "ann_test",
#                       cols = c("lightgrey", "brown")) +
#   scale_x_discrete(labels = names_motif) +
#   theme(axis.text.x = element_text(angle = 90))
# tiff(glue("{out_path}dot_raw_motif_sub3_5_test.tiff"),
#      height = 10, width = 8, units = "in", res = 300, compression = "lzw")
# print(dot_motif1)
# dev.off()

## temp anotation
sc_obj_sub3_5$ann_level4_temp3_5 <- "CD4_CD103_TRM_Th17"
sc_obj$ann_level4_temp3_5 <- sc_obj$ann_level4_temp3_4
sc_obj$ann_level4_temp3_5[match(colnames(sc_obj_sub3_5), colnames(sc_obj))] <- sc_obj_sub3_5$ann_level4_temp3_5

##
tiff(file = glue("{out_path}UMAP_sub3_5_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_5, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_5",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp3_5_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_5",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 3_6 #####
sel_ct <- c("CD4_CD103_TRM_Th17.1")
sc_obj_sub3_6 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp3_5 %in% sel_ct])
DefaultAssay(sc_obj_sub3_6) <- "RNA"
sc_obj_sub3_6 <- NormalizeData(sc_obj_sub3_6)
#
sc_obj_sub3_6 <- scsub.renorm(sc_obj = sc_obj_sub3_6,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub3_6 <- dr.cl.wnn(sc_obj = sc_obj_sub3_6,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_6_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub3_6, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3_6.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_6, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub3_6.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_6, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub3_6.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_6, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub3_6, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub3_6.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

# ## see in all-cell scale
# sc_obj_sub_test <- sc_obj
# sc_obj_sub_test$ann_test <- sc_obj_sub_test$ann_level4_temp3_5
# sc_obj_sub_test$ann_test[match(colnames(sc_obj_sub3_6),
#                                colnames(sc_obj_sub_test))] <- paste0("a", sc_obj_sub3_6$seurat_clusters)
# dot_motif1 <- DotPlot(sc_obj_sub_test,
#                       assay = "chromvar",
#                       col.min = -2,
#                       col.max = 2,
#                       features = sel_motif,
#                       group.by = "ann_test",
#                       cols = c("lightgrey", "brown")) +
#   scale_x_discrete(labels = names_motif) +
#   theme(axis.text.x = element_text(angle = 90))
# tiff(glue("{out_path}dot_raw_motif_sub3_6_test.tiff"),
#      height = 10, width = 8, units = "in", res = 300, compression = "lzw")
# print(dot_motif1)
# dev.off()

## temp anotation
sc_obj_sub3_6$ann_level4_temp3_6 <- "CD4_CD103_TRM_Th17.1"
sc_obj$ann_level4_temp3_6 <- sc_obj$ann_level4_temp3_5
sc_obj$ann_level4_temp3_6[match(colnames(sc_obj_sub3_6), colnames(sc_obj))] <- sc_obj_sub3_6$ann_level4_temp3_6

##
tiff(file = glue("{out_path}UMAP_sub3_6_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_6, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_6",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp3_6_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_6",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 3_7 #####
sel_ct <- c("CD4_CD103_TRM_Th1")
sc_obj_sub3_7 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level4_temp3_6 %in% sel_ct])
DefaultAssay(sc_obj_sub3_7) <- "RNA"
sc_obj_sub3_7 <- NormalizeData(sc_obj_sub3_7)
#
sc_obj_sub3_7 <- scsub.renorm(sc_obj = sc_obj_sub3_7,
                              do_ADT = F,
                              do_ATAC = T,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub3_7 <- dr.cl.wnn(sc_obj = sc_obj_sub3_7,
                           redc_list = list("harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:20, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           sprd = 0.6,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_7_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub3_7, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3_7.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_7, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub3_7.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_7, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD8", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "rna_TRDC", "rna_EGR1", "rna_HSPH1",
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
tiff(glue("{out_path}dot_rna_sub3_7.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_7, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub3_7, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub3_7.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

# ## see in all-cell scale
# sc_obj_sub_test <- sc_obj
# sc_obj_sub_test$ann_test <- sc_obj_sub_test$ann_level4_temp3_6
# sc_obj_sub_test$ann_test[match(colnames(sc_obj_sub3_7),
#                                colnames(sc_obj_sub_test))] <- paste0("a", sc_obj_sub3_7$seurat_clusters)
# dot_motif1 <- DotPlot(sc_obj_sub_test,
#                       assay = "chromvar",
#                       col.min = -2,
#                       col.max = 2,
#                       features = sel_motif,
#                       group.by = "ann_test",
#                       cols = c("lightgrey", "brown")) +
#   scale_x_discrete(labels = names_motif) +
#   theme(axis.text.x = element_text(angle = 90))
# tiff(glue("{out_path}dot_raw_motif_sub3_7_test.tiff"),
#      height = 10, width = 8, units = "in", res = 300, compression = "lzw")
# print(dot_motif1)
# dev.off()

## temp anotation
sc_obj_sub3_7$ann_level4_temp3_7 <- "CD4_CD103_TRM_Th1"
sc_obj$ann_level4_temp3_7 <- sc_obj$ann_level4_temp3_6
sc_obj$ann_level4_temp3_7[match(colnames(sc_obj_sub3_7), colnames(sc_obj))] <- sc_obj_sub3_7$ann_level4_temp3_7

##
tiff(file = glue("{out_path}UMAP_sub3_7_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_7, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_7",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp3_7_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_7",
        raster = F,
        label = T) + NoLegend()
dev.off()

