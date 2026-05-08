## module load r/4.5.0
library(Seurat)
library(bigreadr)
library(dplyr)
library(stringr)
library(glue)
library(harmony)

## set work dir
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
code_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/code/FUNCTION/"
source(glue("{code_path}process/PROCESS_FUN.R"))
##
adt_mk <- c("CD3","CD4", "CD8", "CD161", "CD127", "TCR-ab","TCR-Va7.2", "TCR-Vd2",
            "HLA-DR", "CD64", "CD123", "CD19", "CD16", "CD56", "CD1c", "CD11c", "CD11b")
feature_list <- list(
  "T_features_list" = c("CD2", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "TRDC", "IFNG"),
  "NK_features_list" = c("KLRF1","KLRD1","NCAM1","CD160", "FCGR3A", "KIT"),
  "MP_features_list" = c("CD14", "CD163", "CD68", "CSF1R", "TPSAB1") ,
  "DC_features_list" = c("CD1C", "LILRA4", "CD1E", "FCER1A"),
  "B_features_list" = c("CD79A", "MS4A1", "JSRP1", "MZB1", "CD38"),
  "Gra_features_list" = c("ENPP3", "ITGAM", "FUT4"),
  "PLT_features_list" = c("PF4", "PPBP", "GNG11")
)
##### load ADT and RNA data #####
data_path <- "03_output/03_clustering/all/"
out_path <- "03_output/03_clustering/all/WNN_ADT_RNA/"
scRNA_obj <- readRDS(glue("{data_path}scRNA_obj_test.rds"))
sc_obj <- scRNA_obj
scADT_obj <- readRDS(glue("{data_path}scADT_obj.rds"))
sc_obj[["ADT"]] <- scADT_obj[["ADT"]]
sc_obj@reductions$harmony_adt <- scADT_obj@reductions$harmony_adt

##### DRCL on all immune cells #####
seed_use <- 20250528
pc_clust_adt <- 20
pc_clust_gex <- 30
resx <- 1
sc_obj <- dr.cl.wnn(sc_obj = sc_obj,
                    redc_list = list("harmony_adt", "harmony_SCT"),
                    dim_list = list(1:pc_clust_adt, 1:pc_clust_gex),
                    k_nn = 20,
                    prune_SNN = 1/20,
                    n_iter = 500,
                    res = resx,
                    cl_method = 1,
                    run_umap = T,
                    n_neig = 30L,
                    n_epochs = 300,
                    neg_rate = 10L,
                    min_dist = 0.4,
                    sprd = 0.6,
                    seed_use = seed_use)

##
tiff(file = glue("{out_path}UMAP_raw.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(file = glue("{out_path}UMAP_condition.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "Condition", 
        label = T) + NoLegend()
dev.off()
##
tiff(file = glue("{out_path}UMAP_section.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "Section", 
        label = T) + NoLegend()
dev.off()

##### Cell type annotation #####
tiff(glue("{out_path}ft.tiff"), 
     height = 12, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", "adt_CD161", 
                         "adt_HLA-DR", "adt_CD64", "adt_CD123", "adt_CD19", 
                         "rna_CD3E", "rna_TRDC", "rna_KIT", "rna_TPSAB1"), 
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99",
            raster = F)
dev.off()
##
tiff(glue("{out_path}dot_rna.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text = element_text(angle = 90))
dev.off()

##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

##
sc_obj$ann_level1_temp1 <- "T"
sc_obj$ann_level1_temp1[sc_obj$wsnn_res.1 %in% c(10, 12, 25, 26, 28)] <- "B"
sc_obj$ann_level1_temp1[sc_obj$wsnn_res.1 %in% c(15, 22, 24, 27)] <- "ILC"
##
sc_obj$ann_level2_temp1 <- sc_obj$ann_level1_temp1
##
tiff(file = glue("{out_path}UMAP_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        group.by = "ann_level1_temp1",
        pt.size = 1E-5, 
        label = T) + 
  NoLegend()
dev.off()
##
tiff(file = glue("{out_path}UMAP_ann2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        group.by = "ann_level2_temp1",
        pt.size = 1E-5, 
        label = T) + 
  NoLegend()
dev.off()
#########################refine annotation#########################
##### round 1 #####
##### round 1.1 #####
sc_obj_sub1_1 <- subset(sc_obj, ann_level1_temp1 %in% c("ILC"))
DefaultAssay(sc_obj_sub1_1) <- "RNA"
sc_obj_sub1_1 <- NormalizeData(sc_obj_sub1_1)
#
sc_obj_sub1_1 <- scsub.renorm(sc_obj = sc_obj_sub1_1,
                            do_ADT = T,
                            do_ATAC = F,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")
sc_obj_sub1_1 <- dr.cl.wnn(sc_obj = sc_obj_sub1_1,
                         redc_list = list("harmony_adt", "harmony_SCT"),
                         dim_list = list(1:15, 1:10),
                         k_nn = 20,
                         prune_SNN = 1/20,
                         n_iter = 300,
                         res = 0.8,
                         run_umap = T,
                         n_neig = 30L,
                         n_epochs = 500,
                         neg_rate = 5L,
                         min_dist = 0.4,
                         seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub1_1.tiff"), 
     height = 12, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_1, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", "adt_CD161", 
                         "adt_HLA-DR", "adt_CD64", "adt_CD123", "adt_CD19", 
                         "rna_CD1C", "rna_CLEC9A", "rna_KIT", "rna_CD68"), 
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub1_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_1, 
        reduction = "wnn.umap", 
        # group.by = "ann_bulk",
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_1.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_1, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub1_1,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub1_1.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 3 7 8 11 13 19
sc_obj_sub1_1$ann_level1_temp2 <- sc_obj_sub1_1$ann_level1_temp1 %>% as.character()
sc_obj_sub1_1$ann_level1_temp2[sc_obj_sub1_1$wsnn_res.0.8 %in% c(9, 12, 15)] <- "T"
sc_obj_sub1_1$ann_level1_temp2[sc_obj_sub1_1$wsnn_res.0.8 %in% c(16)] <- "B"
sc_obj$ann_level1_temp2 <- sc_obj$ann_level2_temp1 %>% as.character()
sc_obj$ann_level1_temp2[match(colnames(sc_obj_sub1_1), colnames(sc_obj))] <- sc_obj_sub1_1$ann_level1_temp2

##
tiff(file = glue("{out_path}UMAP_sub1_1_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_1, 
        reduction = "wnn.umap", 
        group.by = "ann_level1_temp2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level1_temp2",
        label = T,
        raster = F) + NoLegend()
dev.off()
##### round 1_2 #####
sc_obj_sub1_2 <- subset(sc_obj, ann_level1_temp2 %in% c("B"))
DefaultAssay(sc_obj_sub1_2) <- "RNA"
sc_obj_sub1_2 <- NormalizeData(sc_obj_sub1_2)
#
sc_obj_sub1_2 <- scsub.renorm(sc_obj = sc_obj_sub1_2,
                            do_ADT = T,
                            do_ATAC = F,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")##
sc_obj_sub1_2 <- dr.cl.wnn(sc_obj = sc_obj_sub1_2,
                         redc_list = list("harmony_adt", "harmony_SCT"),
                         dim_list = list(1:30, 1:15),
                         k_nn = 20,
                         prune_SNN = 1/20,
                         n_iter = 300,
                         res = 0.6,
                         run_umap = T,
                         n_neig = 30L,
                         n_epochs = 500,
                         neg_rate = 20L,
                         min_dist = 0.4,
                         seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub1_2.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_2, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", "adt_CD161", 
                         "adt_HLA-DR", "adt_CD64", "adt_CD123", "adt_CD19", 
                         "rna_CD3E", "rna_MZB1", "rna_KIT", "rna_CD79A"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub1_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_2, 
        reduction = "wnn.umap", 
        # group.by = "ann_bulk",
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_2, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub1_2,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt1_2.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 3 7 8 11 13 19
sc_obj_sub1_2$ann_level1_temp3 <- sc_obj_sub1_2$ann_level1_temp2 %>% as.character()
sc_obj$ann_level1_temp3 <- sc_obj$ann_level1_temp2 %>% as.character()
sc_obj$ann_level1_temp3[match(colnames(sc_obj_sub1_2), colnames(sc_obj))] <- sc_obj_sub1_2$ann_level1_temp3

##
tiff(file = glue("{out_path}UMAP_sub1_2_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_2, 
        reduction = "wnn.umap", 
        group.by = "ann_level1_temp3",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level1_temp3",
        label = T,
        raster = F) + NoLegend()
dev.off()


##### round 1_3 #####
sc_obj_sub1_3 <- subset(sc_obj, ann_level1_temp3 %in% c("T"))
DefaultAssay(sc_obj_sub1_3) <- "RNA"
sc_obj_sub1_3 <- NormalizeData(sc_obj_sub1_3)
#
sc_obj_sub1_3 <- scsub.renorm(sc_obj = sc_obj_sub1_3,
                            do_ADT = T,
                            do_ATAC = F,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")##
sc_obj_sub1_3 <- dr.cl.wnn(sc_obj = sc_obj_sub1_3,
                         redc_list = list("harmony_adt", "harmony_SCT"),
                         dim_list = list(1:15, 1:25),
                         k_nn = 20,
                         prune_SNN = 1/20,
                         n_iter = 300,
                         res = 1,
                         run_umap = T,
                         n_neig = 30L,
                         n_epochs = 100,
                         neg_rate = 5L,
                         min_dist = 0.3,
                         seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub1_3.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_3, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", "adt_CD161", 
                         "adt_HLA-DR", "adt_CD64", "adt_CD123", "adt_CD19", 
                         "rna_CD3E", "rna_TRDC", "rna_KIT", "rna_TPSAB1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub1_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_3, 
        reduction = "wnn.umap", 
        # group.by = "ann_bulk",
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_3.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_3, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub1_3,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt1_3.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 3 7 8 11 13 19
sc_obj_sub1_3$ann_level1_temp4 <- sc_obj_sub1_3$ann_level1_temp3 %>% as.character()
sc_obj$ann_level1_temp4 <- sc_obj$ann_level1_temp3 %>% as.character()

##
tiff(file = glue("{out_path}UMAP_sub1_3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_3, 
        reduction = "wnn.umap", 
        group.by = "ann_level1_temp4",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level1_temp4",
        label = T,
        raster = F) + NoLegend()
dev.off()

#####round 2: determine ILC/Myeloid#####
sc_obj_sub2 <- subset(sc_obj, ann_level1_temp4 %in% c("ILC"))
DefaultAssay(sc_obj_sub2) <- "RNA"
sc_obj_sub2 <- NormalizeData(sc_obj_sub2)
#
sc_obj_sub2 <- scsub.renorm(sc_obj = sc_obj_sub2,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2 <- dr.cl.wnn(sc_obj = sc_obj_sub2,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:15, 1:10),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 0.8,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 500,
                           neg_rate = 5L,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub2.tiff"), 
     height = 12, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", "adt_CD161", 
                         "adt_HLA-DR", "adt_CD64", "adt_CD123", "adt_CD19", 
                         "rna_CD1C", "rna_CLEC9A", "rna_KIT", "rna_CD68"), 
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2, 
        reduction = "wnn.umap", 
        # group.by = "ann_bulk",
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub2.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 3 7 8 11 13 19
sc_obj_sub2$ann_level2_temp2 <- sc_obj_sub2$ann_level1_temp4 %>% as.character()
sc_obj_sub2$ann_level2_temp2[sc_obj_sub2$wsnn_res.0.8 %in% c(3)] <- "cDC"
sc_obj_sub2$ann_level2_temp2[sc_obj_sub2$wsnn_res.0.8 %in% c(9)] <- "pDC"
sc_obj_sub2$ann_level2_temp2[sc_obj_sub2$wsnn_res.0.8 %in% c(5, 11)] <- "Macrophage"
sc_obj_sub2$ann_level2_temp2[sc_obj_sub2$wsnn_res.0.8 %in% c(4, 8, 12)] <- "ILC"
sc_obj_sub2$ann_level2_temp2[sc_obj_sub2$wsnn_res.0.8 %in% c(0:2, 6, 7, 10)] <- "NK"
sc_obj$ann_level2_temp2 <- sc_obj$ann_level1_temp4 %>% as.character()
sc_obj$ann_level2_temp2[match(colnames(sc_obj_sub2), colnames(sc_obj))] <- sc_obj_sub2$ann_level2_temp2

##
tiff(file = glue("{out_path}UMAP_sub2_1_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp2",
        label = T,
        raster = F) + NoLegend()
dev.off()

## save ILC obj for plot as no further DRCL
sc_obj_sub2$ann_level2_final <- sc_obj_sub2$ann_level2_temp2
saveRDS(sc_obj_sub2, file = glue("{out_path}sc_obj_sub_ILC.rds"))

##### round 3: determine B/Plasma cells (just keep this level) #####
sc_obj_sub3 <- sc_obj_sub1_2
tiff(glue("{out_path}ft_sub3.tiff"), 
     height = 12, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", "adt_CD161", 
                         "adt_HLA-DR", "adt_CD64", "adt_CD123", "adt_CD19", 
                         "rna_CD1C", "rna_CLEC9A", "rna_KIT", "rna_CD68"), 
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3, 
        reduction = "wnn.umap", 
        # group.by = "ann_bulk",
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}dot_rna_sub3.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub3,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub3.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 3 7 8 11 13 19
sc_obj_sub3$ann_level2_temp3 <- "B"
sc_obj_sub3$ann_level2_temp3[sc_obj_sub3$wsnn_res.0.6 %in% c(3)] <- "Plasma"
sc_obj$ann_level2_temp3 <- sc_obj$ann_level2_temp2 %>% as.character()
sc_obj$ann_level2_temp3[match(colnames(sc_obj_sub3), colnames(sc_obj))] <- sc_obj_sub3$ann_level2_temp3

##
tiff(file = glue("{out_path}UMAP_sub2_2_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp3",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp3",
        label = T,
        raster = F) + NoLegend()
dev.off()


##### combine annotation #####
sc_obj@meta.data$ann_level2 <- sc_obj@meta.data$ann_level2_temp3 %>%
  factor(., levels = c("T", 
                       "NK", "ILC", 
                       "B", "Plasma", 
                       "Macrophage", "cDC", "pDC"))
##
sc_obj@meta.data$ann_level1 <- sc_obj@meta.data$ann_level2 %>% as.character()
sc_obj@meta.data$ann_level1[sc_obj@meta.data$ann_level2 %in% c("T")] <- "T"
sc_obj@meta.data$ann_level1[sc_obj@meta.data$ann_level2 %in% c("NK", "ILC")] <- "ILC"
sc_obj@meta.data$ann_level1[sc_obj@meta.data$ann_level2 %in% c("Macrophage", "cDC", "pDC")] <- "Myeloid"
sc_obj@meta.data$ann_level1 <- factor(sc_obj@meta.data$ann_level1, levels = c("T", "ILC", "B", "Plasma", "Myeloid"))
sc_meta <- cbind(sc_obj@meta.data,
                 Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta_ann_level2.rds"))

sc_obj_test <- sc_obj
DefaultAssay(sc_obj_test) <- "ADT"
sc_obj_test[["RNA"]] <- NULL
sc_obj_test[["SCT"]] <- NULL
saveRDS(sc_obj_test, file = glue("{out_path}scWNN_obj_drcl.rds"))
saveRDS(sc_obj, file = glue("{out_path}scWNN_obj.rds"))

##### round 4: dertemine T cells (ext round) #####
source("code/drcl/T/drcl_level2_T_adt.R")
source("code/drcl/T/drcl_level2_T_rna.R")
source("code/drcl/T/drcl_level2_T_wnn.R")
##
sc_obj <- readRDS(glue("{out_path}scWNN_obj_drcl.rds"))
sc_obj@meta.data <- readRDS(glue("{out_path}sc_meta_ann_level2.rds"))
sc_meta_t <- readRDS("03_output/03_clustering/T/WNN_ADT_RNA/sc_meta.rds")
sc_obj@meta.data$ann_level2_refine <- sc_obj@meta.data$ann_level2_temp3
sc_obj@meta.data$ann_level2_refine[match(rownames(sc_meta_t),
                                         colnames(sc_obj))] <- as.character(sc_meta_t$ann_level2_final)
sc_obj@meta.data$ann_level2_refine <- factor(sc_obj@meta.data$ann_level2_refine,
                                             levels = c("CD4T", "CD8T", "gdT", 
                                                        "B", "Plasma", 
                                                        "NK", "ILC", 
                                                        "Macrophage", "cDC", "pDC"))
##
tiff(file = glue("{out_path}UMAP_ann_level2_refine.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        group.by = "ann_level2_refine",
        reduction = "wnn.umap", 
        raster = F,
        label = T) + NoLegend()
dev.off()

##### save #####
sc_meta <- cbind(sc_obj@meta.data,
                 Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta_ann_level2_refine.rds"))
