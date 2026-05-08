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
adt_mk <- c("CD3","CD4", "CD8", "CD161", 
            "HLA-DR", "CD64", "CD123", "CD19")
rna_mk_list <- list(
  "T_features_list" = c("CD2", "CD3D", "CD3E", "CD3G"),
  "CD4_features_list" = c("CD4","FOXP3","CTLA4"),
  "CD8_features_list" = c("CD8A","CD8B", "LTB", "S100B"),
  "gdT_features_list" = c("TRGV9","TRDV2", "TRGC1","TRDC")
)
##### load ADT and RNA data #####
data_path <- "03_output/03_clustering/T/"
out_path <- "03_output/03_clustering/T/WNN_ADT_RNA/"
# scRNA_obj <- readRDS(glue("{data_path}scRNA_obj_test.rds"))
sc_obj <- scRNA_obj
sc_obj$wnnUMAP_1 <- NULL
sc_obj$wnnUMAP_2 <- NULL
scADT_obj <- readRDS(glue("{data_path}scADT_obj.rds"))
sc_obj[["ADT"]] <- scADT_obj[["ADT"]]
sc_obj@reductions$harmony_adt <- scADT_obj@reductions$harmony_adt
sc_obj@reductions$umap_adt <- scADT_obj@reductions$umap

##### DRCL on T immune cells #####
seed_use <- 20250528
n_adt <- 15
n_rna <- 20
res_use <- 1
sc_obj <- dr.cl.wnn(sc_obj = sc_obj,
                    redc_list = list("harmony_adt", "harmony_SCT"),
                    dim_list = list(1:n_adt, 1:n_rna),
                    k_nn = 20,
                    prune_SNN = 1/20,
                    n_iter = 300,
                    res = res_use,
                    run_umap = T,
                    n_neig = 30L,
                    n_epochs = 300,
                    neg_rate = 30L,
                    min_dist = 0.35,
                    sprd = 0.7,
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
out_path <- "03_output/03_clustering/T/WNN_ADT_RNA/"
# sc_obj <- readRDS(glue("{out_path}scWNN_obj_test.rds"))
# sc_meta <- readRDS(glue("{out_path}sc_meta.rds"))
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)

##
tiff(glue("{out_path}ft.tiff"), 
     height = 12, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", "adt_CD161", 
                         "adt_HLA-DR", "adt_CD64", "adt_CD123", "rna_CD8A", 
                         "rna_CD8B", "rna_TRDC", "rna_KIT", "rna_TPSAB1"), 
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", group.by = "seurat_clusters",
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
# annotation
sc_obj$ann_level2_temp1 <- "Mixed"
sc_obj$ann_level2_temp1[sc_obj$wsnn_res.1 %in% c(0, 2, 6, 10, 15, 17, 18)] <- "CD8T"
sc_obj$ann_level2_temp1[sc_obj$wsnn_res.1 %in% c(1, 3, 5, 7, 8, 9, 12, 19, 20)] <- "CD4T"
sc_obj$ann_level2_temp1[sc_obj$wsnn_res.1 %in% c(4, 23, 24)] <- "gdT"

tiff(file = glue("{out_path}UMAP_temp1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp1",
        label = T) + NoLegend()
dev.off()

#########################refine annotation#########################
##### round 1: extract mixed cells #####
##### round 1.1: mixed cells in CD4+T #####
sc_obj_sub1_1 <- subset(sc_obj, ann_level2_temp1 == "CD4T")
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
#
sc_obj_sub1_1 <- dr.cl.wnn(sc_obj = sc_obj_sub1_1,
                        redc_list = list("harmony_adt", "harmony_SCT"),
                        dim_list = list(1:20, 1:30),
                        k_nn = 20,
                        prune_SNN = 1/20,
                        n_iter = 300,
                        res = 1,
                        run_umap = T,
                        n_neig = 30L,
                        n_epochs = 300,
                        neg_rate = 10L,
                        min_dist = 0.3,
                        seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_1, 
        reduction = "wnn.umap",
        raster = F,
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub1_1.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_1, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            raster = F,
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_1.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_1, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
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
# mix_cells
mix_cell_1_1 <- c(colnames(sc_obj_sub1_1)[sc_obj_sub1_1$wsnn_res.1 %in% c(6, 14)],
                  colnames(sc_obj_sub1_1)[sc_obj_sub1_1[["ADT"]]$data["CD4", ] <= 
                                            sc_obj_sub1_1[["ADT"]]$data["CD8", ]]) %>% unique()

##### round 1.2: mixed cells in CD8+T #####
sc_obj_sub1_2 <- subset(sc_obj, ann_level2_temp1 == "CD8T")
DefaultAssay(sc_obj_sub1_2) <- "RNA"
sc_obj_sub1_2 <- NormalizeData(sc_obj_sub1_2)
#
sc_obj_sub1_2 <- scsub.renorm(sc_obj = sc_obj_sub1_2,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
#
sc_obj_sub1_2 <- dr.cl.wnn(sc_obj = sc_obj_sub1_2,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:20, 1:30),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 1,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 200,
                           neg_rate = 5L,
                           min_dist = 0.2,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_2, 
        reduction = "wnn.umap",
        raster = F,
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub1_2.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_2, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            raster = F,
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_2, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub1_2,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub1_2.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# mix_cells
mix_cell_1_2 <- c(colnames(sc_obj_sub1_2)[sc_obj_sub1_2$wsnn_res.1 %in% c(9, 12, 15)],
               colnames(sc_obj_sub1_2)[sc_obj_sub1_2[["ADT"]]$data["CD4", ] >= 
                                         sc_obj_sub1_2[["ADT"]]$data["CD8", ]]) %>% unique()
##### round 1.3: mixed cells in gdT #####
sc_obj_sub1_3 <- subset(sc_obj, ann_level2_temp1 == "gdT")
DefaultAssay(sc_obj_sub1_3) <- "RNA"
sc_obj_sub1_3 <- NormalizeData(sc_obj_sub1_3)
#
sc_obj_sub1_3 <- scsub.renorm(sc_obj = sc_obj_sub1_3,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
#
sc_obj_sub1_3 <- dr.cl.wnn(sc_obj = sc_obj_sub1_3,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:15, 1:20),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 1.5,
                           run_umap = T,
                           n_neig = 20L,
                           n_epochs = 500,
                           neg_rate = 10L,
                           min_dist = 0.3,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_3, 
        reduction = "wnn.umap",
        raster = F,
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub1_3.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_3, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            raster = F,
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_3.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_3, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub1_3,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub1_3.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# mix_cells
mix_cell_1_3 <- colnames(sc_obj_sub1_3)[sc_obj_sub1_3$wsnn_res.1.5 %in% c(5, 7, 12, 18)]

##### round 2: allocate for mixed cells  #####
## (Divede into Confirmed abT and Potential gdT cells)
## combined mixed cells
mix_cell1 <- Reduce("union", list(mix_cell_1_1, mix_cell_1_2, mix_cell_1_3))
sc_obj$ann_level2_temp2 <- sc_obj$ann_level2_temp1
sc_obj$ann_level2_temp2[match(mix_cell1, colnames(sc_obj))] <- "Mixed"
##
tiff(glue("{out_path}dot_rna_temp2.tiff"), 
     height = 4, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", 
        group.by = "ann_level2_temp2",
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    clus_col = "ann_level2_temp2",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_temp2.tiff"),
     height = 4, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

##### round 2.1 #####
sc_obj_sub2_1 <- subset(sc_obj, ann_level2_temp2 == "Mixed")
sc_obj_sub2_1 <- scsub.renorm(sc_obj = sc_obj_sub2_1,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_1 <- dr.cl.wnn(sc_obj = sc_obj_sub2_1,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:15, 1:20),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           cl_method = 1,
                           n_iter = 300,
                           res = 1,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.2,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_1, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub2_1.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_1, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            ncol = 4, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_1.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_1, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_1,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub2_1.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
# annotation
sc_obj_sub2_1$ann_level2_temp2_1 <- "Mixed_abT"
gd_clus <- c(3, 4, 7, 8, 12, 17, 19)
sc_obj_sub2_1$ann_level2_temp2_1[sc_obj_sub2_1$wsnn_res.1 %in% gd_clus] <- "Mixed_gdT"
sc_obj$ann_level2_temp2_1 <- sc_obj$ann_level2_temp2
sc_obj$ann_level2_temp2_1[match(colnames(sc_obj_sub2_1), colnames(sc_obj))] <- sc_obj_sub2_1$ann_level2_temp2_1
##
tiff(file = glue("{out_path}UMAP_mix_sub2_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_1, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp2_1",
        label = T) + NoLegend()
dev.off()

##### round 2.2 #####
sc_obj_sub2_2 <- subset(sc_obj, ann_level2_temp2_1 == "Mixed_gdT")
sc_obj_sub2_2 <- scsub.renorm(sc_obj = sc_obj_sub2_2,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_2 <- dr.cl.wnn(sc_obj = sc_obj_sub2_2,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           cl_method = 1,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.2,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_2, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub2_2.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_2, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            ncol = 4, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_2, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_2,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub2_2.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
# annotation
sc_obj_sub2_2$ann_level2_temp2_2 <- "Mixed_gdT"
sc_obj_sub2_2$ann_level2_temp2_2[sc_obj_sub2_2$wsnn_res.2 %in% c(0, 1, 5, 7, 8, 9, 
                                                                 11, 13, 14, 16, 18, 
                                                                 20:24, 26, 27, 28)] <- "Mixed_abT"
sc_obj$ann_level2_temp2_2 <- sc_obj$ann_level2_temp2_1
sc_obj$ann_level2_temp2_2[match(colnames(sc_obj_sub2_2), colnames(sc_obj))] <- sc_obj_sub2_2$ann_level2_temp2_2
##
tiff(file = glue("{out_path}UMAP_mix_sub2_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_2, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp2_2",
        label = T) + NoLegend()
dev.off()

##### round 2.3 #####
sc_obj_sub2_3 <- subset(sc_obj, ann_level2_temp2_2 == "Mixed_abT")
sc_obj_sub2_3 <- scsub.renorm(sc_obj = sc_obj_sub2_3,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_3 <- dr.cl.wnn(sc_obj = sc_obj_sub2_3,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           cl_method = 1,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.4,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_3, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub2_3.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_3, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            ncol = 4, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_3.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_3, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_3,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub2_3.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
# annotation
sc_obj_sub2_3$ann_level2_temp2_3 <- "Mixed_abT"
sc_obj_sub2_3$ann_level2_temp2_3[sc_obj_sub2_3$wsnn_res.2 %in% c(3, 9, 31, 34)] <- "Mixed_gdT"
sc_obj$ann_level2_temp2_3 <- sc_obj$ann_level2_temp2_2
sc_obj$ann_level2_temp2_3[match(colnames(sc_obj_sub2_3), colnames(sc_obj))] <- sc_obj_sub2_3$ann_level2_temp2_3
##
tiff(file = glue("{out_path}UMAP_mix_sub2_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_3, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp2_3",
        label = T) + NoLegend()
dev.off()

##### round 2.4 #####
sc_obj_sub2_4 <- subset(sc_obj, ann_level2_temp2_3 == "Mixed_gdT")
sc_obj_sub2_4 <- scsub.renorm(sc_obj = sc_obj_sub2_4,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_4 <- dr.cl.wnn(sc_obj = sc_obj_sub2_4,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           cl_method = 1,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.3,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_4, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub2_4.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_4, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            ncol = 4, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_4.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_4, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_4,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub2_4.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
# annotation
sc_obj_sub2_4$ann_level2_temp2_4 <- "Mixed_gdT"
sc_obj_sub2_4$ann_level2_temp2_4[sc_obj_sub2_4$wsnn_res.2 %in% c(1, 12, 20)] <- "Mixed_abT"
sc_obj$ann_level2_temp2_4 <- sc_obj$ann_level2_temp2_3
sc_obj$ann_level2_temp2_4[match(colnames(sc_obj_sub2_4), colnames(sc_obj))] <- sc_obj_sub2_4$ann_level2_temp2_4
##
tiff(file = glue("{out_path}UMAP_mix_sub2_4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_4, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp2_4",
        label = T) + NoLegend()
dev.off()


############## refine annotation 3 ############## 
sc_obj$ann_level2_temp3 <- sc_obj$ann_level2_temp2_4
sc_obj$ann_level2_temp3[sc_obj$ann_level2_temp3 == "Mixed_gdT"] <- "gdT"
sc_obj$ann_level2_temp3[sc_obj$ann_level2_temp3 == "Mixed_abT" & 
                          sc_obj[["ADT"]]$data["CD4", ] >= 
                          sc_obj[["ADT"]]$data["CD8", ]] <- "CD4T"
sc_obj$ann_level2_temp3[sc_obj$ann_level2_temp3 == "Mixed_abT" & 
                          sc_obj[["ADT"]]$data["CD4", ] < 
                          sc_obj[["ADT"]]$data["CD8", ]] <- "CD8T"
tiff(file = glue("{out_path}UMAP_sub3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp3",
        label = T) + NoLegend()
dev.off()
##### round 3.1 #####
sc_obj_sub3_1 <- subset(sc_obj, ann_level2_temp3 == "gdT")
DefaultAssay(sc_obj_sub3_1) <- "RNA"
sc_obj_sub3_1 <- NormalizeData(sc_obj_sub3_1)
#
sc_obj_sub3_1 <- scsub.renorm(sc_obj = sc_obj_sub3_1,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
#
sc_obj_sub3_1 <- dr.cl.wnn(sc_obj = sc_obj_sub3_1,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.4,
                           seed_use = seed_use)

##
tiff(file = glue("{out_path}UMAP_sub3_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_1, 
        reduction = "wnn.umap", 
        raster = F,
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub3_1.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_1, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            raster = F,
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub3_1.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_1, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub3_1,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub3_1.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation
sc_obj_sub3_1$ann_level2_temp3_1 <- "gdT"
ab_clus <- c(2, 14)
sc_obj_sub3_1$ann_level2_temp3_1[sc_obj_sub3_1$wsnn_res.2 %in% ab_clus & 
                                   sc_obj_sub3_1[["ADT"]]$data["CD4", ] >= 
                                   sc_obj_sub3_1[["ADT"]]$data["CD8", ]] <- "CD4T"
sc_obj_sub3_1$ann_level2_temp3_1[sc_obj_sub3_1$wsnn_res.2 %in% ab_clus & 
                                   sc_obj_sub3_1[["ADT"]]$data["CD4", ] < 
                                   sc_obj_sub3_1[["ADT"]]$data["CD8", ]] <- "CD8T"
sc_obj$ann_level2_temp3_1 <- sc_obj$ann_level2_temp3
sc_obj$ann_level2_temp3_1[match(colnames(sc_obj_sub3_1), colnames(sc_obj))] <- sc_obj_sub3_1$ann_level2_temp3_1
#
tiff(file = glue("{out_path}UMAP_temp3_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_1, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp3_1",
        label = T) + NoLegend()
dev.off()

##### round 3.2 #####
sc_obj_sub3_2 <- subset(sc_obj, ann_level2_temp3_1 == "CD4T")
DefaultAssay(sc_obj_sub3_2) <- "RNA"
sc_obj_sub3_2 <- NormalizeData(sc_obj_sub3_2)
#
sc_obj_sub3_2 <- scsub.renorm(sc_obj = sc_obj_sub3_2,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
#
sc_obj_sub3_2 <- dr.cl.wnn(sc_obj = sc_obj_sub3_2,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.3,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_2, 
        reduction = "wnn.umap", 
        raster = F,
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub3_2.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_2, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            raster = F,
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub3_2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_2, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub3_2,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub3_2.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation
sc_obj_sub3_2$ann_level2_temp3_2 <- "CD4T"
sc_obj_sub3_2$ann_level2_temp3_2[sc_obj_sub3_2$wsnn_res.2 %in% c(29)] <- "gdT"
sc_obj$ann_level2_temp3_2 <- sc_obj$ann_level2_temp3_1
sc_obj$ann_level2_temp3_2[match(colnames(sc_obj_sub3_2), colnames(sc_obj))] <- sc_obj_sub3_2$ann_level2_temp3_2

tiff(file = glue("{out_path}UMAP_temp3_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp3_2",
        label = T) + NoLegend()
dev.off()

##### round 3.3 #####
sc_obj_sub3_3 <- subset(sc_obj, ann_level2_temp3_2 == "CD8T")
DefaultAssay(sc_obj_sub3_3) <- "RNA"
sc_obj_sub3_3 <- NormalizeData(sc_obj_sub3_3)
#
sc_obj_sub3_3 <- scsub.renorm(sc_obj = sc_obj_sub3_3,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
#
sc_obj_sub3_3 <- dr.cl.wnn(sc_obj = sc_obj_sub3_3,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.2,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub3_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_3, 
        reduction = "wnn.umap", 
        raster = F,
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub3_3.tiff"), 
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_3, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8", 
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"), 
            raster = F,
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub3_3.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_3, 
        assay = "RNA", 
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub3_3,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub3_3.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation
sc_obj_sub3_3$ann_level2_temp3_3 <- sc_obj_sub3_3$ann_level2_temp3_2
sc_obj$ann_level2_temp3_3 <- sc_obj$ann_level2_temp3_2

tiff(file = glue("{out_path}UMAP_temp3_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp3_3",
        label = T) + NoLegend()
dev.off()

##### round 3.4 recheck gdT #####
sc_obj_sub3_4 <- subset(sc_obj, ann_level2_temp3_3 == "gdT")
DefaultAssay(sc_obj_sub3_4) <- "RNA"
sc_obj_sub3_4 <- NormalizeData(sc_obj_sub3_4)
#
sc_obj_sub3_4 <- scsub.renorm(sc_obj = sc_obj_sub3_4,
                              do_ADT = T,
                              do_ATAC = F,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
#
sc_obj_sub3_4 <- dr.cl.wnn(sc_obj = sc_obj_sub3_4,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.4,
                           seed_use = seed_use)

##
tiff(file = glue("{out_path}UMAP_sub3_4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_4,
        reduction = "wnn.umap",
        raster = F,
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub3_4.tiff"),
     height = 8, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_4,
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_CD4", "adt_CD8",
                         "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"),
            raster = F,
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub3_4.tiff"),
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_4,
        assay = "RNA",
        features = unlist(rna_mk_list) %>% unique) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub3_4,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub3_4.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

##
# annotation
sc_obj_sub3_4$ann_level2_temp3_4 <- "gdT"
sc_obj_sub3_4$ann_level2_temp3_4[sc_obj_sub3_4$wsnn_res.2 %in% c(17) & 
                                   sc_obj_sub3_4[["ADT"]]$data["CD4", ] >= 
                                   sc_obj_sub3_4[["ADT"]]$data["CD8", ]] <- "CD4T"
sc_obj_sub3_4$ann_level2_temp3_4[sc_obj_sub3_4$wsnn_res.2 %in% c(17) & 
                                   sc_obj_sub3_4[["ADT"]]$data["CD4", ] < 
                                   sc_obj_sub3_4[["ADT"]]$data["CD8", ]] <- "CD8T"
sc_obj$ann_level2_temp3_4 <- sc_obj$ann_level2_temp3_3
sc_obj$ann_level2_temp3_4[match(colnames(sc_obj_sub3_4), colnames(sc_obj))] <- sc_obj_sub3_4$ann_level2_temp3_4

#
tiff(file = glue("{out_path}UMAP_temp3_4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_4, 
        reduction = "wnn.umap", 
        group.by = "ann_level2_temp3_4",
        label = T) + NoLegend()
dev.off()

# ##### round 3.5 (recheck gdT and CD4T, PASSED) #####
# nn = 3
# test_ct <- c("gdT", "CD4T", "CD8T")[nn]
# n_adt_test <- c(20, 20, 20)[nn]
# n_rna_test <- c(15, 15, 15)[nn]
# #
# sc_obj_sub3_5 <- subset(sc_obj, ann_level2_temp3_4 == test_ct)
# DefaultAssay(sc_obj_sub3_5) <- "RNA"
# sc_obj_sub3_5 <- NormalizeData(sc_obj_sub3_5)
# #
# sc_obj_sub3_5 <- scsub.renorm(sc_obj = sc_obj_sub3_5,
#                               do_ADT = T,
#                               do_ATAC = F,
#                               do_RNA = F,
#                               do_SCT = T,
#                               do_harmony = T,
#                               batch_col = "Batch")
# #
# sc_obj_sub3_5 <- dr.cl.wnn(sc_obj = sc_obj_sub3_5,
#                            redc_list = list("harmony_adt", "harmony_SCT"),
#                            dim_list = list(1:n_adt_test, 1:n_rna_test),
#                            k_nn = 20,
#                            prune_SNN = 1/20,
#                            n_iter = 300,
#                            res = 2,
#                            run_umap = T,
#                            n_neig = 30L,
#                            n_epochs = 300,
#                            neg_rate = 10L,
#                            min_dist = 0.4,
#                            seed_use = seed_use)
# 
# ##
# tiff(file = glue("{out_path}UMAP_sub3_5_{test_ct}.tiff"),
#      width = 6, height = 6, units = "in", res = 600, compression = "lzw")
# DimPlot(sc_obj_sub3_5,
#         reduction = "wnn.umap",
#         raster = F,
#         label = T) + NoLegend()
# dev.off()
# ##
# tiff(glue("{out_path}ft_sub3_5_{test_ct}.tiff"),
#      height = 8, width = 16, units = "in", res = 300, compression = "lzw")
# FeaturePlot(sc_obj_sub3_5,
#             reduction = "wnn.umap",
#             features = c("adt_CD3", "adt_CD4", "adt_CD8",
#                          "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_TRDC", "rna_TPSAB1"),
#             raster = F,
#             ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
# dev.off()
# ##
# tiff(glue("{out_path}dot_rna_sub3_5_{test_ct}.tiff"),
#      height = 6, width = 12, units = "in", res = 300, compression = "lzw")
# DotPlot(sc_obj_sub3_5,
#         assay = "RNA",
#         features = unlist(rna_mk_list) %>% unique) +
#   theme(axis.text.x = element_text(angle = 90))
# dev.off()
# ##
# dot_adt <- heat.adt(sc_obj = sc_obj_sub3_5,
#                     clus_col = "seurat_clusters",
#                     assay_use = "ADT",
#                     adt_use = adt_mk)
# tiff(glue("{out_path}dot_adt_sub3_5_{test_ct}.tiff"),
#      height = 6, width = 4, units = "in", res = 300, compression = "lzw")
# dot_adt
# dev.off()


##### save #####
sc_obj@meta.data$ann_level2_final <- sc_obj@meta.data$ann_level2_temp3_4
Idents(sc_obj) <- sc_obj@meta.data$ann_level2_final
tiff(file = glue("{out_path}UMAP_ann_bulk2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()

tiff(glue("{out_path}dot_rna_final2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", 
        group.by = "ann_level2_final",
        features = unlist(rna_mk_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    clus_col = "ann_level2_final",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_final2.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

sc_meta <- cbind(sc_obj@meta.data,
                 Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta.rds"))
saveRDS(sc_obj, file = glue("{out_path}scWNN_obj.rds"))
