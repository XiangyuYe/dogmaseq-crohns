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
  "B_features_list" = c('MS4A1', 'CD79A'),
  "Bnaive_features_list" = c("IGHD", "FCER2", "TCL1A", "IL4R", "CD83", "CD69", "CREM"),
  "Bmem_features_list" = c("CD27", "AIM2", "TNFRSF13B", "TXNIP", "MT1X", "MT2A", "SLC30A1", "DUSP4"),
  "GC_features_list" = c("AICDA", "GCSAM",
                         'LRMP', 'SUGCT', 'MME', 'MKI67', 'BCL6', "RGS13", 
                         "CCL3", "CCL4"),
  "ISG_features_list" = c("ISG15","IFI6","IFI44L","LY6E"),
  "HSP_features_list" = c("HSPA1A", "HSPA1B", "IFNG", "NR4A2"),
  "Act_features_list" = c("JUN", "EGR1", "TOX", "FAS"),
  "ABC_features_list" = c("FCRLA", "FCRL3", "FCRL4", "FCRL5", "ITGAX", "TBX21"),
  "Breg_features_list" = c("IL10", "CD1D", "CD5", "TGFB1"),
  "Plasma_features_list" = c("CD38", "SDC1", "MZB1", 
                             "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHM", "IGHV3-7", "IGKC"),
  "MBC_features_list" = c("CD44", "GPR183")
)
adt_mk <- c("CD19", "CD161", "CD38", "CD24", "CD71", "CD44", "CD62L",
            "IgD", "IgE", "IgM", "CD11c", "CD1c", "CD86", "CD83",
            "HLA-DR", "HLA-A-B-C", "CD27", "CD25", "CD123", "CD11b", "CD57", "CD20")
seed_use <- 20250528

##### load ADT, RNA, and ATAC data #####
data_path <- "03_output/03_clustering/B/"
out_path <- "03_output/03_clustering/B/WNN_ADT_RNA/"
sc_obj <- readRDS(glue("{data_path}scRNA_obj_test.rds"))
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)
sc_obj$wnnUMAP_1 <- NULL
sc_obj$wnnUMAP_2 <- NULL
#
scADT_obj <- readRDS(glue("{data_path}scADT_obj.rds"))
sc_obj[["ADT"]] <- scADT_obj[["ADT"]]
sc_obj@reductions$harmony_adt <- scADT_obj@reductions$harmony_adt
sc_obj@reductions$umap_adt <- scADT_obj@reductions$umap

##### DRCL on Plasma #####
sc_obj <- dr.cl.wnn(sc_obj = sc_obj,
                    redc_list = list(
                      "harmony_adt",
                      "harmony_SCT"),
                    dim_list = list(
                      1:25,
                      1:15),
                    k_nn = 20,
                    prune_SNN = 1/20,
                    n_iter = 300,
                    res = 1,
                    run_umap = T,
                    n_neig = 20L,
                    n_epochs = 500,
                    neg_rate = 30L,
                    min_dist = 0.35,
                    seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_test.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F, 
        label = T) +
  theme(title = element_blank())
dev.off()
## annotation
sc_obj$ann_level4_temp1 <- sc_obj$ann_level3_temp1 <- sc_obj$ann_level2_refine
##
tiff(file = glue("{out_path}UMAP_level3_temp1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        group.by = "ann_level3_temp1",
        reduction = "wnn.umap", 
        label = T,
        raster = F) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_level4_temp1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        group.by = "ann_level4_temp1",
        reduction = "wnn.umap", 
        label = T,
        raster = F) + NoLegend()
dev.off()

############## refine annotation ############## 
##### round 1 #####
sc_obj_sub1 <- subset(sc_obj, ann_level2_refine %in% c("B"))
DefaultAssay(sc_obj_sub1) <- "RNA"
sc_obj_sub1 <- NormalizeData(sc_obj_sub1)
#
sc_obj_sub1 <- scsub.renorm(sc_obj = sc_obj_sub1,
                            do_ADT = T,
                            do_ATAC = T,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")
sc_obj_sub1 <- dr.cl.wnn(sc_obj = sc_obj_sub1,
                         redc_list = list("harmony_adt", "harmony_SCT"),
                         dim_list = list(1:30, 1:15),
                         k_nn = 20,
                         prune_SNN = 1/20,
                         n_iter = 300,
                         res = 1,
                         run_umap = T,
                         n_neig = 30L,
                         n_epochs = 500,
                         neg_rate = 20L,
                         min_dist = 0.4,
                         seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub1.tiff"), 
     height = 12, width = 16, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_IgM", "adt_CD25", "adt_IgD", 
                         "adt_CD27", "adt_CD19", "rna_IGHG1", "rna_IGHD", 
                         "rna_IGHA1", "rna_TXNIP", "rna_MZB1", "rna_GPR183"), 
            ncol = 4, min.cutoff = "q1", max.cutoff = "q99")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1, 
        reduction = "wnn.umap", 
        # group.by = "ann_bulk",
        label = T) + NoLegend()
dev.off()
##
sc_meta_ref <- readRDS("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/03_clustering_test/B/sc_meta.rds")
sc_obj_sub1$ann_ref <- sc_meta_ref[colnames(sc_obj_sub1) %>%
                                gsub("_DOGMAseq\\-", "", .) %>%
                                gsub("Duerr_", "", .),
                              "ann_level3_final"]
##
tiff(file = glue("{out_path}UMAP_sub1_ann_ref.tiff"),
     width = 9, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub1, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        cols = cols4all::c4a("rainbow", 13),
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub1,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_sub1.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 3 7 8 11 13 19
sc_obj_sub1$ann_level3_temp2 <- "B_memory"
sc_obj_sub1$ann_level3_temp2[sc_obj_sub1$wsnn_res.1 %in% c(9)] <- "GC"
sc_obj_sub1$ann_level3_temp2[sc_obj_sub1$wsnn_res.1 %in% c(10)] <- "ABC"
sc_obj_sub1$ann_level3_temp2[sc_obj_sub1$wsnn_res.1 %in% c(1:3, 6, 7, 14, 15)] <- "B_naive"
#
sc_obj_sub1$ann_level4_temp2 <- sc_obj_sub1$ann_level3_temp2
sc_obj_sub1$ann_level4_temp2[sc_obj_sub1$wsnn_res.1 %in% c(3, 13)] <- "B_naive_Activated"
sc_obj_sub1$ann_level4_temp2[sc_obj_sub1$wsnn_res.1 %in% c(8)] <- "B_memory_CD4"
sc_obj_sub1$ann_level4_temp2[sc_obj_sub1$wsnn_res.1 %in% c(4)] <- "B_memory_TOX"
#
sc_obj$ann_level3_temp2 <- sc_obj$ann_level3_temp1 %>% as.character()
sc_obj$ann_level3_temp2[match(colnames(sc_obj_sub1), colnames(sc_obj))] <- sc_obj_sub1$ann_level3_temp2
sc_obj$ann_level4_temp2 <- sc_obj$ann_level4_temp1 %>% as.character()
sc_obj$ann_level4_temp2[match(colnames(sc_obj_sub1), colnames(sc_obj))] <- sc_obj_sub1$ann_level4_temp2
##
tiff(file = glue("{out_path}UMAP_sub1_ann_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2",
        label = T,
        raster = F) + NoLegend()
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub1_ann_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp2",
        label = T,
        raster = F) + NoLegend()
dev.off()

##### round 2 #####
sc_obj_sub2 <- subset(sc_obj, ann_level3_temp2 %in% c("Plasma"))
DefaultAssay(sc_obj_sub2) <- "RNA"
sc_obj_sub2 <- NormalizeData(sc_obj_sub2)
#
sc_obj_sub2 <- scsub.renorm(sc_obj = sc_obj_sub2,
                            do_ADT = T,
                            do_ATAC = T,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")##
sc_obj_sub2 <- dr.cl.wnn(sc_obj = sc_obj_sub2,
                         redc_list = list("harmony_adt", "harmony_SCT"),
                         dim_list = list(1:15, 1:8),
                         k_nn = 20,
                         prune_SNN = 1/15,
                         n_iter = 300,
                         res = 1,
                         run_umap = T,
                         n_neig = 30L,
                         n_epochs = 500,
                         neg_rate = 10L,
                         min_dist = 0.4,
                         seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub2.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2, 
            reduction = "wnn.umap",
            features = c("adt_CD3", "adt_IgM", "adt_CD25", "adt_IgD", 
                         "adt_CD27", "adt_CD19", "rna_IGHG1", "rna_IGHD", 
                         "rna_IGHA1", "rna_CD27", "rna_EGR1", "rna_IGKC"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
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
sc_meta_ref <- readRDS("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/03_clustering_test/B/sc_meta.rds")
sc_obj_sub2$ann_ref <- sc_meta_ref[colnames(sc_obj_sub2) %>%
                                     gsub("_DOGMAseq\\-", "", .) %>%
                                     gsub("Duerr_", "", .),
                                   "ann_level3_final"]
##
tiff(file = glue("{out_path}UMAP_sub2_ann_ref.tiff"),
     width = 9, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub2, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        cols = cols4all::c4a("rainbow", 5),
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()

##
tiff(glue("{out_path}dot_rna_sub2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt2.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation
sc_obj_sub2$ann_level3_temp3 <- "Plasma_IgA"
sc_obj_sub2$ann_level3_temp3[sc_obj_sub2$wsnn_res.1 %in% c(3)] <- "Plasma_IgG"
sc_obj$ann_level3_temp3 <- sc_obj$ann_level3_temp2 %>% as.character()
sc_obj$ann_level3_temp3[match(colnames(sc_obj_sub2), colnames(sc_obj))] <- sc_obj_sub2$ann_level3_temp3
#
sc_obj_sub2$ann_level4_temp3 <- sc_obj_sub2$ann_level3_temp3
sc_obj$ann_level4_temp3 <- sc_obj$ann_level4_temp2 %>% as.character()
sc_obj$ann_level4_temp3[match(colnames(sc_obj_sub2), colnames(sc_obj))] <- sc_obj_sub2$ann_level4_temp3

##
tiff(file = glue("{out_path}UMAP_sub2_ann_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp3_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3",
        label = T,
        raster = F) + NoLegend()
dev.off()

##### combine annotation #####
ct_level <- list("B" = list("ann_level4_final" = c("B_naive", "B_naive_Activated",
                                                   "B_memory", "B_memory_TOX", "B_memory_CD4", 
                                                   "ABC",
                                                   "GC",
                                                   "Plasma_IgA", "Plasma_IgG"),
                            "ann_level3_final" = c("B_naive", "B_memory",
                                                   "ABC",
                                                   "GC",
                                                   "Plasma_IgA", "Plasma_IgG"),
                            "ann_level2_final" = c("B", "Plasma")))
##### level 3 #####
sc_obj$ann_level3_final <- sc_obj$ann_level3_temp3
sc_obj$ann_level3_final <- factor(sc_obj$ann_level3_final,
                                  levels = ct_level$B$ann_level3_final)
Idents(sc_obj) <- sc_obj$ann_level3_final
tiff(file = glue("{out_path}UMAP_ann_level3.tiff"),
     width = 9, height = 6.5, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        cols = cols4all::c4a("rainbow", nlevels(sc_obj$ann_level3_final)),
        label = F)
dev.off()
##
tiff(glue("{out_path}dot_rna_level3_final.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", 
        group.by = "ann_level3_final",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    clus_col = "ann_level3_final",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_level3_final.tiff"),
     height = 6, width = 6, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

##### level 4 #####
sc_obj$ann_level4_final <- sc_obj$ann_level4_temp3
sc_obj$ann_level4_final <- factor(sc_obj$ann_level4_final,
                                  levels = ct_level$B$ann_level4_final)
Idents(sc_obj) <- sc_obj$ann_level4_final
tiff(file = glue("{out_path}UMAP_ann_level4.tiff"),
     width = 9, height = 6.5, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        cols = cols4all::c4a("rainbow", nlevels(sc_obj$ann_level4_final)),
        label = F)
dev.off()
##
tiff(glue("{out_path}dot_rna_level4_final.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", 
        group.by = "ann_level4_final",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    clus_col = "ann_level4_final",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_level4_final.tiff"),
     height = 6, width = 6, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

##### output #####
sc_meta <- cbind(sc_obj@meta.data,
                 Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta.rds"))
fwrite2(sc_meta, file = glue("{out_path}sc_meta.txt"), row.names = T)
#
saveRDS(sc_obj@reductions, file = glue("{out_path}reducWNN_test.rds"))
saveRDS(sc_obj@commands, file = glue("{out_path}cmdWNN_test.rds"))

sc_obj_test <- sc_obj
DefaultAssay(sc_obj_test) <- "ADT"
sc_obj_test[["peaks"]] <- NULL
sc_obj_test[["RNA"]] <- NULL
sc_obj_test[["SCT"]] <- NULL
saveRDS(sc_obj_test, file = glue("{out_path}scWNN_obj_drcl.rds"))
saveRDS(sc_obj, file = glue("{out_path}scWNN_obj.rds"))

