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
adt_mk <- c("CD4", "CD8","CD45RA", "CD45RO", "CD62L", "CD49a", "CD69", "CX3CR1",
            "CD183", "CD195", "CD103", "KLRG1","CD38",
            "CD194", "CD196", "CD39", "CD161", "CD26",
            "CD279", "CD185", "CD25", "CD127", "NKp80", "CD196", "CD162")
feature_list <- list(
  "T_features_list" = c("CD3D", "CD4", "CD8A", "CD8B","TRDC", "S1PR1", "CX3CR1"),
  "NK_features_list" = c("KLRD1", "KLRC2", "KLRC3", "NCR1"),
  "MAIT_features_list" = c("SLC4A10","TRAV1-2", 
                           "TRAV12-2", "TRAV8","TRAV21",
                           "KLRB1", "IL26", "IL23R", "IL22", "KIT"),
  "gdT_features_list" = c("TRGV9","TRDV2", "TRGC1","TRDC", "TYROBP"),
  "ISG_features_list" = c("ISG15","IFI6","IFI44L","LY6E", "EGR1", "EGR2", "CREM"),
  "Prolif_features_list" = c("MKI67","TYMS","PCNA"),
  "Act_features_list" = c("JUN", "EGR1", "TOX", "FAS"),
  "Naive_features_list" = c("CCR7","SELL","LEF1","TCF7", "IL7R"),
  "Effect_features_list" = c("GZMK","GZMA","GNLY","GZMB","FGFBP2","NKG7", 
                             "TNF", "IL2", "IL2RA", "IFNG", 'KLRG1'),
  "Memery_features_list" = c("GPR183","S100A4"),
  "Exhaust_features_list" = c("HAVCR2", "KLRK1", "LAG3", "CAV1"),
  "HSP_features_list" = c("HSPA1A", "HSPA1B", "HSPH1", "NR4A2"),
  "Th1_features_list" = c("CCL5", "CCR5"),
  "Th17_features_list" = c("RORC", "CCR6", "IL17A", "IL21", "IL4", "IL5", "IL9", "IL13"),
  "Tfh_features_list" = c("BCL6", "CXCR5"),
  "Treg_features_list" = c("FOXP3", "IKZF2", "CTLA4", "TIGIT")
)
atac_mk <- c("CCR7", "CCL5", "TBX21", "RORC", "IL17A", "BCL6", "CXCR5", "FOXP3", "IKZF2")
names_motif <- c("TBX21", "RORC", "RORA.1", "RORA.2", "GATA3", 
                 "FOS::JUND", "FOS", "KLF14", "SP1")
sel_motif <- lapply(names_motif, function(x){
  getMatrixSet(
    JASPAR2020, 
    opts = list(name = x, 
                tax_group = "vertebrates", collection = "CORE")
  ) %>% names()
}) %>% unlist
seed_use <- 20250528

##### load ADT, RNA, and ATAC data #####
data_path <- "03_output/03_clustering/CD8T/"
out_path <- "03_output/03_clustering/CD8T/WNN_ADT_RNA/"
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
#
cutoff_q <- "q70"
scATAC_obj <- readRDS(glue("{data_path}scATAC_obj_test_{cutoff_q}.rds"))
sc_obj[["peaks"]] <- scATAC_obj[["peaks"]]
sc_obj@reductions$harmony_lsi <- scATAC_obj@reductions$harmony_lsi
sc_obj@reductions$umap_lsi <- scATAC_obj@reductions$umap_lsi
#
sc_chromvar <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_scchromvar_assay.rds") %>%
  CreateSeuratObject(assay = "chromvar",
                     min.cells = 0,
                     min.features = 0)
sc_obj[["chromvar"]] <- subset(sc_chromvar, cells = colnames(sc_obj))[["chromvar"]]

##
sc_obj <- dr.cl.wnn(sc_obj = sc_obj,
                    redc_list = list("harmony_adt", "harmony_SCT", "harmony_lsi"),
                    dim_list = list(1:30, 1:15, 1:20),
                    k_nn = 20,
                    prune_SNN = 1/20,
                    n_iter = 300,
                    res = 3,
                    cl_method = 1,
                    run_umap = T,
                    n_neig = 20L,
                    n_epochs = 500,
                    neg_rate = 20L,
                    sprd = 0.5,
                    min_dist = 0.4,
                    seed_use = seed_use)
# sc_obj <- FindClusters(sc_obj,
#                        n.iter = 300,
#                        algorithm = 1,
#                        graph.name = "wsnn",
#                        resolution = 2,
#                        verbose = T,
#                        random.seed = seed_use)
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
sc_meta_ref <- readRDS("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/03_clustering_test/CD8T/test153000/test302015/sc_meta.rds")
sc_obj$ann_ref <- sc_meta_ref[colnames(sc_obj) %>% 
                                gsub("_DOGMAseq\\-", "", .) %>%
                                gsub("Duerr_", "", .), 
                              "ann_level3_final"]
##
tiff(file = glue("{out_path}UMAP_test_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj, !is.na(ann_ref)), 
        reduction = "wnn.umap", 
        group.by = "ann_ref",
        cols = cols4all::c4a("rainbow", 14),
        raster = F, 
        label = F) +
  theme(title = element_blank())
dev.off()

##
tiff(glue("{out_path}ft.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj, 
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
#
tiff(glue("{out_path}dot_rna_test.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", 
        col.min = -2, 
        col.max = 2, 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_test.tiff"),
     height = 10, width = 6, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
##
dot_motif <- DotPlot(sc_obj, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_test.tiff"), 
     height = 10, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## annotation
sc_obj$ann_level3_temp1 <- "CD8_TRM"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.3 %in% c(13, 38)] <- "CD8_naiveT"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.3 %in% c(10, 29, 30, 32)] <- "CD8_TEM"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.3 %in% c(24)] <- "CD8_TCM"

sc_obj$ann_level4_temp1 <- sc_obj$ann_level3_temp1
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

############## refine annotation ############## 
##### round 1 #####
sc_obj_sub1 <- subset(sc_obj, cells = colnames(sc_obj)[sc_obj$ann_level3_temp1 == "CD8_TRM"])
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
                         redc_list = list("harmony_SCT", "harmony_lsi"),
                         dim_list = list(1:15, 1:30),
                         k_nn = 20,
                         prune_SNN = 1/20,
                         n_iter = 300,
                         res = 1,
                         run_umap = T,
                         n_neig = 20L,
                         n_epochs = 300,
                         neg_rate = 20L,
                         sprd = 0.6,
                         min_dist = 0.4,
                         seed_use = seed_use)

##
sc_obj_sub1$ann_ref <- sc_meta_ref[colnames(sc_obj_sub1) %>%
                                     gsub("_DOGMAseq\\-", "", .) %>%
                                     gsub("Duerr_", "", .),
                                   "ann_level4_final"]
##
tiff(file = glue("{out_path}UMAP_sub1_ann_ref.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj_sub1, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        cols = cols4all::c4a("rainbow", 14),
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()
##
sc_obj_sub1 <- FindClusters(sc_obj_sub1,
                            n.iter = 300,
                            graph.name = "wsnn",
                            algorithm = 1,
                            resolution = 1,
                            verbose = T,
                            random.seed = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1, 
        reduction = "wnn.umap",
        label = T)
dev.off()

##
tiff(glue("{out_path}ft_sub1.tiff"), 
     height = 16, width = 26, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1, 
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
tiff(glue("{out_path}dot_rna_sub1.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1, 
        assay = "RNA", 
        group.by = "seurat_clusters",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_motif <- DotPlot(sc_obj_sub1, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "seurat_clusters",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_sub1.tiff"), 
     height = 8, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

# annotation 3 7 8 11 13 19
sc_obj_sub1$ann_level4_temp1 <- "CD8_TRM"
sc_obj_sub1$ann_level4_temp1[sc_obj_sub1$wsnn_res.1 %in% c(6, 8, 13, 15)] <- "CD8_TRM_Tc17_IL26"
sc_obj_sub1$ann_level4_temp1[sc_obj_sub1$wsnn_res.1 %in% c(10)] <- "CD8_TRM_IKZF2"
sc_obj_sub1$ann_level4_temp1[sc_obj_sub1$wsnn_res.1 %in% c(7)] <- "CD8_TRM_SP"
sc_obj_sub1$ann_level4_temp1[sc_obj_sub1$wsnn_res.1 %in% c(1)] <- "CD8_TRM_HSP"
sc_obj_sub1$ann_level4_temp1[sc_obj_sub1$wsnn_res.1 %in% c(11, 12)] <- "CD8_TRM_EGR1"
sc_obj_sub1$ann_level4_temp1[sc_obj_sub1$wsnn_res.1 %in% c(2)] <- "CD8_TRM_IL2"
sc_obj_sub1$ann_level4_temp1[sc_obj_sub1$wsnn_res.1 %in% c(3, 9)] <- "CD8_TRM_AP1"
##
tiff(file = glue("{out_path}UMAP_sub1_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1",
        label = T) + NoLegend()
dev.off()

## 
sc_obj$ann_level4_temp1 <- sc_obj$ann_level3_temp1
sc_obj$ann_level4_temp1[sc_obj$ann_level3_temp1 == "CD8_TEM"] <- "CD8_TEM_Tc1"
sc_obj$ann_level4_temp1[match(colnames(sc_obj_sub1), colnames(sc_obj))] <- sc_obj_sub1$ann_level4_temp1
##
tiff(file = glue("{out_path}UMAP_temp2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp1",
        label = T,
        raster = F) + NoLegend()
dev.off()

sc_obj_sub1$ann_level4_final <- factor(sc_obj_sub1$ann_level4_temp1,
                                       levels = c("CD8_TRM_Tc17_IL26",
                                                  "CD8_TRM",
                                                  "CD8_TRM_SP",
                                                  "CD8_TRM_EGR1",
                                                  "CD8_TRM_HSP",
                                                  "CD8_TRM_IL2",
                                                  "CD8_TRM_AP1",
                                                  "CD8_TRM_IKZF2"))
saveRDS(sc_obj_sub1, file = glue("{out_path}scWNN_obj_sub1.rds"))

##### save #####
sc_obj@meta.data$ann_level4_final <- sc_obj@meta.data$ann_level4_temp1
Idents(sc_obj) <- sc_obj@meta.data$ann_level4_final
sc_obj$ann_level4_final <- factor(sc_obj$ann_level4_final,
                                  levels = c("CD8_naiveT", 
                                             "CD8_TCM",
                                             "CD8_TEM_Tc1", 
                                             "CD8_TRM_Tc17_IL26", 
                                             "CD8_TRM",
                                             "CD8_TRM_SP",
                                             "CD8_TRM_EGR1", 
                                             "CD8_TRM_HSP",
                                             "CD8_TRM_IL2", 
                                             "CD8_TRM_AP1", 
                                             "CD8_TRM_IKZF2"))
Idents(sc_obj) <- sc_obj@meta.data$ann_level4_final
tiff(file = glue("{out_path}UMAP_ann_level4.tiff"),
     width = 9, height = 6.5, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        cols = cols4all::c4a("rainbow", 12),
        label = F)
dev.off()
##
sc_obj@meta.data$ann_level3_final <- factor(sc_obj@meta.data$ann_level3_temp1,
                                            levels = c("CD8_naiveT", 
                                                       "CD8_TCM",
                                                       "CD8_TEM", 
                                                       "CD8_TRM"))
##
Idents(sc_obj) <- sc_obj@meta.data$ann_level3_final
tiff(file = glue("{out_path}UMAP_ann_level3.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        cols = cols4all::c4a("rainbow", 4),
        label = F)
dev.off()
##
tiff(glue("{out_path}dot_rna_level4_final.tiff"), 
     height = 6, width = 15, units = "in", res = 300, compression = "lzw")
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
##
dot_motif <- DotPlot(sc_obj, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "ann_level4_final",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_motif_level4_final.tiff"), 
     height = 6, width = 6, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

##
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

