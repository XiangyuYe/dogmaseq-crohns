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
            "CD183", "CD195", "CD103", "CD161", "KLRG1",
            "CD194", "CD196", "HLA-DR",
            "CD279", "CD278", "CD185", "CD25", "CD127")
feature_list <- list(
  "T_features_list" = c("CD3D", "CD4", "CD8A", "CD8B","TRDC", "S1PR1", "CX3CR1"),
  "Prolif_features_list" = c("MKI67","TYMS","PCNA"),
  "Naive_features_list" = c("CCR7","SELL","LEF1","TCF7", "KLF2"),
  "Memery_features_list" = c("GPR183","S100A4"),
  "Effect_features_list" = c("GZMK","GZMA", "GZMB", "GNLY","NKG7", "RUNX3"),
  "Exhaust_features_list" = c("HAVCR2", "KLRK1", "LAG3", "CAV1"),
  "Th1_features_list" = c("CCL5", "CCR5", "IFNG"),
  "Th17_features_list" = c("RORC", "CCR6", "IL17A", "IL4", "IL5", "IL9", "IL13"),
  "Tfh_features_list" = c("IL21", "BCL6", "CXCR5", "ICOS", "CXCL13"),
  "Treg_features_list" = c("FOXP3", "IKZF2")
)
atac_mk <- c("CCR7", "CCL5", "TBX21", "RORC", "IL17A", "BCL6", "CXCR5", "FOXP3", "IKZF2")
sel_motif <- c("MA0690.1", "MA1151.1", "MA0071.1", "MA0072.1", "MA0037.3", "MA0462", "MA0850.1")
names_motif <- c("TBX21", "RORC", "RORA.1", "RORA.2", "GATA3", "BCL6", "FOXP3")
seed_use <- 20250528

##### load ADT, RNA, and ATAC data #####
data_path <- "03_output/03_clustering/CD4T/"
out_path <- "03_output/03_clustering/CD4T/WNN_ADT_RNA/"
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
# scATAC_obj <- readRDS(glue("{data_path}scATAC_obj_test_{cutoff_q}.rds"))
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
                    redc_list = list(
                      "harmony_adt",
                      "harmony_SCT", 
                      "harmony_lsi"),
                    dim_list = list(
                      1:30,
                      1:20, 
                      1:15),
                    k_nn = 20,
                    prune_SNN = 1/20,
                    n_iter = 300,
                    res = 2,
                    cl_method = 1,
                    run_umap = T,
                    n_neig = 30L,
                    n_epochs = 500,
                    neg_rate = 20L,
                    sprd = 0.5,
                    min_dist = 0.4,
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
##
tiff(glue("{out_path}ft.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
ft_mk <- FeaturePlot(sc_obj, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD127", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "adt_CD278", "adt_CD279", "rna_S100A4",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
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
dot_adt <- heat.adt(sc_obj = sc_obj,
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_test.tiff"),
     height = 10, width = 5, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
##
atac_cd4 <- c("BCL6", "CXCR5", "TBX21", "CCL5", "IFNG", 
              "IL17A", "IL17F", "RORA", "RORC", "BATF", "FOXP3")
cv_plt <- peak.set.plt(sc_obj = sc_obj,
                       assay_use = "peaks",
                       mk_list = atac_cd4,
                       extend_kb = 3000)
ggsave(file = glue("{out_path}dot_plt_test.png"),
       cv_plt,
       height = 18, width = 15,units = "in", dpi = 300, limitsize = F)

## annotation 2 5 7 12 18 20 23 26 30
sc_obj$ann_level3_temp1 <- "Mixed"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(5, 21, 29, 30)] <- "CD4_naiveT"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(7, 9, 10, 23)] <- "CD4_Treg"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(0, 14, 15, 17, 25)] <- "CD4_TCM"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(27)] <- "CD4_TEM"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(1, 2, 3, 13, 16, 22, 26)] <- "CD4_CD103_TRM"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(6, 11, 12, 31)] <- "CD4_TRM"

## annotation 2
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
tiff(file = glue("{out_path}UMAP_level4_temp1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        group.by = "ann_level4_temp1",
        reduction = "wnn.umap", 
        label = T,
        raster = F) + NoLegend()
dev.off()

plt4x <- DotPlot(sc_obj, 
                 assay = "RNA", 
                 col.min = -2, 
                 col.max = 2, 
                 group.by = "ann_level4_temp1",
                 features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_rna_level3_temp1.tiff"), 
     height = 5, width = 15, units = "in", res = 300, compression = "lzw")
print(plt4x)
dev.off()

tiff(glue("{out_path}dot_raw_adt_level3_temp1.tiff"),
     height = 5, width = 6, units = "in", res = 300, compression = "lzw")
heat.adt(sc_obj,
         clus_col = "ann_level4_temp1",
         assay_use = "ADT",
         adt_use = adt_mk) %>% print
dev.off()

############## refine annotation ############## 
##### round 1: Mixed cells #####
sc_obj_sub1_1 <- subset(sc_obj, ann_level3_temp1 == "Mixed")
DefaultAssay(sc_obj_sub1_1) <- "RNA"
sc_obj_sub1_1 <- NormalizeData(sc_obj_sub1_1)
#
sc_obj_sub1_1 <- scsub.renorm(sc_obj = sc_obj_sub1_1,
                              do_ADT = T,
                              do_ATAC = T,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub1_1 <- dr.cl.wnn(sc_obj = sc_obj_sub1_1,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:30, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 3,
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
        label = T) + NoLegend()
dev.off()
##
ft_plt <- FeaturePlot(sc_obj_sub1_1, 
                      reduction = "wnn.umap",
                      features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                                   "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                                   "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                                   "chromvar_MA0690.1", "chromvar_MA1151.1"), 
                      ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
ft_plt[[14]] <- ft_plt[[14]] + labs(title = "chromvar_TBX21")
ft_plt[[15]] <- ft_plt[[15]] + labs(title = "chromvar_RORC")
tiff(glue("{out_path}ft_sub1_1.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
ft_plt
dev.off()
##
tiff(glue("{out_path}dot_rna_sub1_1.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub1_1, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub1_1,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt1_1.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 
# 37      36 4.703655
# 33      32 4.722765
# 28      27 4.762657
# 11      10 4.840630
# 36      35 4.886868
# 8        7 5.188337
# 20      19 5.342027
# 7        6 5.411723
# 25      24 5.610231
# 14      13 5.731013
# 9        8 5.959906
# 4        3 6.061103
# 3        2 6.570429
# 31      30 7.107240
# 12      11 7.533789
sc_obj_sub1_1$ann_level3_temp1_2 <- sc_obj_sub1_1$ann_level3_temp1
sc_obj_sub1_1$ann_level3_temp1_2[sc_obj_sub1_1$wsnn_res.3 %in% c(2, 3, 6, 7, 8, 10, 11, 13, 19, 
                                                                 24, 27, 30, 32, 35, 36)] <- "CD4_CD103_TRM"
sc_obj_sub1_1$ann_level3_temp1_2[sc_obj_sub1_1$wsnn_res.3 %in% c(0, 9, 15, 17, 20, 25, 26, 34)] <- "CD4_TRM"
sc_obj_sub1_1$ann_level3_temp1_2[sc_obj_sub1_1$wsnn_res.3 %in% c(14, 22, 28, 31)] <- "CD4_TCM"
sc_obj$ann_level3_temp1_2 <- sc_obj$ann_level3_temp1
sc_obj$ann_level3_temp1_2[match(colnames(sc_obj_sub1_1), colnames(sc_obj))] <- sc_obj_sub1_1$ann_level3_temp1_2
##
tiff(file = glue("{out_path}UMAP_sub1_1_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_1, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp1_2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_2_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp1_2",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 1.2: Mixed cells2 #####
sc_obj_sub1_2 <- subset(sc_obj, ann_level3_temp1_2 == "Mixed")
DefaultAssay(sc_obj_sub1_2) <- "RNA"
sc_obj_sub1_2 <- NormalizeData(sc_obj_sub1_2)
#
sc_obj_sub1_2 <- scsub.renorm(sc_obj = sc_obj_sub1_2,
                              do_ADT = T,
                              do_ATAC = F,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub1_2 <- dr.cl.wnn(sc_obj = sc_obj_sub1_2,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:15, 1:20),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 3,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.3,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_2, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub1_2.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_2, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                         "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
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
# 4        3 1.403102
# 11      10 1.499011
# 10       9 1.699234
# 20      19 1.712214
# 1        0 1.736601
# 21      20 1.749400
# 25      24 1.879192

# annotation 
sc_obj_sub1_2$ann_level3_temp1_3 <- "Mixed"
sc_obj_sub1_2$ann_level3_temp1_3[sc_obj_sub1_2$wsnn_res.3 %in% c(1, 2, 15, 23)] <- "CD4_CD103_TRM"
sc_obj_sub1_2$ann_level3_temp1_3[sc_obj_sub1_2$wsnn_res.3 %in% c(0, 3, 9, 10, 20, 24)] <- "CD4_TRM"
sc_obj_sub1_2$ann_level3_temp1_3[sc_obj_sub1_2$wsnn_res.3 %in% c(19)] <- "CD4_TEM"
sc_obj$ann_level3_temp1_3 <- sc_obj$ann_level3_temp1_2
sc_obj$ann_level3_temp1_3[match(colnames(sc_obj_sub1_2), colnames(sc_obj))] <- sc_obj_sub1_2$ann_level3_temp1_3
##
tiff(file = glue("{out_path}UMAP_sub1_2_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_2, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp1_3",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_3_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp1_3",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 1.3: Mixed cells4 #####
sc_obj_sub1_3 <- subset(sc_obj, ann_level3_temp1_3 == "Mixed")
DefaultAssay(sc_obj_sub1_3) <- "RNA"
sc_obj_sub1_3 <- NormalizeData(sc_obj_sub1_3)
#
sc_obj_sub1_3 <- scsub.renorm(sc_obj = sc_obj_sub1_3,
                              do_ADT = T,
                              do_ATAC = F,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub1_3 <- dr.cl.wnn(sc_obj = sc_obj_sub1_3,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:10, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/15,
                           n_iter = 300,
                           res = 3,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.3,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub1_3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_3, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub1_3.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub1_3, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                         "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
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

# annotation 
sc_obj_sub1_3$ann_level3_temp1_4 <- "CD4_TRM"
sc_obj_sub1_3$ann_level3_temp1_4[sc_obj_sub1_3$wsnn_res.3 %in% c(8, 12, 20)] <- "CD4_CD103_TRM"
sc_obj$ann_level3_temp1_4 <- sc_obj$ann_level3_temp1_3
sc_obj$ann_level3_temp1_4[match(colnames(sc_obj_sub1_3), colnames(sc_obj))] <- sc_obj_sub1_3$ann_level3_temp1_4
##
tiff(file = glue("{out_path}UMAP_sub1_3_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub1_3, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp1_4",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp1_4_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp1_4",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2: refine level3 #####
##### round 2.1: CD4_naiveT #####
sc_obj_sub2_1 <- subset(sc_obj, ann_level3_temp1_4 == "CD4_naiveT")
DefaultAssay(sc_obj_sub2_1) <- "RNA"
sc_obj_sub2_1 <- NormalizeData(sc_obj_sub2_1)
##
sc_obj_sub2_1 <- scsub.renorm(sc_obj = sc_obj_sub2_1,
                            do_ADT = T,
                            do_ATAC = T,
                            q_atac = cutoff_q,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")
## 4
sc_obj_sub2_1 <- dr.cl.wnn(sc_obj = sc_obj_sub2_1,
                           redc_list = list("harmony_adt", "harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:30, 1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 1,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 100,
                           neg_rate = 5L,
                           min_dist = 0.2,
                           seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub2_1.tiff"), 
     height = 12, width = 20, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_1, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD62L", "adt_CD69", 
                         "adt_CD103", "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_FOXP3", "rna_CCL5", "rna_IKZF2", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub2_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_1, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_1.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_1, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
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
sc_obj_sub2_1$ann_level3_temp2_1 <- sc_obj_sub2_1$ann_level3_temp1_4
sc_obj$ann_level3_temp2_1 <- sc_obj$ann_level3_temp1_4
sc_obj$ann_level3_temp2_1[match(colnames(sc_obj_sub2_1), colnames(sc_obj))] <- sc_obj_sub2_1$ann_level3_temp2_1

##
tiff(file = glue("{out_path}UMAP_sub2_1_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_1, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_1",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_1_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_1",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2.2: CD4_Treg #####
sc_obj_sub2_2 <- subset(sc_obj, ann_level3_temp2_1 == "CD4_Treg")
DefaultAssay(sc_obj_sub2_2) <- "RNA"
sc_obj_sub2_2 <- NormalizeData(sc_obj_sub2_2)
##
sc_obj_sub2_2 <- scsub.renorm(sc_obj = sc_obj_sub2_2,
                            do_ADT = T,
                            do_ATAC = T,
                            q_atac = cutoff_q,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")
##
sc_obj_sub2_2 <- dr.cl.wnn(sc_obj = sc_obj_sub2_2,
                           redc_list = list("harmony_adt", "harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:30, 1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 1,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 100,
                           neg_rate = 5L,
                           min_dist = 0.2,
                           seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub2_2.tiff"), 
     height = 12, width = 20, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_2, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD62L", "adt_CD69", 
                         "adt_CD103", "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_FOXP3", "rna_CCL5", "rna_IKZF2", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub2_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_2, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_2, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
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
sc_obj_sub2_2$ann_level3_temp2_2 <- sc_obj_sub2_2$ann_level3_temp2_1
sc_obj$ann_level3_temp2_2 <- sc_obj$ann_level3_temp2_1
sc_obj$ann_level3_temp2_2[match(colnames(sc_obj_sub2_2), colnames(sc_obj))] <- sc_obj_sub2_2$ann_level3_temp2_2

##
tiff(file = glue("{out_path}UMAP_sub2_2_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_2, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_2",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_2_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_2",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2.3: CD4_TCM #####
sc_obj_sub2_3 <- subset(sc_obj, ann_level3_temp2_2 == "CD4_TCM")
DefaultAssay(sc_obj_sub2_3) <- "RNA"
sc_obj_sub2_3 <- NormalizeData(sc_obj_sub2_3)
#
sc_obj_sub2_3 <- scsub.renorm(sc_obj = sc_obj_sub2_3,
                            do_ADT = T,
                            do_ATAC = T,
                            q_atac = cutoff_q,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")
sc_obj_sub2_3 <- dr.cl.wnn(sc_obj = sc_obj_sub2_3,
                           redc_list = list("harmony_adt", "harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:30, 1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 2,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 100,
                           neg_rate = 5L,
                           min_dist = 0.2,
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
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_3, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                         "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_3.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_3, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_3,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt2_3.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 
sc_obj_sub2_3$ann_level3_temp2_3 <- sc_obj_sub2_3$ann_level3_temp2_2
sc_obj$ann_level3_temp2_3 <- sc_obj$ann_level3_temp2_2
sc_obj$ann_level3_temp2_3[match(colnames(sc_obj_sub2_3), colnames(sc_obj))] <- sc_obj_sub2_3$ann_level3_temp2_3
##
tiff(file = glue("{out_path}UMAP_sub2_3_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_3, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_3",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_3_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_3",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2.4: CD4_TRM #####
sc_obj_sub2_4 <- subset(sc_obj, ann_level3_temp2_3 == "CD4_TRM")
DefaultAssay(sc_obj_sub2_4) <- "RNA"
sc_obj_sub2_4 <- NormalizeData(sc_obj_sub2_4)
#
sc_obj_sub2_4 <- scsub.renorm(sc_obj = sc_obj_sub2_4,
                              do_ADT = T,
                              do_ATAC = T,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_4 <- dr.cl.wnn(sc_obj = sc_obj_sub2_4,
                           redc_list = list("harmony_adt", "harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:30, 1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 3,
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
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_4, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                         "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_4.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_4, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_4,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt2_4.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 
sc_obj_sub2_4$ann_level3_temp2_4 <- sc_obj_sub2_4$ann_level3_temp2_3
sc_obj_sub2_4$ann_level3_temp2_4[sc_obj_sub2_4$wsnn_res.3 %in% c(6, 17, 19, 20, 22)] <- "CD4_TCM"
sc_obj_sub2_4$ann_level3_temp2_4[sc_obj_sub2_4$wsnn_res.3 %in% c(7)] <- "CD4_CD103_TRM"
sc_obj_sub2_4$ann_level3_temp2_4[sc_obj_sub2_4[["ADT"]]$data["CD103",] > 8] <- "CD4_CD103_TRM"
sc_obj$ann_level3_temp2_4 <- sc_obj$ann_level3_temp2_3
sc_obj$ann_level3_temp2_4[match(colnames(sc_obj_sub2_4), colnames(sc_obj))] <- sc_obj_sub2_4$ann_level3_temp2_4
##
tiff(file = glue("{out_path}UMAP_sub2_4_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_4, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_4",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_4_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_4",
        raster = F,
        label = T) + NoLegend()
dev.off()


##### round 2.5: CD4_CD103_TRM #####
sc_obj_sub2_5 <- subset(sc_obj, ann_level3_temp2_4 == "CD4_CD103_TRM")
DefaultAssay(sc_obj_sub2_5) <- "RNA"
sc_obj_sub2_5 <- NormalizeData(sc_obj_sub2_5)
#
sc_obj_sub2_5 <- scsub.renorm(sc_obj = sc_obj_sub2_5,
                              do_ADT = T,
                              do_ATAC = T,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_5 <- dr.cl.wnn(sc_obj = sc_obj_sub2_5,
                           redc_list = list("harmony_adt", "harmony_SCT", "harmony_lsi"),
                           dim_list = list(1:30, 1:20, 1:15),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 3,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.3,
                           seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_5.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_5, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub2_5.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_5, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                         "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_5.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_5, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_5,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt2_5.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 
sc_obj_sub2_5$ann_level3_temp2_5 <- sc_obj_sub2_5$ann_level3_temp2_4
sc_obj_sub2_5$ann_level3_temp2_5[sc_obj_sub2_5$wsnn_res.3 %in% c(1, 3, 12, 20 ,21, 29)] <- "CD4_TRM"
sc_obj_sub2_5$ann_level3_temp2_5[sc_obj_sub2_5[["ADT"]]$data["CD103",] < 1] <- "CD4_TRM"
sc_obj$ann_level3_temp2_5 <- sc_obj$ann_level3_temp2_4
sc_obj$ann_level3_temp2_5[match(colnames(sc_obj_sub2_5), colnames(sc_obj))] <- sc_obj_sub2_5$ann_level3_temp2_5
##
tiff(file = glue("{out_path}UMAP_sub2_5_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_5, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_5",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_5_level3.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_5",
        raster = F,
        label = T)
dev.off()

##### round 2.6: CD4_TRM2 #####
sc_obj_sub2_6 <- subset(sc_obj, ann_level3_temp2_5 == "CD4_TRM")
#
sc_obj_sub2_6 <- scsub.renorm(sc_obj = sc_obj_sub2_6,
                              do_ADT = T,
                              do_ATAC = F,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = F,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_6 <- FindNeighbors(sc_obj_sub2_6, 
                               reduction = "harmony_adt",
                               n.trees = 300,
                               k.param = 20, 
                               dims = 1:20) %>%
  FindClusters(., 
               n.iter = 300,
               resolution = 1,
               random.seed = seed_use) %>% 
  RunUMAP(., 
          dims = 1:20, 
          reduction = 'harmony_adt',
          n.neighbors = 30L,
          umap.method = "uwot",
          n.epochs = 300,
          negative.sample.rate = 10L,
          min.dist = 0.4,
          seed.use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_6.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_6, 
        reduction = "umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub2_6.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_6, 
            reduction = "umap",
            features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                         "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_6.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_6, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_6,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt2_6.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 
sc_obj_sub2_6$ann_level3_temp2_6 <- sc_obj_sub2_6$ann_level3_temp2_5
sc_obj_sub2_6$ann_level3_temp2_6[sc_obj_sub2_6$ADT_snn_res.1 %in% c(4)] <- "CD4_CD103_TRM"
sc_obj_sub2_6$ann_level3_temp2_6[sc_obj_sub2_6$ADT_snn_res.1 %in% c(0, 10) & 
                                   sc_obj_sub2_6[["ADT"]]$data["CD103",] > 5] <- "CD4_CD103_TRM"
sc_obj_sub2_6$ann_level3_temp2_6[sc_obj_sub2_6$ADT_snn_res.1 %in% c(12)] <- "CD4_TCM"
sc_obj$ann_level3_temp2_6 <- sc_obj$ann_level3_temp2_5
sc_obj$ann_level3_temp2_6[match(colnames(sc_obj_sub2_6), colnames(sc_obj))] <- sc_obj_sub2_6$ann_level3_temp2_6
##
tiff(file = glue("{out_path}UMAP_sub2_6_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_6, 
        reduction = "umap", 
        group.by = "ann_level3_temp2_6",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_6_level3.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_6",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### round 2.7: CD4_CD103_TRM2 #####
sc_obj_sub2_7 <- subset(sc_obj, ann_level3_temp2_6 == "CD4_CD103_TRM")
#
sc_obj_sub2_7 <- scsub.renorm(sc_obj = sc_obj_sub2_7,
                              do_ADT = T,
                              do_ATAC = F,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = F,
                              do_harmony = T,
                              batch_col = "Batch")
sc_obj_sub2_7 <- FindNeighbors(sc_obj_sub2_7, 
                               reduction = "harmony_adt",
                               n.trees = 300,
                               k.param = 20, 
                               dims = 1:20) %>%
  FindClusters(., 
               n.iter = 300,
               resolution = 1,
               random.seed = seed_use) %>% 
  RunUMAP(., 
          dims = 1:20, 
          reduction = 'harmony_adt',
          n.neighbors = 30L,
          umap.method = "uwot",
          n.epochs = 300,
          negative.sample.rate = 10L,
          min.dist = 0.4,
          seed.use = seed_use)
##
tiff(file = glue("{out_path}UMAP_sub2_7.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_7, 
        reduction = "umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub2_7.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_7, 
            reduction = "umap",
            features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                         "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_7.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_7, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_7,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt2_7.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 
sc_obj_sub2_7$ann_level3_temp2_7 <- sc_obj_sub2_7$ann_level3_temp2_6
sc_obj$ann_level3_temp2_7 <- sc_obj$ann_level3_temp2_6
sc_obj$ann_level3_temp2_7[match(colnames(sc_obj_sub2_7), colnames(sc_obj))] <- sc_obj_sub2_7$ann_level3_temp2_7
##
tiff(file = glue("{out_path}UMAP_sub2_7_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_7, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_7",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_7_level3.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_7",
        raster = F,
        label = T)
dev.off()

##### round 2.8: CD4_TCM2 #####
sc_obj_sub2_8 <- subset(sc_obj, ann_level3_temp2_7 == "CD4_TCM")
DefaultAssay(sc_obj_sub2_8) <- "RNA"
sc_obj_sub2_8 <- NormalizeData(sc_obj_sub2_8)
##
sc_obj_sub2_8 <- scsub.renorm(sc_obj = sc_obj_sub2_8,
                              do_ADT = T,
                              do_ATAC = T,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
##
sc_obj_sub2_8 <- dr.cl.wnn(sc_obj = sc_obj_sub2_8,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:15, 1:30),
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
tiff(file = glue("{out_path}UMAP_sub2_8.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_8, 
        reduction = "umap", 
        label = T) + NoLegend()
dev.off()
##
tiff(glue("{out_path}ft_sub2_8.tiff"), 
     height = 12, width = 22, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub2_8, 
            reduction = "umap",
            features = c("adt_CD4", "adt_CD103", "adt_CD62L", "adt_CD69", 
                         "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_CXCR5", "rna_CCL5", "rna_FOXP3", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(glue("{out_path}dot_rna_sub2_8.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub2_8, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj_sub2_8,
                    clus_col = "seurat_clusters",
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt2_8.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

# annotation 
sc_obj_sub2_8$ann_level3_temp2_8 <- sc_obj_sub2_8$ann_level3_temp2_7
sc_obj_sub2_8$ann_level3_temp2_8[sc_obj_sub2_8$wsnn_res.2 %in% c(20)] <- "CD4_Treg"
sc_obj$ann_level3_temp2_8 <- sc_obj$ann_level3_temp2_7
sc_obj$ann_level3_temp2_8[match(colnames(sc_obj_sub2_8), colnames(sc_obj))] <- sc_obj_sub2_8$ann_level3_temp2_8

##
tiff(file = glue("{out_path}UMAP_sub2_8_level3_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub2_8, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_8",
        label = T) + NoLegend()
dev.off()

##
tiff(file = glue("{out_path}UMAP_temp2_8_level3.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level3_temp2_8",
        raster = F,
        label = T)
dev.off()

# ##### recheck TRM (PASSED) #####
# test_ct <- c("CD4_TRM", "CD4_CD103_TRM")[2]
# sc_obj_sub1_test <- subset(sc_obj, ann_level3_temp2_7 == test_ct)
# DefaultAssay(sc_obj_sub1_test) <- "RNA"
# sc_obj_sub1_test <- NormalizeData(sc_obj_sub1_test)
# ##
# sc_obj_sub1_test <- scsub.renorm(sc_obj = sc_obj_sub1_test,
#                                  do_ADT = T,
#                                  do_ATAC = T,
#                                  q_atac = cutoff_q,
#                                  do_RNA = F,
#                                  do_SCT = T,
#                                  do_harmony = T,
#                                  batch_col = "Batch")
# ##
# sc_obj_sub1_test <- dr.cl.wnn(sc_obj = sc_obj_sub1_test,
#                               redc_list = list("harmony_adt", "harmony_SCT"),
#                               dim_list = list(1:15, 1:30),
#                               k_nn = 20,
#                               prune_SNN = 1/20,
#                               n_iter = 300,
#                               res = 2,
#                               run_umap = T,
#                               n_neig = 30L,
#                               n_epochs = 300,
#                               neg_rate = 10L,
#                               min_dist = 0.3,
#                               seed_use = seed_use)
# ##
# tiff(glue("{out_path}ft_sub1_test_{test_ct}.tiff"),
#      height = 12, width = 20, units = "in", res = 300, compression = "lzw")
# FeaturePlot(sc_obj_sub1_test,
#             reduction = "wnn.umap",
#             features = c("adt_CD4", "adt_CD62L", "adt_CD69",
#                          "adt_CD103", "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
#                          "rna_CCR7", "rna_FOXP3", "rna_CCL5", "rna_IKZF2", "rna_CCR6",
#                          "chromvar_MA0690.1", "chromvar_MA1151.1"),
#             ncol = 5, min.cutoff = "q10", max.cutoff = "q90")
# dev.off()
# ##
# tiff(file = glue("{out_path}UMAP_sub1_test_{test_ct}.tiff"),
#      width = 6, height = 6, units = "in", res = 600, compression = "lzw")
# DimPlot(sc_obj_sub1_test,
#         reduction = "wnn.umap",
#         label = T) + NoLegend()
# dev.off()
# ##
# tiff(glue("{out_path}dot_rna_sub1_test_{test_ct}.tiff"),
#      height = 6, width = 12, units = "in", res = 300, compression = "lzw")
# DotPlot(sc_obj_sub1_test,
#         assay = "RNA",
#         features = unlist(feature_list) %>% unique) +
#   theme(axis.text.x = element_text(angle = 90))
# dev.off()
# ##
# dot_adt <- heat.adt(sc_obj = sc_obj_sub1_test,
#                     clus_col = "seurat_clusters",
#                     assay_use = "ADT",
#                     adt_use = adt_mk)
# tiff(glue("{out_path}dot_adt_sub1_test_{test_ct}.tiff"),
#      height = 6, width = 4, units = "in", res = 300, compression = "lzw")
# dot_adt
# dev.off()

##### round 3: define level4 #####
sc_obj$ann_level3_final <- sc_obj$ann_level3_temp2_8

##### round 3.1: Tregs #####
sc_obj_sub3_1 <- subset(sc_obj, ann_level3_final == "CD4_Treg")
DefaultAssay(sc_obj_sub3_1) <- "RNA"
sc_obj_sub3_1 <- NormalizeData(sc_obj_sub3_1)
##
sc_obj_sub3_1 <- scsub.renorm(sc_obj = sc_obj_sub3_1,
                            do_ADT = F,
                            do_ATAC = T,
                            q_atac = cutoff_q,
                            do_RNA = F,
                            do_SCT = T,
                            do_harmony = T,
                            batch_col = "Batch")
##
sc_obj_sub3_1 <- dr.cl.wnn(sc_obj = sc_obj_sub3_1,
                         redc_list = list("harmony_SCT", "harmony_lsi"),
                         dim_list = list(1:15, 1:20),
                         k_nn = 20,
                         prune_SNN = 1/20,
                         n_iter = 300,
                         res = 0.5,
                         run_umap = T,
                         n_neig = 30L,
                         n_epochs = 200,
                         neg_rate = 10L,
                         min_dist = 0.3,
                         seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub3_1.tiff"), 
     height = 12, width = 20, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_1, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD62L", "adt_CD69", 
                         "adt_CD103", "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_FOXP3", "rna_CCL5", "rna_IKZF2", "rna_CCR6", 
                         "chromvar_MA0690.1", "chromvar_MA1151.1"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
#
tiff(file = glue("{out_path}UMAP_sub3_1.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_1, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()

##
tiff(glue("{out_path}dot_rna_sub3_1.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_1, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
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
sc_obj_sub3_1$ann_level4_temp3_1 <- sc_obj_sub3_1$ann_level3_final
sc_obj_sub3_1$ann_level4_temp3_1[sc_obj_sub3_1$wsnn_res.0.5 %in% c(1, 5)] <- "CD4_IKZF2low_Treg"
sc_obj_sub3_1$ann_level4_temp3_1[sc_obj_sub3_1$wsnn_res.0.5 %in% c(0, 4)] <- "CD4_Treg"
sc_obj_sub3_1$ann_level4_temp3_1[sc_obj_sub3_1$wsnn_res.0.5 %in% c(2, 3, 6)] <- "CD4_Treg_naive"

sc_obj$ann_level4_temp3_1 <- sc_obj$ann_level3_final
sc_obj$ann_level4_temp3_1[match(colnames(sc_obj_sub3_1), colnames(sc_obj))] <- sc_obj_sub3_1$ann_level4_temp3_1

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

##### round 3.2: TCM #####
sc_obj_sub3_2 <- subset(sc_obj, ann_level3_final == "CD4_TCM")
DefaultAssay(sc_obj_sub3_2) <- "RNA"
sc_obj_sub3_2 <- NormalizeData(sc_obj_sub3_2)
##
sc_obj_sub3_2 <- scsub.renorm(sc_obj = sc_obj_sub3_2,
                              do_ADT = T,
                              do_ATAC = T,
                              q_atac = cutoff_q,
                              do_RNA = F,
                              do_SCT = T,
                              do_harmony = T,
                              batch_col = "Batch")
##
sc_obj_sub3_2 <- dr.cl.wnn(sc_obj = sc_obj_sub3_2,
                           redc_list = list("harmony_adt", "harmony_SCT"),
                           dim_list = list(1:15, 1:20),
                           k_nn = 20,
                           prune_SNN = 1/20,
                           n_iter = 300,
                           res = 0.5,
                           run_umap = T,
                           n_neig = 30L,
                           n_epochs = 300,
                           neg_rate = 10L,
                           min_dist = 0.3,
                           seed_use = seed_use)
##
tiff(glue("{out_path}ft_sub3_2.tiff"), 
     height = 12, width = 20, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj_sub3_2, 
            reduction = "wnn.umap",
            features = c("adt_CD4", "adt_CD62L", "adt_CD69", 
                         "adt_CD103", "adt_CD185", "adt_CD183", "adt_CD196", "adt_CD279",
                         "rna_CCR7", "rna_FOXP3", "rna_CCL5", "rna_IL21", "rna_CCR6", 
                         "rna_CXCR5", "rna_CXCL13"), 
            ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
dev.off()
##
tiff(file = glue("{out_path}UMAP_sub3_2.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_2, 
        reduction = "wnn.umap", 
        label = T) + NoLegend()
dev.off()

##
tiff(glue("{out_path}dot_rna_sub3_2.tiff"), 
     height = 6, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj_sub3_2, 
        assay = "RNA", 
        features = unlist(feature_list) %>% unique) + 
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
sc_obj_sub3_2$ann_level4_temp3_2 <- "CD4_TCM"
sc_obj_sub3_2$ann_level4_temp3_2[sc_obj_sub3_2$wsnn_res.0.5 %in% c(1, 4, 5)] <- "CD4_CD69low_TCM"
sc_obj_sub3_2$ann_level4_temp3_2[sc_obj_sub3_2$wsnn_res.0.5 %in% c(3)] <- "CD4_TCM_Tfh_like"

sc_obj$ann_level4_temp3_2 <- sc_obj$ann_level4_temp3_1
sc_obj$ann_level4_temp3_2[match(colnames(sc_obj_sub3_2), colnames(sc_obj))] <- sc_obj_sub3_2$ann_level4_temp3_2

##
tiff(file = glue("{out_path}UMAP_sub3_2_level4_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj_sub3_2, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_2",
        label = T) + NoLegend()
dev.off()

# Idents(sc_obj_sub3_22) <- sc_obj_sub3_2$ann_level4_temp3_2
# tiff(file = glue("{out_path}UMAP_sub3_22.tiff"),
#      width = 6, height = 6, units = "in", res = 600, compression = "lzw")
# DimPlot(sc_obj_sub3_22, 
#         reduction = "wnn.umap", 
#         label = T) + NoLegend()
# dev.off()
##
tiff(file = glue("{out_path}UMAP_temp3_2_level4.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        group.by = "ann_level4_temp3_2",
        raster = F,
        label = T) + NoLegend()
dev.off()

##### output annotation level 3.5 #####
sc_obj@meta.data$ann_level3.5_final <- sc_obj@meta.data$ann_level4_temp3_2
sc_obj$ann_level3.5_final <- factor(sc_obj$ann_level3.5_final,
                                    levels = c("CD4_naiveT", 
                                               "CD4_CD69low_TCM", "CD4_TCM", "CD4_TCM_Tfh_like", 
                                               "CD4_TEM",
                                               "CD4_TRM", "CD4_CD103_TRM",
                                               "CD4_Treg_naive", "CD4_IKZF2low_Treg", "CD4_Treg"))
##
Idents(sc_obj) <- sc_obj@meta.data$ann_level3.5_final
tiff(file = glue("{out_path}UMAP_ann_level3.5.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        cols = cols4all::c4a("rainbow", 10),
        label = F)
dev.off()
##
sc_obj$ann_level3_final <- factor(sc_obj$ann_level3_final,
                                  levels = c("CD4_naiveT", 
                                             "CD4_TCM", "CD4_TEM",
                                             "CD4_TRM", "CD4_CD103_TRM", "CD4_Treg"))
##
Idents(sc_obj) <- sc_obj@meta.data$ann_level3_final
tiff(file = glue("{out_path}UMAP_ann_level3.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        cols = cols4all::c4a("rainbow", 6),
        label = F)
dev.off()
##
sc_meta <- cbind(sc_obj@meta.data,
                 Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta_ann_level3.5_final.rds"))

##### Subclustering on TRM #####
# sc_obj <- readRDS(glue("{out_path}scWNN_obj.rds"))
sc_obj@meta.data <- readRDS(glue("{out_path}sc_meta_ann_level3.5_final.rds"))
sc_meta_trm <- readRDS("03_output/03_clustering/CD4_TRM/WNN_RNA_ATAC/sc_meta.rds")
sc_meta_cd103trm <- readRDS("03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/sc_meta.rds")
sc_obj$ann_level4_final <- sc_obj$ann_level3.5_final %>% as.character()
sc_obj$ann_level4_final[sc_obj$ann_level4_final == "CD4_TEM"] <- "CD4_TEM_Th1"
sc_obj$ann_level4_final[match(rownames(sc_meta_trm), colnames(sc_obj))] <- sc_meta_trm$ann_level4_final %>% as.character()
sc_obj$ann_level4_final[match(rownames(sc_meta_cd103trm), colnames(sc_obj))] <- sc_meta_cd103trm$ann_level4_final %>% as.character()
sc_obj$ann_level4_final <- factor(sc_obj$ann_level4_final,
                                   levels = c("CD4_naiveT",
                                              "CD4_CD69low_TCM", "CD4_TCM", "CD4_TCM_Tfh_like",
                                              "CD4_TEM_Th1",
                                              "CD4_TRM", "CD4_TRM_Th1", "CD4_CD103_TRM_Th1",
                                              "CD4_TRM_Th17", "CD4_CD103_TRM_Th17",
                                              "CD4_TRM_Th17.1", "CD4_CD103_TRM_Th17.1",
                                              "CD4_Treg_naive", "CD4_IKZF2low_Treg", "CD4_Treg"))

##
Idents(sc_obj) <- sc_obj@meta.data$ann_level4_final
tiff(file = glue("{out_path}UMAP_ann_level4_final.tiff"),
     width = 9, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj,
        reduction = "wnn.umap",
        raster = F,
        cols = cols4all::c4a("rainbow", 15),
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
sc_obj_test <- sc_obj
DefaultAssay(sc_obj_test) <- "ADT"
sc_obj_test[["peaks"]] <- NULL
sc_obj_test[["chromvar"]] <- NULL
sc_obj_test[["RNA"]] <- NULL
sc_obj_test[["SCT"]] <- NULL
saveRDS(sc_obj_test, file = glue("{out_path}scWNN_obj_drcl.rds"))
saveRDS(sc_obj, file = glue("{out_path}scWNN_obj.rds"))
saveRDS(sc_obj@reductions, file = glue("{out_path}reducWNN_test.rds"))
saveRDS(sc_obj@commands, file = glue("{out_path}cmdWNN_test.rds"))

# sc_obj <- readRDS(glue("{out_path}scWNN_obj_drcl.rds"))
# colnames(sc_obj@meta.data)[grep("refine", colnames(sc_obj@meta.data))] <- "ann_level4_final"
# sc_meta <- cbind(sc_obj@meta.data,
#                  Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
# saveRDS(sc_meta, file = glue("{out_path}sc_meta.rds"))
# fwrite2(sc_meta, file = glue("{out_path}sc_meta.txt"), row.names = T)
# 
# saveRDS(sc_obj, file = glue("{out_path}scWNN_obj_drcl.rds"))
