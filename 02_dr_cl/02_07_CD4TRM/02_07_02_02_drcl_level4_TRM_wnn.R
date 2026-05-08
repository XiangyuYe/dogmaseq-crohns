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
  "Th17_features_list" = c("RORC", "CCR6", "IL17A", "IL4", "IL5", "IL9", "IL13", "IL2"),
  "Tfh_features_list" = c("IL21", "BCL6", "CXCR5", "ICOS", "CXCL13"),
  "Treg_features_list" = c("FOXP3", "IKZF2")
)
atac_mk <- c("CCR7", "CCL5", "TBX21", "RORC", "IL17A", "BCL6", "CXCR5", "FOXP3", "IKZF2")
sel_motif <- c("MA0690.1", "MA1151.1", "MA0071.1", "MA0072.1", "MA0037.3", "MA0463.2", "MA0850.1")
names_motif <- c("TBX21", "RORC", "RORA.1", "RORA.2", "GATA3", "BCL6", "FOXP3")
seed_use <- 20250528

##### load ADT, RNA, and ATAC data #####
data_path <- "03_output/03_clustering/CD4_TRM/"
out_path <- "03_output/03_clustering/CD4_TRM/WNN_RNA_ATAC/"
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
sc_chromvar <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_scchromvar_assay.rds") %>%
  CreateSeuratObject(assay = "chromvar",
                     min.cells = 0,
                     min.features = 0)
sc_obj[["chromvar"]] <- subset(sc_chromvar, cells = colnames(sc_obj))[["chromvar"]]

##### DRCL #####
sc_obj <- dr.cl.wnn(sc_obj = sc_obj,
                    redc_list = list("harmony_SCT", "harmony_lsi"),
                    dim_list = list(1:15, 1:25),
                    k_nn = 20,
                    prune_SNN = 1/15,
                    n_iter = 300,
                    res = 0.5,
                    cl_method = 1,
                    run_umap = T,
                    n_neig = 20L,
                    n_epochs = 500,
                    neg_rate = 30L,
                    sprd = 0.8,
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

# ##
# sc_meta_ref <- readRDS("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/03_clustering_test/CD4T/test153000/test302015/sc_meta_refine0726.rds")
# sc_obj$ann_ref <- sc_meta_ref[colnames(sc_obj) %>%
#                                 gsub("_DOGMAseq\\-", "", .) %>%
#                                 gsub("Duerr_", "", .),
#                               "ann_level4_final"]
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
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj, 
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

##### bulk heatmap #####
DefaultAssay(sc_obj) <- "RNA"
Idents(sc_obj) <- sc_obj$seurat_clusters
ct_markers <- FindAllMarkers(sc_obj, 
                             only.pos = T)
##
ct_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
##
sc_obj_bulk <- Seurat::AverageExpression(sc_obj,
                                         features = unique(top5$gene),
                                         group.by = "seurat_clusters",
                                         assays = "RNA",
                                         return.seurat = T,
                                         layer = "data")
tiff(glue("{out_path}pseudo_heat.tiff"),
     height = 7, width = 5, units = "in", res = 300, compression = "lzw")
pheatmap::pheatmap(sc_obj_bulk[["RNA"]]$data[top5$gene, ],
                   scale = "row",
                   border_color = NA,
                   fontsize = 8,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()

##### anotation #####
sc_obj$ann_level4_final <- "CD4_TRM_Th17"
sc_obj$ann_level4_final[sc_obj$wsnn_res.0.6 %in% c(0, 8, 9)] <- "CD4_TRM_Th1"
sc_obj$ann_level4_final[sc_obj$wsnn_res.0.6 %in% c(2, 10, 11)] <- "CD4_TRM_Th17.1"
## 
tiff(file = glue("{out_path}UMAP_level4_final.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        group.by = "ann_level4_final",
        reduction = "wnn.umap", 
        label = T,
        raster = F) + NoLegend()
dev.off()
##
sc_obj$ann_level5_temp <- paste0(sc_obj$ann_level4_final , 
                                 "_ct", sc_obj$seurat_clusters)
## 
tiff(file = glue("{out_path}UMAP_level5_temp.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        group.by = "ann_level5_temp",
        reduction = "wnn.umap", 
        label = T,
        raster = F) + NoLegend()
dev.off()

##### save data #####
sc_meta <- cbind(sc_obj@meta.data,
                 Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta.rds"))
fwrite2(sc_meta, file = glue("{out_path}sc_meta.txt"), row.names = T)

##
saveRDS(sc_obj@reductions, file = glue("{out_path}reducWNN_test.rds"))
saveRDS(sc_obj@commands, file = glue("{out_path}cmdWNN_test.rds"))

sc_obj_test <- sc_obj
DefaultAssay(sc_obj_test) <- "ADT"
sc_obj_test[["ATAC"]] <- NULL
sc_obj_test[["RNA"]] <- NULL
sc_obj_test[["SCT"]] <- NULL
saveRDS(sc_obj_test, file = glue("{out_path}scWNN_obj_drcl.rds"))
saveRDS(sc_obj, file = glue("{out_path}scWNN_obj.rds"))

# sc_obj <- readRDS(glue("{out_path}scWNN_obj_drcl.rds"))
# colnames(sc_obj@meta.data)[grep("refine", colnames(sc_obj@meta.data))] <- "ann_level4_final"
# sc_meta <- cbind(sc_obj@meta.data,
#                  Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
# saveRDS(sc_meta, file = glue("{out_path}sc_meta.rds"))
# fwrite2(sc_meta, file = glue("{out_path}sc_meta.txt"), row.names = T)
# 
# saveRDS(sc_obj, file = glue("{out_path}scWNN_obj_drcl.rds"))

