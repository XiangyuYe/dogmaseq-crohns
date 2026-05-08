#
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
library(cols4all)
##
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/"
setwd(project_path)
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/process_ATAC.R")

adt_mk <- c("CD4", "CD8","CD45RA", "CD45RO", "CD62L", "CD49a", "CD69", "CX3CR1",
            "CD183", "CD195", "CD103", "CD161", "KLRG1", "CD16", "CD39",
            "CD194", "CD196", "FceRIa",
            "CD279", "CD185", "CD25", "CD127", "TCR-Vd2")
# FceRIa
# [161] "TCR-ab"                   "TCR-Va7.2"               
# [163] "TCR-Vd2"                  "TIGIT"    
feature_list <- list(
  "T_features_list" = c("CD8A", "CD8B", "TRDC", "TRAC", "TRGC1", 
                        "KIT", "SLC4A10", "IL26", "IL23R", "SOX4"),
  "Prolif_features_list" = c("MKI67","TYMS","PCNA"),
  "Naive_features_list" = c("CCR7","SELL","LEF1","TCF7", "KLF2"),
  "Memery_features_list" = c("GPR183","S100A4", "TNFRSF9", "TGFBR3", 
                             "FCRL3", "FCRL4", "FCRL5",
                             "CCL3", "CCL4",  "CCL20", "IL7R", "CCL5", "CCR6", "RORA"),
  "Effect_features_list" = c("GZMK", "GZMA", "GZMB", "GZMH", "KLRF1", "IFNG", "IL17A", "IL17F"),
  "Exhaust_features_list" = c("HAVCR2", "LAG3", "TIGIT", "IKZF2", "PDCD1")
)
atac_mk <- c("CCR7", "CCL5", "TBX21", "RORC", "IL17A", "BCL6", "CXCR5", "FOXP3", "IKZF2")
sel_motif <- c("MA0690.1", "MA1151.1", "MA0071.1", "MA0072.1", "MA0050.2", "MA0517.1", "MA0471.2")
names_motif <- c("TBX21", "RORC", "RORA.1", "RORA.2", "IRF1", "STAT1_2", "E2F6")
seed_use <- 20250528

# MA0050.2     0   4.296137 0.686 0.299         0
# MA0772.1     0   2.976167 0.766 0.421         0
# MA0652.1     0   2.988289 0.735 0.391         0
# MA0517.1     0   3.970383 0.614 0.274         0
# MA0051.1     0   2.230429 0.686 0.360         0
# MA0653.1     0   3.043729 0.747 0.427         0
# MA1596.1     0   3.234491 0.718 0.398         0
# MA1419.1     0   2.581840 0.719 0.407         0
# MA0471.2     0   3.788488 0.645 0.335         0
# MA1125.1     0   2.474209 0.723 0.416         0

####################################
data_path <- "03_output/03_clustering_test/gdT/"
out_path <- "03_output/03_clustering_test/gdT/"
sc_obj <- readRDS(glue("{data_path}scRNA_obj_test.rds"))
sc_obj$wnnUMAP_1 <- NULL
sc_obj$wnnUMAP_2 <- NULL
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)
#
scADT_obj <- readRDS(glue("{data_path}scADT_obj_test.rds"))
sc_obj[["ADT"]] <- scADT_obj[["ADT"]]
sc_obj@reductions$harmony_adt <- scADT_obj@reductions$harmony_adt
sc_obj@reductions$umap_adt <- scADT_obj@reductions$umap
#
scATAC_obj <- readRDS(glue("{data_path}scATAC_obj_test.rds"))
sc_obj[["ATAC"]] <- scATAC_obj[["ATAC"]]
sc_obj@reductions$harmony_lsi <- scATAC_obj@reductions$harmony_lsi
sc_obj@reductions$umap_lsi <- scATAC_obj@reductions$umap_lsi
#
sc_chromvar <- readRDS("03_output/03_clustering_test/scchromvar_assay.rds") %>%
  CreateSeuratObject(assay = "chromvar",
                     min.cells = 0, 
                     min.features = 0)
sc_obj[["chromvar"]] <- subset(sc_chromvar, cells = colnames(sc_obj))[["chromvar"]]
##
sc_obj <- dr.cl.wnn(sc_obj = sc_obj,
                    redc_list = list("harmony_adt", "harmony_SCT", "harmony_lsi"),
                    dim_list = list(1:10, 1:20, 1:30),
                    k_nn = 20,
                    prune_SNN = 1/20,
                    n_iter = 300,
                    res = 2,
                    cl_method = 1,
                    run_umap = T,
                    n_neig = 20L,
                    n_epochs = 500,
                    neg_rate = 20L,
                    sprd = 0.5,
                    min_dist = 0.2,
                    seed_use = seed_use)
sc_obj <- FindClusters(sc_obj,
                       n.iter = 300,
                       algorithm = 1,
                       graph.name = "wsnn",
                       resolution = 2,
                       verbose = T,
                       random.seed = seed_use)
##
tiff(file = glue("{out_path}UMAP_test.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        # group.by = "wsnn_res.1",
        raster = F, 
        label = T) +
  theme(title = element_blank())
dev.off()
##
sc_obj$section_group <- paste0(sc_obj$section_comb, ":", sc_obj$condition)
tiff(file = glue("{out_path}UMAP_test_split.tiff"),
     width = 12, height = 8, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        split.by = "section_group",
        raster = F, 
        label = T, ncol = 3) +
  NoLegend()
dev.off()

####
tiff(glue("{out_path}ft.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj, 
            reduction = "wnn.umap",
            features = c("adt_CD45RA", "adt_CD127", "adt_CD103", "adt_CD62L", 
                         "adt_CD69", "adt_CD49a", "adt_TCR-Vd2", "rna_IKZF2",
                         "rna_TCF7", "rna_KLRF1", "rna_TRDV2", "rna_GZMH", "rna_GZMB", 
                         "rna_GZMK", "chromvar_MA1151.1"), 
            raster = F,
            ncol = 5, 
            min.cutoff = "q5", 
            max.cutoff = "q95")
dev.off()



tiff(glue("{out_path}dot_rna_test.tiff"), 
     height = 10, width = 12, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        assay = "RNA", 
        col.min = -2, 
        col.max = 2, 
        group.by = "wsnn_res.2",
        features = unlist(feature_list) %>% unique) + 
  theme(axis.text.x = element_text(angle = 90))
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    assay_use = "ADT",
                    adt_use = adt_mk)
tiff(glue("{out_path}dot_adt_test.tiff"),
     height = 6, width = 6, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
##
dot_motif <- DotPlot(sc_obj, 
                     assay = "chromvar", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_motif,
                     group.by = "wsnn_res.2",
                     cols = c("lightgrey", "brown")) + 
  scale_x_discrete(labels = names_motif) +
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_motif_test.tiff"), 
     height = 10, width = 5, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

## annotation
sc_obj$ann_level3_temp1 <- "gdT_CD103_Vd1"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(14)] <- "gdT_CD103_Vd2"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(7, 18)] <- "gdT_Vd2"
sc_obj$ann_level3_temp1[sc_obj$wsnn_res.2 %in% c(15)] <- "gdT_Vd1"

## annotation 2
sc_obj$ann_level4_temp1 <- sc_obj$ann_level3_temp1
#
sc_obj$ann_level4_temp1[sc_obj$wsnn_res.2 %in% c(0, 5, 10, 17, 21)] <- "gdT_CD103_Vd1_IRF1"
sc_obj$ann_level4_temp1[sc_obj$wsnn_res.2 %in% c(8, 9)] <- "gdT_CD103_Vd1_GZMA"
sc_obj$ann_level4_temp1[sc_obj$wsnn_res.2 %in% c(15)] <- "gdT_Vd1_GZMH"
sc_obj$ann_level4_temp1[sc_obj$wsnn_res.2 %in% c(12, 13, 19)] <- "gdT_CD103_Vd1_TCF7"
sc_obj$ann_level4_temp1[sc_obj$wsnn_res.2 %in% c(7, 18)] <- "gdT_Vd2_GZMK"

##### save #####
##
sc_obj@meta.data$ann_level4_final <- sc_obj@meta.data$ann_level4_temp1
Idents(sc_obj) <- sc_obj@meta.data$ann_level4_final
sc_obj$ann_level4_final <- factor(sc_obj$ann_level4_final,
                                  levels = c("gdT_Vd1_GZMH", 
                                             "gdT_CD103_Vd1", "gdT_CD103_Vd1_TCF7", 
                                             "gdT_CD103_Vd1_IRF1", "gdT_CD103_Vd1_GZMA", 
                                             "gdT_Vd2_GZMK", "gdT_CD103_Vd2"))
Idents(sc_obj) <- sc_obj@meta.data$ann_level4_final
tiff(file = glue("{out_path}UMAP_ann_level4.tiff"),
     width = 9, height = 6.5, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F,
        cols = cols4all::c4a("rainbow", 7),
        label = F)
dev.off()
##
sc_obj@meta.data$ann_level3_final <- factor(sc_obj@meta.data$ann_level3_temp1,
                                            levels = c("gdT_Vd1", "gdT_CD103_Vd1",
                                                       "gdT_Vd2", "gdT_CD103_Vd2"))
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

# ##
# Idents(sc_obj) <- sc_obj@meta.data$ann_level4_final
# sc_obj[["peaks"]] <- sc_obj[["recall_peaks"]]
# cr_plt <- CoveragePlot(
#   object = sc_obj,
#   assay = "peaks",
#   region = "CCL5",
#   annotation = T,
#   peaks = T, 
#   extend.downstream = 1000,
#   extend.upstream = 1000
# )
# tiff(glue("{out_path}cr_plt_level4_final.tiff"), 
#      height = 10, width = 6, units = "in", res = 300, compression = "lzw")
# print(cr_plt)
# dev.off()

##
sc_meta <- cbind(sc_obj@meta.data,
                 Embeddings(sc_obj, reduction = "wnn.umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta.rds"))
#
saveRDS(sc_obj@reductions, file = glue("{out_path}reducWNN_test.rds"))
saveRDS(sc_obj@commands, file = glue("{out_path}cmdWNN_test.rds"))

sc_obj_test <- sc_obj
DefaultAssay(sc_obj_test) <- "ADT"
sc_obj_test[["ATAC"]] <- NULL
sc_obj_test[["RNA"]] <- NULL
sc_obj_test[["SCT"]] <- NULL
saveRDS(sc_obj_test, file = glue("{out_path}scWNN_obj_drcl.rds"))
saveRDS(sc_obj, file = glue("{out_path}scWNN_obj.rds"))

##
DefaultAssay(sc_obj) <- "RNA"
Idents(sc_obj) <- sc_obj$ann_level4_temp1
ct_markers <- FindAllMarkers(sc_obj, 
                             only.pos = T)
##
ct_markers <- ct_markers[order(ct_markers$cluster),]
ct_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  # dplyr::filter(pct.1 > 0.3) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

# sc_obj <- ScaleData(sc_obj)
heat_plt <- DoHeatmap(sc_obj, features = top5$gene) + NoLegend()
ggsave(glue("{out_path}heat_test_ann_level4_temp1.png"), 
       heat_plt, width = 15, height = 10, units = "in", dpi = 300)


sc_obj_bulk <- Seurat::AverageExpression(sc_obj,
                                         features = unique(top5$gene),
                                         group.by = "ann_level4_temp1",
                                         assays = "RNA",
                                         return.seurat = T,
                                         layer = "data")
tiff(glue("{out_path}pseudo_heat_test_ann_level4_temp1.tiff"),
     height = 10, width = 5, units = "in", res = 300, compression = "lzw")
pheatmap::pheatmap(sc_obj_bulk[["RNA"]]$data[unique(top5$gene),],
                   scale = "row",
                   border_color = NA,
                   fontsize = 8,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()

