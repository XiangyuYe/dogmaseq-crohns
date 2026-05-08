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
  "Prolif_features_list" = c("MKI67","TYMS","PCNA", "TOP2A"),
  "MAIT_features_list" = c("SLC4A10", "KIT"),
  "Naive_features_list" = c("CCR7","SELL","LEF1","TCF7"),
  "Memery_features_list" = c("GPR183","S100A4"),
  "Effect_features_list" = c("GZMK","GZMA", "GZMB", "GNLY","NKG7"),
  "Exhaust_features_list" = c("HAVCR2", "KLRK1", "LAG3"),
  "Th1_features_list" = c("CCL5", "CCR5", "IFNG"),
  "Th17_features_list" = c("RORC", "CCR6", "IL17A", "IL17F"),
  "HSP_features_list" = c("HSPA1A", "HSPA1B", "HSPH1"),
  "Tfh_features_list" = c("IL21", "BCL6", "CXCR5", "ICOS", "CXCL13"),
  "Treg_features_list" = c("FOXP3", "IKZF2", "IKZF1"),
  "IL_features_list" = c("IL4", "IL5", "IL9", "IL13", "IL2", "IL2RA", "IL2RB", "IL12RB2"),
  "TF_features_list" = c("CREM", "EGR1", "EGR2", "EGR3", "RUNX1", "RUNX2", "RUNX3",
                         "STAT1", "STAT3", "STAT4", "STAT5A", "STAT5B")
)
atac_mk <- c("CCR7", "CCL5", "TBX21", "RORC", "IL17A", "BCL6", "CXCR5", "FOXP3", "IKZF2")
names_motif <- c("TBX21", "RORC", "RORA", "STAT1", "STAT3")
sel_motif <- lapply(names_motif, function(x){
  namex <- getMatrixSet(
    JASPAR2020, 
    opts = list(name = x, 
                tax_group = "vertebrates", collection = "CORE")
  ) %>% names()
  return(namex)
}) %>% unlist

seed_use <- 20250528

##### load ADT, RNA, and ATAC data #####
data_path <- "03_output/03_clustering/CD4_CD103_TRM/"
out_path <- "03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/"
scRNA_obj <- readRDS(glue("{data_path}scRNA_obj_test.rds"))
sc_obj <- scRNA_obj
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

##### DRCL #####
eig_val <- (sc_obj@reductions[["harmony_SCT"]]@stdev)^2
var_explained <- eig_val / sum(eig_val)
cumsum(var_explained)

eig_val <- (sc_obj@reductions[["harmony_lsi"]]@stdev)^2
var_explained <- eig_val / sum(eig_val)
cumsum(var_explained)


sc_obj <- dr.cl.wnn(sc_obj = sc_obj,
                    redc_list = list("harmony_SCT", "harmony_lsi"),
                    dim_list = list(1:15, 1:30),
                    k_nn = 20,
                    prune_SNN = 1/20,
                    n_iter = 300,
                    res = 1,
                    cl_method = 1,
                    run_umap = T,
                    n_neig = 20L,
                    n_epochs = 500,
                    neg_rate = 50L,
                    sprd = 0.7,
                    min_dist = 0.4,
                    seed_use = seed_use)
##
tiff(file = glue("{out_path}UMAP_test2.tiff"),
     width = 7, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(sc_obj, 
        reduction = "wnn.umap", 
        raster = F, 
        label = T) +
  theme(title = element_blank())
dev.off()
##
sc_meta_ref <- readRDS("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/03_output/03_clustering_test/CD4_CD103_TRM/test302015/sc_meta.rds")
sc_obj$ann_ref <- sc_meta_ref[colnames(sc_obj) %>% 
                                gsub("_DOGMAseq\\-", "", .) %>%
                                gsub("Duerr_", "", .), 
                              "ann_level5_final"]
##
tiff(file = glue("{out_path}UMAP_test_ann_ref.tiff"),
     width = 9, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(subset(sc_obj, !is.na(ann_ref)),
        reduction = "wnn.umap",
        group.by = "ann_ref",
        cols = cd4cd103trm_cols,
        raster = F,
        label = F) +
  theme(title = element_blank())
dev.off()

## 
tiff(glue("{out_path}ft.tiff"), 
     height = 18, width = 34, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj, 
            reduction = "wnn.umap",
            features = c("rna_EGR1", "rna_EGR2", "rna_IKZF1", "rna_IL2", "IL12RB1", "IL12RB2", "IL17A", "TNF",
                         "rna_FOXP3", "rna_CCL5", "rna_CCR6", "rna_KIT", "rna_CREM",
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


sp_outpath <- "03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/"
e_regulon_name_trans <- readRDS(glue("{sp_outpath}e_regulon_name_trans.rds"))
scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
scAUC_obj_CD103 <- subset(scAUC_obj, cells = colnames(sc_obj))
sc_obj[["AUC"]] <- scAUC_obj_CD103[["AUC"]]
sel_er <- rownames(sc_obj[["AUC"]])[grep("ROR|IRF|CREM", rownames(sc_obj[["AUC"]]))]

dot_motif <- DotPlot(sc_obj, 
                     assay = "AUC", 
                     col.min = -2,
                     col.max = 2,
                     features = sel_er,
                     cols = c("lightgrey", "brown")) + 
  theme(axis.text.x = element_text(angle = 90))
tiff(glue("{out_path}dot_raw_AUC_test.tiff"), 
     height = 10, width = 7, units = "in", res = 300, compression = "lzw")
print(dot_motif)
dev.off()

tiff(glue("{out_path}ft_auc.tiff"), 
     height = 8, width = 17, units = "in", res = 300, compression = "lzw")
FeaturePlot(sc_obj, 
            reduction = "wnn.umap",
            features = sel_er, 
            raster = F,
            ncol = 4, 
            min.cutoff = "q1", 
            max.cutoff = "q99")
dev.off()
##### temp anotation #####
sc_obj$ann_level4_final <- "CD4_CD103_TRM_Th17"
sc_obj$ann_level4_final[sc_obj$wsnn_res.1 %in% c(5)] <- "CD4_CD103_TRM_Th1"
sc_obj$ann_level4_final[sc_obj$wsnn_res.1 %in% c(9, 11)] <- "CD4_CD103_TRM_Th17.1"

##### bulk heatmap #####
#### RNA ####
exclude_gene <- readRDS("03_output/02_clean/exclude_gene.rds")
DefaultAssay(sc_obj) <- "RNA"
Idents(sc_obj) <- sc_obj$seurat_clusters
ct_markers <- FindAllMarkers(sc_obj, 
                             features = setdiff(rownames(sc_obj), exclude_gene), 
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
pheatmap::pheatmap(sc_obj_bulk[["RNA"]]$data[unique(top5$gene), ],
                   scale = "row",
                   border_color = NA,
                   fontsize = 8,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()

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
sc_obj$ann_level5_temp <- paste0(sc_obj$ann_level4_final, 
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

##### save temp data #####
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
sc_obj$ann_level5_final <- sc_obj$ann_level4_final
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(0)] <- "CD4_CD103_TRM_Th17-HSP"
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(1)] <- "CD4_CD103_TRM_Th17-RUNX1"
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(2)] <- "CD4_CD103_TRM_Th17-NFKB"
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(3)] <- "CD4_CD103_TRM_Th17-CREM"
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(6)] <- "CD4_CD103_TRM_Th17-IKZF1"
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(7)] <- "CD4_CD103_TRM_Th17-MAIT"
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(8)] <- "CD4_CD103_TRM_Th17-EGR2"
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(10)] <- "CD4_CD103_TRM_Th17-EGR1"
sc_obj$ann_level5_final[sc_obj$wsnn_res.1 %in% c(11)] <- "CD4_CD103_TRM_Th17.1-KLR"
sc_obj$ann_level5_final <- factor(sc_obj$ann_level5_final,
                                  levels = c("CD4_CD103_TRM_Th1", 
                                             "CD4_CD103_TRM_Th17",
                                             "CD4_CD103_TRM_Th17-CREM",
                                             "CD4_CD103_TRM_Th17-EGR1",
                                             "CD4_CD103_TRM_Th17-EGR2",
                                             "CD4_CD103_TRM_Th17-HSP",
                                             "CD4_CD103_TRM_Th17-IKZF1",
                                             "CD4_CD103_TRM_Th17-MAIT",
                                             "CD4_CD103_TRM_Th17-NFKB",
                                             "CD4_CD103_TRM_Th17-RUNX1",
                                             "CD4_CD103_TRM_Th17.1",
                                             "CD4_CD103_TRM_Th17.1-KLR"))
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


