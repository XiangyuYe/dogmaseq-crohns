## module load r/4.5.0
library(Seurat)
library(Signac)
library(dplyr)
library(glue)

## set path
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

## load data
scATAC_obj <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_addmotif.rds")
sc_meta_CD4 <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta.rds")
scATAC_obj <- subset(scATAC_obj, cells = rownames(sc_meta_CD4))
scATAC_obj@meta.data <- sc_meta_CD4

## load DAPs
DefaultAssay(scATAC_obj) <- "peaks"
Idents(scATAC_obj) <- scATAC_obj$ann_level4_final
dap_cd4 <- readRDS("03_output/04_Diff/DAR/TI/all/DESeq2/pseudo_bulk_ann_level2_refine_II_vs_NN.rds")[["CD4T"]]
dap_cd4_nu <- readRDS("03_output/04_Diff/DAR/TI/all/DESeq2/pseudo_bulk_ann_level2_refine_NU_vs_NN.rds")[["CD4T"]]

#### II vs NN ####
sig_up_peak <- dap_cd4$Term[which(dap_cd4$padj < 0.05 & dap_cd4$log2FC > 0.25)]
sig_down_peak <- dap_cd4$Term[which(dap_cd4$padj < 0.05 & dap_cd4$log2FC < -0.25)]

# find peaks open
open_peaks <- AccessiblePeaks(scATAC_obj)

# match the overall GC content in the peak set
meta_feature <- GetAssayData(scATAC_obj, 
                             assay = "peaks", 
                             layer = "meta.features")
peaks_matched_up <- MatchRegionStats(
  meta.feature = meta_feature[open_peaks, ],
  query.feature = meta_feature[sig_up_peak, ],
  features.match = c("percentile"),
  n = 50000)
enriched_motif_up <- FindMotifs(scATAC_obj,
                                background = peaks_matched_up,
                                features = sig_up_peak)
rownames(enriched_motif_up) <- NULL
enriched_motif_up$cluster <- "CD4T"
enriched_motif_up$dir <- "Up"

# match the overall GC content in the peak set
peaks_matched_down <- MatchRegionStats(
  meta.feature = meta_feature[open_peaks, ],
  query.feature = meta_feature[sig_down_peak, ],
  features.match = c("percentile"),
  n = 50000)
enriched_motif_down <- FindMotifs(scATAC_obj,
                                   background = peaks_matched_down,
                                   features = sig_down_peak)
rownames(enriched_motif_down) <- NULL
enriched_motif_down$cluster <- "CD4T"
enriched_motif_down$dir <- "Down"

##
motif_enrich_df <- rbind(enriched_motif_up, enriched_motif_down)
motif_enrich_df$Term <- motif_enrich_df$motif
motif_enrich_df$motif <- motif_enrich_df$motif.name
saveRDS(motif_enrich_df, file = "03_output/04_Diff/DAR/TI/all/DESeq2/motif_enrich_df_CD4T_II_vs_NN.rds")

#### NU vs NN ####
sig_up_peak <- dap_cd4_nu$Term[which(dap_cd4_nu$padj < 0.05 & dap_cd4_nu$log2FC > 0.25)]
sig_down_peak <- dap_cd4_nu$Term[which(dap_cd4_nu$padj < 0.05 & dap_cd4_nu$log2FC < -0.25)]

# match the overall GC content in the peak set
peaks_matched_up <- MatchRegionStats(
  meta.feature = meta_feature[open_peaks, ],
  query.feature = meta_feature[sig_up_peak, ],
  features.match = c("percentile"),
  n = 50000)
enriched_motif_up <- FindMotifs(scATAC_obj,
                                background = peaks_matched_up,
                                features = sig_up_peak)
rownames(enriched_motif_up) <- NULL
enriched_motif_up$cluster <- "CD4T"
enriched_motif_up$dir <- "Up"

# match the overall GC content in the peak set
peaks_matched_down <- MatchRegionStats(
  meta.feature = meta_feature[open_peaks, ],
  query.feature = meta_feature[sig_down_peak, ],
  features.match = c("percentile"),
  n = 50000)
enriched_motif_down <- FindMotifs(scATAC_obj,
                                  background = peaks_matched_down,
                                  features = sig_down_peak)
rownames(enriched_motif_down) <- NULL
enriched_motif_down$cluster <- "CD4T"
enriched_motif_down$dir <- "Down"

##
motif_enrich_df <- rbind(enriched_motif_up, enriched_motif_down)
motif_enrich_df$Term <- motif_enrich_df$motif
motif_enrich_df$motif <- motif_enrich_df$motif.name
saveRDS(motif_enrich_df, file = "03_output/04_Diff/DAR/TI/all/DESeq2/motif_enrich_df_CD4T_NU_vs_NN.rds")
