## load r/4.5.0
library(Seurat)
library(reticulate)
library(DropletUtils)
library(bigreadr)
library(glue)
library(dplyr)

##### data input and select parameters
wkdt <- "/ix1/rduerr/shared/rduerr_wchen/CD_annotation/"
setwd(wkdt)
# 1. DOGMA-seq CD4T
load(glue("CD4_final/CD4T_all_transferred_ann_removed.RData"))
count_dat <- GetAssayData(data,
                          assay = "RNA",
                          slot = "counts") %>%
  as.matrix() %>% as(., "CsparseMatrix")

write10xCounts(
  glue("CD4_final/gex"),
  count_dat,
  gene.type = "Gene Expression"
)
system(glue("mkdir -p {wkdt}CD4_final/TCAT/"))
system(glue("mkdir -p {wkdt}CD4_final/SCENIC+/"))
system(glue("mkdir -p {wkdt}CD4_final/meta/"))
saveRDS(data@reductions, file = glue("{wkdt}CD4_final/meta/reduction.rds"))
saveRDS(data@meta.data, file = glue("{wkdt}CD4_final/meta/meta_data.rds"))

#------ STARCAT RUNNING --------#
wkdt <- "/ix1/rduerr/shared/rduerr_wchen/CD_annotation/"

typex <- "CD4_final"
reduc_use <- "wnn.umap" 

tcat_path <- glue("{wkdt}/{typex}/TCAT/")

meta_data <- readRDS(glue("{wkdt}/{typex}/meta/meta_data.rds"))
reduct_obj <- readRDS(glue("{wkdt}/{typex}/meta/reduction.rds"))

usage_df <- fread2(glue("{tcat_path}usage.txt"))
score_df <- fread2(glue("{tcat_path}scores.txt"))
usage_df <- tibble::column_to_rownames(usage_df, "V1")
score_df <- tibble::column_to_rownames(score_df, "V1")

sc_test <- CreateSeuratObject(counts = t(usage_df), 
                              assay = "score", 
                              meta.data = cbind(meta_data, 
                                                score_df[rownames(usage_df), ]))
sc_test[["score"]]$data <- sc_test[["score"]]$counts
sc_test@reductions <- reduct_obj
saveRDS(sc_test, file = glue("{tcat_path}scScore_obj.rds"))
##
tiff(glue("{tcat_path}feature_score1.tiff"), 
     height = 15, width = 23, units = "in", compression = "lzw", res = 300)
FeaturePlot(sc_test, 
            reduction = reduc_use, raster = F, 
            features = c("CD4-Naive", "CD4-CM", 
                         "Th1-Like", "Th17-Resting", "Th17-Activated", 
                         "Th2-Activated", "Th2-Resting",
                         "Tfh-1", "Tfh-2", "Th22", "Tph", "Treg"),
            max.cutoff = "q99", min.cutoff = "q1",
            ncol = 4)
dev.off()

##
tiff(glue("{tcat_path}feature_score3.tiff"), 
     height = 15, width = 18, units = "in", compression = "lzw", res = 300)
FeaturePlot(sc_test, 
            reduction = reduc_use,raster = F, 
            features = c("CellCycle-G2M", "CellCycle-Late-S", "CellCycle-S", 
                         "Cytotoxic", "Exhaustion", "ISG",
                         "HLA", "Translation", "Mito"),
            max.cutoff = "q99", min.cutoff = "q1",
            ncol = 3)
dev.off()




