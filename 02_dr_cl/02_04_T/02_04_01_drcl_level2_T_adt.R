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

##### load ADT data #####
scADT_obj <- readRDS("03_output/02_clean/scADT_obj_immune.rds")
sc_meta_all <- readRDS("03_output/03_clustering/all/WNN_ADT_RNA/sc_meta_ann_level2.rds")
scADT_obj <- subset(scADT_obj, cells = rownames(sc_meta_all)[sc_meta_all$ann_level2 == "T"])
##### DRCL on all T cells #####
seed_use <- 20250528
pc_clust_adt <- 20
resx <- 0.6
ctrl_pt <- rownames(scADT_obj)[grep("Ctrl", rownames(scADT_obj))]
VariableFeatures(scADT_obj) <- setdiff(rownames(scADT_obj), ctrl_pt)
scADT_obj <- ScaleData(scADT_obj) %>% 
  RunPCA() %>%
  RunHarmony(.,
             group.by.vars = "Batch", 
             reduction.use = "pca",
             reduction.save = "harmony_adt",
             assay.use = "ADT",
             project.dim = F)

scADT_obj <- FindNeighbors(scADT_obj, 
                           reduction = "harmony_adt",
                           n.trees = 300,
                           k.param = 20, 
                           dims = 1:pc_clust_adt)
scADT_obj <- FindClusters(scADT_obj, 
                          n.iter = 300,
                          resolution = resx,
                          random.seed = seed_use)
scADT_obj <- RunUMAP(scADT_obj, 
                     dims = 1:pc_clust_adt, 
                     reduction = 'harmony_adt',
                     n.neighbors = 30L,
                     umap.method = "uwot",
                     n.epochs = 300,
                     negative.sample.rate = 10L,
                     min.dist = 0.4,
                     seed.use = seed_use)
## plot
out_path <- "03_output/03_clustering/T/"
system(glue("mkdir -p {out_path}"))
tiff(file = glue("{out_path}UMAP_ADT.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(scADT_obj, 
        reduction = "umap", 
        raster = F,
        label = T) + NoLegend()
dev.off()

##
adt_mk <- c("CD45", "CD3","CD4", "CD8", "CD161", 
            "HLA-DR", "CD64", "CD123", "CD19")
heat_adt <- heat.adt(sc_obj = scADT_obj,
                     clus_col = "seurat_clusters",
                     assay_use = "ADT",
                     adt_use = adt_mk)
tiff(glue("{out_path}dot_adt.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
heat_adt
dev.off()

##### output #####
saveRDS(scADT_obj, file = glue("{out_path}scADT_obj.rds"))
saveRDS(scADT_obj@reductions, file = glue("{out_path}reduc_adt.rds"))
saveRDS(scADT_obj@commands, file = glue("{out_path}cmd_adt.rds"))

