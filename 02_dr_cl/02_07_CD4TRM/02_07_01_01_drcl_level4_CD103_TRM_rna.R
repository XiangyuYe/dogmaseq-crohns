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

##### load RNA data #####
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
exclude_gene <- readRDS(file = "03_output/02_clean/exclude_gene.rds")
sc_meta_cd4 <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta_ann_level3.5_final.rds")
scRNA_obj <- subset(scRNA_obj, cells = rownames(sc_meta_cd4)[sc_meta_cd4$ann_level3_final == "CD4_CD103_TRM"])

##### DRCL on all immune cells #####
seed_use <- 20250528
pc_clust_gex <- 20
resx <- 0.6
scRNA_obj <- SCT.harmony(scRNA_obj = scRNA_obj,
                         batch_col = "Batch",
                         n_higvar = 3000,
                         exclude_gene = exclude_gene,
                         seed_use = seed_use)
##
scRNA_obj <- FindNeighbors(scRNA_obj, 
                           reduction = "harmony_SCT",
                           n.trees = 300,
                           k.param = 20, 
                           dims = 1:pc_clust_gex)
scRNA_obj <- FindClusters(scRNA_obj, 
                          n.iter = 300,
                          resolution = resx,
                          random.seed = seed_use)
scRNA_obj <- RunUMAP(scRNA_obj, 
                     dims = 1:pc_clust_gex, 
                     reduction = "harmony_SCT",
                     n.neighbors = 30L,
                     umap.method = "uwot",
                     n.epochs = 300,
                     negative.sample.rate = 10L,
                     min.dist = 0.4,
                     seed.use = seed_use)

out_path <- "03_output/03_clustering/CD4_CD103_TRM/"
system(glue("mkdir -p {out_path}"))
tiff(file = glue("{out_path}UMAP_RNA.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
print(DimPlot(scRNA_obj, 
              reduction = "umap", 
              raster = F,
              label = T) + NoLegend())
dev.off()

##
saveRDS(scRNA_obj, file = glue("{out_path}scRNA_obj_test.rds"))
saveRDS(scRNA_obj@meta.data, file = glue("{out_path}scRNA_meta_test.rds"))
saveRDS(scRNA_obj[["SCT"]], file = glue("{out_path}scSCT_obj_test.rds"))
saveRDS(scRNA_obj@reductions, file = glue("{out_path}reducRNA_test.rds"))
saveRDS(scRNA_obj@commands, file = glue("{out_path}cmdRNA_test.rds"))

