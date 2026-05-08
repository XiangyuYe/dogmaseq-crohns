## module load r/4.5.0
library(Seurat)
library(Signac)
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
scATAC_obj <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_scATAC_obj.rds")
sc_meta_cd4 <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta_ann_level3.5_final.rds")
scATAC_obj <- subset(scATAC_obj, cells = rownames(sc_meta_cd4)[sc_meta_cd4$ann_level3_final == "CD4_TRM"])

##### DRCL on all immune cells #####
seed_use <- 20250528
cutoff_q <- "q70"
scATAC_obj <- dr.cl.ATAC(scATAC_obj = scATAC_obj,
                         cutoff_q = cutoff_q,
                         use_assay = "peaks",
                         run_harmony = T,
                         batch_var = "Batch",
                         var_exp = 0.8,
                         res = 0.6,
                         run_umap = T,
                         seed_use = seed_use)
##
out_path <- "03_output/03_clustering/CD4_TRM/"
system(glue("mkdir -p {out_path}"))
tiff(file = glue("{out_path}UMAP_ATAC.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
print(DimPlot(scATAC_obj, 
              reduction = "umap", 
              label = T) + 
        NoLegend())
dev.off()

####
saveRDS(scATAC_obj, file = glue("{out_path}scATAC_obj_test_{cutoff_q}.rds"))
saveRDS(scATAC_obj@reductions, file = glue("{out_path}reduc_ATAC_test_{cutoff_q}.rds"))
saveRDS(scATAC_obj@commands, file = glue("{out_path}cmd_ATAC_test_{cutoff_q}.rds"))

