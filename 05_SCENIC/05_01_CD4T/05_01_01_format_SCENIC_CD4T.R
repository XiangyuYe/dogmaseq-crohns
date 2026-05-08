## module load r/4.5.0
library(Seurat)
library(dplyr)
library(glue)
## set work dir
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

##### export matrix for SCENIC+ #####
ctx <- "CD4T"
data_path <- glue("03_output/03_clustering/{ctx}/WNN_ADT_RNA/")
scenic_path <- glue("03_output/07_SCENIC/{ctx}/")
# sc_obj <- readRDS(glue("{data_path}scWNN_obj.rds"))
##
source("code/FUNCTION/scenicplus/export_scenicplus.R")
export.scenic(seurat_obj = sc_obj,
              assay_rna = "RNA",
              assay_atac = "peaks",
              outpath = scenic_path)

##### DA analysis #####
## load data
data_path <- "03_output/03_clustering/CD4T/"
region_set_path <- "03_output/07_SCENIC/CD4T/outs/region_sets/"
#
cutoff_q <- "q70"
scATAC_obj <- readRDS(glue("{data_path}scATAC_obj_test_{cutoff_q}.rds"))
scATAC_obj@meta.data <- readRDS(glue("{data_path}WNN_ADT_RNA/sc_meta.rds"))

## level3
DefaultAssay(scATAC_obj) <- "peaks"
Idents(scATAC_obj) <- scATAC_obj$ann_level3_final
dap_level3 <- FindAllMarkers(scATAC_obj, only.pos = T)
saveRDS(dap_level3, file = glue("{data_path}dap_level3.rds"))
## level4 group by level3
dap_level4_by3 <- SplitObject(scATAC_obj, split.by = "ann_level3_final") %>%
  lapply(., function(scATAC_objx){
    
    Idents(scATAC_objx) <- scATAC_objx$ann_level4_final
    if (length(unique(scATAC_objx$ann_level4_final)) > 1) {
      dap_level4by3x <- FindAllMarkers(scATAC_objx, only.pos = T)
    } else {
      dap_level4by3x <- NULL
    }
    return(dap_level4by3x)
  }) %>% Reduce("rbind", .)
saveRDS(dap_level4_by3, file = glue("{data_path}dap_level4_by3.rds"))

##### output as region set for SCENIC+ #####
dap_level3 <- readRDS(glue("{data_path}dap_level3.rds"))
dap_level4_by3 <- readRDS(glue("{data_path}dap_level4_by3.rds"))
#
system(glue("mkdir -p {region_set_path}DAP_level3/"))
lapply(unique(dap_level3$cluster), function(ctx){
  
  dap_level3x <- subset(dap_level3, cluster == ctx & 
                          p_val_adj < 0.05 & 
                          avg_log2FC > 0.25)
  dap_level3x <- dap_level3x[order(dap_level3x$avg_log2FC, -dap_level3x$p_val_adj, decreasing = T),]
  if (nrow(dap_level3x) > 0) {
    dapx <- stringr::str_split(dap_level3x$gene, "\\-", simplify = T) %>% as.data.frame()
    write.table(dapx, file = glue("{region_set_path}DAP_level3/{ctx}.bed"), 
                sep = "\t", row.names = F, col.names = F, quote = F)
  }
  return(glue("{ctx}: {nrow(dap_level3x)}"))
})
##
system(glue("mkdir -p {region_set_path}DAP_level4_by3/"))
lapply(unique(dap_level4_by3$cluster), function(ctx){
  
  dap_level4_by3x <- subset(dap_level4_by3, cluster == ctx & 
                              p_val_adj < 0.05 & 
                              avg_log2FC > 0.25)
  dap_level4_by3x <- dap_level4_by3x[order(dap_level4_by3x$avg_log2FC, -dap_level4_by3x$p_val_adj, decreasing = T),]
  if (nrow(dap_level4_by3x) > 0) {
    dapx <- stringr::str_split(dap_level4_by3x$gene, "\\-", simplify = T) %>% as.data.frame()
    write.table(dapx, file = glue("{region_set_path}DAP_level4_by3/{ctx}.bed"), 
                sep = "\t", row.names = F, col.names = F, quote = F)
  }
  return(glue("{ctx}: {nrow(dap_level4_by3x)}"))
})

