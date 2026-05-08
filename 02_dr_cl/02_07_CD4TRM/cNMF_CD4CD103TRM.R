##### cNMF #####
library(Seurat)
library(bigreadr)
library(dplyr)
library(glue)
library(stringr)
library(DropletUtils)
library(Matrix)

project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

data_path <- glue("03_output/03_clustering/CD4_CD103_TRM/")
##### load RNA data #####
sc_obj <- readRDS(glue("{data_path}WNN_RNA_ATAC/scWNN_obj.rds"))
expr_data <- sc_obj[["RNA"]]$counts
expr_data <- expr_data[rownames(sc_obj[["SCT"]]),]
fwrite2(sc_obj@meta.data, file = glue("{data_path}cNMF/metadata.tsv"), row.names = T)
system(glue("rm -rf {data_path}cNMF/gex/"))
write10xCounts(
  glue("{data_path}cNMF/gex/"),
  expr_data,
  gene.type = "Gene Expression"
)
##
sc_obj_th17 <- subset(sc_obj, ann_level4_final == "CD4_CD103_TRM_Th17")
expr_data <- sc_obj_th17[["RNA"]]$counts
expr_data <- expr_data[rownames(sc_obj_th17[["SCT"]]),]
fwrite2(sc_obj_th17@meta.data, file = glue("{data_path}cNMF_Th17/metadata.tsv"), row.names = T)
system(glue("rm -rf {data_path}cNMF_Th17/gex/"))
write10xCounts(
  glue("{data_path}cNMF_Th17/gex/"),
  expr_data,
  gene.type = "Gene Expression"
)
##################################################
library(pheatmap)
library(bigreadr)
library(dplyr)
library(stringr)

project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/"
setwd(project_path)
##
##### check cNMF results #####
nmf_usage11 <- fread2(file = "03_output/03_clustering/CD4_CD103_TRM/cNMF/CD4_CD103_TRM_harmony_usage_k11.txt",
                      header = T) %>%
  tibble::column_to_rownames(., "V1")
sc_obj[["NMF11"]] <- CreateAssayObject(t(nmf_usage11),   
                                       min.cells = 0,
                                       min.features = 0)
DefaultAssay(sc_obj) <- "NMF11"
tiff(glue("{out_path}feature_NMF11_rna.tiff"), 
     height = 12, width = 18, units = "in", compression = "lzw", res = 300)
FeaturePlot(sc_obj,
            reduction = "wnn.umap",
            raster = F,
            features = colnames(nmf_usage11),
            max.cutoff = "q95", min.cutoff = "q5",
            ncol = 4)
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    clus_col = "seurat_clusters",
                    assay_use = "NMF11",
                    adt_use = colnames(nmf_usage11))
tiff(glue("{out_path}dot_score11.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

##
top_gene <- fread2("03_output/03_clustering/CD4_CD103_TRM/cNMF/CD4_CD103_TRM_harmony_top_genes_k11.txt")
top_gene_list <- lapply(1:ncol(top_gene), function(x) top_gene[,x] %>% unlist)
names(top_gene_list) <- colnames(top_gene)

source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/Enrichment.R")
enrich_file_path <- "03_output/03_clustering/CD4_CD103_TRM/cNMF/Enrichment/"
system(glue("mkdir -p {enrich_file_path}"))
lapply(seq_along(top_gene_list), function(x){
  Enrichment.pipline(gene_list = top_gene_list[[x]],
                     list_name = paste0("Enrichment_", names(top_gene_list)[x]),
                     num_show = 10,
                     plot = T,
                     GO = T,
                     KEGG = T,
                     Reactome = T,
                     outpath = enrich_file_path)
})
