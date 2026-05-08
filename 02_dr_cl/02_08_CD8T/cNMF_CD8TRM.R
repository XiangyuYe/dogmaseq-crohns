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

cnmf_path <- glue("03_output/03_clustering/CD8_TRM/cNMF/")
##### load RNA data #####
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
exclude_gene <- readRDS(file = "03_output/02_clean/exclude_gene.rds")
sc_meta_cd8 <- readRDS("03_output/03_clustering/CD8T/WNN_ADT_RNA/sc_meta_ann_level3_final.rds")
scRNA_obj <- subset(scRNA_obj, cells = rownames(sc_meta_cd8)[sc_meta_cd8$ann_level3_final == "CD8_TRM"])

expr_data <- scRNA_obj[["RNA"]]$counts
fwrite2(scRNA_obj@meta.data, file = glue("{cnmf_path}metadata.tsv"), row.names = T)
system(glue("rm -rf {cnmf_path}gex/"))
write10xCounts(
  glue("{cnmf_path}gex/"),
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
     height = 8, width = 13, units = "in", compression = "lzw", res = 300)
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
top_gene <- fread2("03_output/03_clustering/CD4_CD103_TRM/cNMF/CD103TRM_harmony_top_genes_k11.txt")
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
