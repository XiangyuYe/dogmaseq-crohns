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

cnmf_path <- glue("03_output/03_clustering/CD4_TRM/cNMF/")
##### load RNA data #####
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
exclude_gene <- readRDS(file = "03_output/02_clean/exclude_gene.rds")
sc_meta_cd4 <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta_ann_level3.5_final.rds")
scRNA_obj <- subset(scRNA_obj, cells = rownames(sc_meta_cd4)[sc_meta_cd4$ann_level3_final == "CD4_TRM"])

expr_data <- scRNA_obj[["RNA"]]$counts
fwrite2(scRNA_obj@meta.data, file = glue("{cnmf_path}metadata.tsv"), row.names = T)
system(glue("rm -rf {cnmf_path}gex/"))
write10xCounts(
  glue("{cnmf_path}gex/"),
  expr_data,
  gene.type = "Gene Expression"
)

##### enrichment on GEP #####
spectra_tpm_k10 <- fread2(file = "03_output/03_clustering_test/CD4T/test153000/test302015/CD4_TRM/CD4_TRM_harmony_spectra_tpm_k10.txt",
                          header = T) %>%
  tibble::column_to_rownames(., "V1")
gene_bitr <- bitr(rownames(spectra_tpm_k10), 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db, 
                  drop = T)
spectra_tpm_k10 <- spectra_tpm_k10[gene_bitr$SYMBOL,]
##
gsea_reactome_df <- lapply(colnames(spectra_tpm_k10), function(gepx){
  
  geneListx <- spectra_tpm_k10[, gepx, drop = T]
  names(geneListx) <- gene_bitr$ENTREZID
  geneListx <- sort(geneListx, decreasing = T)
  gsea_reactome_dfx <- gsePathway(geneListx,
                                  organism = "human",
                                  minGSSize = 50,
                                  maxGSSize = 500,
                                  eps = 1e-10,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  verbose = T)@result
  if (nrow(gsea_reactome_dfx) > 0) {
    gsea_reactome_dfx$GEP <- gepx
  } else {
    gsea_reactome_dfx <- NULL
  }
  return(gsea_reactome_dfx)
  
}) %>% Reduce("rbind", .) %>% as.data.frame()
gsea_reactome_df <- cbind("DB" = "Reactome",
                          gsea_reactome_df)
saveRDS(gsea_reactome_df, file = "03_output/03_clustering_test/CD4T/test153000/test302015/CD4_TRM/gsea_reactome_k10.rds")
##
gsea_go_df <- lapply(colnames(spectra_tpm_k10), function(gepx){
  #
  geneListx <- spectra_tpm_k10[, gepx, drop = T]
  names(geneListx) <- gene_bitr$ENTREZID
  geneListx <- sort(geneListx, decreasing = T)
  #
  gsea_go_dfx <- gseGO(geneList = geneListx,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       minGSSize = 50,
                       maxGSSize = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       verbose = T)@result
  if (nrow(gsea_go_dfx) > 0) {
    gsea_go_dfx$GEP <- gepx
  } else {
    gsea_go_dfx <- NULL
  }
  return(gsea_go_dfx)
  
}) %>% Reduce("rbind", .) %>% as.data.frame()
# colnames(gsea_go_df)[1] <- "DB"
# gsea_go_df$DB <- paste0("GO-", gsea_go_df$DB)
gsea_go_df <- cbind("DB" = "GO-BP",
                    gsea_go_df)
saveRDS(gsea_go_df, file = "03_output/03_clustering_test/CD4T/test153000/test302015/CD4_TRM/gsea_go_k10.rds")

##
top_gene <- fread2("03_output/03_clustering_test/CD4T/test153000/test302015/CD4_TRM/CD4_TRM_harmony_top_genes_k10.txt")
top_gene_list <- lapply(1:ncol(top_gene), function(x) top_gene[,x] %>% unlist)
names(top_gene_list) <- colnames(top_gene)

source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/Enrichment.R")
enrich_file_path <- "03_output/03_clustering_test/CD4T/test153000/test302015/CD4_TRM/Enrichment/"
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

###
library(pheatmap)
library(bigreadr)
library(dplyr)
library(stringr)
library(Seurat)

project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/"
setwd(project_path)

source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/process_ATAC.R")
out_path <- "03_output/03_clustering_test/CD4_CD103_TRM/test302015/"
sc_obj <- readRDS(glue("{out_path}scWNN_obj.rds"))
nmf_usage7 <- fread2(file = "03_output/03_clustering_test/CD4T/test153000/test302015/CD103TRM/CD103TRM_harmony_usage_k7.txt",
                     header = T) %>%
  tibble::column_to_rownames(., "V1")
sc_obj[["NMF7"]] <- CreateAssayObject(t(nmf_usage7),   
                                      min.cells = 0,
                                      min.features = 0)

DefaultAssay(sc_obj) <- "NMF7"
tiff(glue("{out_path}feature_nmf7.tiff"), 
     height = 12, width = 14, units = "in", compression = "lzw", res = 300)
FeaturePlot(sc_obj,
            reduction = "wnn.umap",
            raster = F,
            features = colnames(nmf_usage7),
            max.cutoff = "q95", min.cutoff = "q5",
            ncol = 3)
dev.off()
##
dot_adt <- heat.adt(sc_obj = sc_obj,
                    clus_col = "seurat_clusters",
                    assay_use = "NMF7",
                    adt_use = colnames(nmf_usage7))
tiff(glue("{out_path}dot_score7.tiff"),
     height = 6, width = 4, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()
##
nmf_usage14 <- fread2(file = "03_output/03_clustering_test/CD4T/test153000/test302015/CD103TRM/CD103TRM_harmony_usage_k14.txt",
                      header = T) %>%
  tibble::column_to_rownames(., "V1")
sc_obj[["NMF14"]] <- CreateAssayObject(t(nmf_usage14),   
                                       min.cells = 0,
                                       min.features = 0)
##
DefaultAssay(sc_obj) <- "NMF14"
tiff(glue("{out_path}feature_nmf10.tiff"), 
     height = 16, width = 18, units = "in", compression = "lzw", res = 300)
FeaturePlot(sc_obj,
            reduction = "wnn.umap",
            raster = F,
            features = colnames(nmf_usage14),
            max.cutoff = "q95", min.cutoff = "q5",
            ncol = 4)
dev.off()
##
dot_adt <- heat.adt(sc_obj,
                    clus_col = "seurat_clusters",
                    assay_use = "NMF14",
                    adt_use = colnames(nmf_usage14))
tiff(glue("{out_path}dot_score10.tiff"),
     height = 6, width = 6, units = "in", res = 300, compression = "lzw")
dot_adt
dev.off()

##
sc_obj$condition <- factor(sc_obj$condition,
                           levels = c("N_N", "N_U", "I_I"))
DefaultAssay(sc_obj) <- "NMF"
tiff(glue("{out_path}dot_score_group.tiff"),
     height = 10, width = 6, units = "in", res = 300, compression = "lzw")
DotPlot(sc_obj, 
        features = colnames(nmf_usage14), 
        cols = c("orange", "blue", "red"), 
        # dot.scale = 8, 
        split.by = "condition") +
  RotatedAxis()
dev.off()


tiff("03_output/03_clustering_test/CD4T/test153000/test302015/CD103TRM/heat_test.tiff",
     height = 15, width = 25, units = "in",
     compression = "lzw", res = 300)
pheatmap::pheatmap(t(nmf_usage), 
                   scale = "column", 
                   show_colnames = F,
                   cluster_rows = T,
                   cluster_cols = T,
                   treeheight_row = 0, 
                   treeheight_col = 0)
dev.off()
