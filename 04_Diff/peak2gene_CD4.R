library(Seurat)
library(Signac)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)

###data input
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

## load data
sc_obj <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/scWNN_obj.rds")
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)
sc_obj@meta.data <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta.rds")

##### peak2gene: CD4T #####
ann_col <- "ann_level2_refine"
deg_path_it <- glue("03_output/04_Diff/DEG/TI/all/DESeq2/")
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4T"]]
deg_cd4 <- deg_df_it1$Term[which(deg_df_it1$padj < 0.05 & abs(deg_df_it1$log2FC) > 0.25)]

# first compute the GC content for each peak
DefaultAssay(sc_obj) <- "peaks"
sc_obj <- RegionStats(sc_obj, 
                      genome = BSgenome.Hsapiens.UCSC.hg38)
# link peaks to genes
sc_obj <- LinkPeaks(
  object = sc_obj,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = deg_cd4
)
links_gr <- Links(sc_obj[["peaks"]])
system(glue("mkdir -p 03_output/04_Diff/peak2gene/TI/"))
saveRDS(links_gr, file = "03_output/04_Diff/peak2gene/TI/LinkPeaks_links_CD4T_DEG_II_vs_NN.rds")

##### peak2gene: CD103_TRM-Th17 #####
use_de <- "DESeq2"
ann_col <- "ann_level4_final"
deg_path_it <- glue("03_output/04_Diff/DEG/TI/CD4T/DESeq2/")
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4_CD103_TRM_Th17"]]
deg_cd103trm17 <- deg_df_it1$Term[which(deg_df_it1$padj < 0.05 & abs(deg_df_it1$log2FC) > 0.25)]

sc_obj_cd103trm17 <- subset(sc_obj, ann_level4_final == "CD4_CD103_TRM_Th17")
# first compute the GC content for each peak
DefaultAssay(sc_obj_cd103trm17) <- "peaks"
sc_obj_cd103trm17 <- RegionStats(sc_obj_cd103trm17, 
                      genome = BSgenome.Hsapiens.UCSC.hg38)
# link peaks to genes
sc_obj_cd103trm17 <- LinkPeaks(
  object = sc_obj_cd103trm17,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = deg_cd103trm17
)
links_gr <- Links(sc_obj_cd103trm17[["peaks"]])
saveRDS(links_gr, file = "03_output/04_Diff/peak2gene/TI/LinkPeaks_links_CD103_TRM_Th17_DEG_II_vs_NN.rds")
