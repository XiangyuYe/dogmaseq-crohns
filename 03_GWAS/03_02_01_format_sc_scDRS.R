##
library(Seurat)
library(dplyr)
library(bigreadr)
library(stringr)
library(SeuratDisk)
library(glue)
library(reticulate)
library(DropletUtils)
library(Matrix)

## set path and parameters
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
scRDS_path <- "03_output/01_scDRS/00_scRNA/"
cov_use = c("Age", "Gender", 
            "use_Aminosalic", "use_Immunomodu", "use_MTX",
            "use_AntiTNF", "use_AntiIL23", "use_AntiIntegrin")
block_col = "Sample"
sample_col = "Sample_ID_exp"
ann_col <- c("ann_level2_refine", "ann_level3_final", "ann_level4_final")
group_col <- "Condition"
subgroup_col <- "Section"
##
sc_meta_all <- readRDS("03_output/03_clustering/sc_meta_ann_comb_all1205.rds")
clinic_df <- readRDS(file = "03_output/02_clean/clinic_df.rds")
sc_meta <- cbind(sc_meta_all[, c("nFeature_RNA", block_col, sample_col, 
                                 ann_col, group_col, subgroup_col)],
                 clinic_df[match(sc_meta_all[[sample_col]], clinic_df[[sample_col]]), cov_use])
sc_meta$sample_ct_level2 <- paste0(sc_meta$Section, ":", sc_meta$ann_level2_refine)
sc_meta$sample_ct_level3 <- paste0(sc_meta$Section, ":", sc_meta$ann_level3_final)
sc_meta$sample_ct_level4 <- paste0(sc_meta$Section, ":", sc_meta$ann_level4_final)

## format h5ad and cov for scDRS
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
scRNA_count <- scRNA_obj[["RNA"]]$counts %>%
  as.matrix() %>% as(., "CsparseMatrix")
system(glue("mkdir -p {scRDS_path}"))
if (file.exists(glue("{scRDS_path}/all/"))) {
  system(glue("rm -rf {scRDS_path}/all"))
}
write10xCounts(
  glue("{scRDS_path}/all/"),
  scRNA_count,
  gene.type = "Gene Expression"
)
##
fwrite2(sc_meta, 
        file = glue("{scRDS_path}all/metadata.tsv"), 
        row.names = T)
# format cov
cov_df <- cbind(data.frame(index = rownames(sc_meta), const = 1),
                sc_meta[,c("nFeature_RNA", cov_use)],
                psych::dummy.code(sc_meta$Sample),
                psych::dummy.code(sc_meta$Condition),
                TI = ifelse(sc_meta$Section == "TI", 1, 0))
fwrite2(cov_df, file = glue("{scRDS_path}all/cov.tsv"), sep = "\t")

