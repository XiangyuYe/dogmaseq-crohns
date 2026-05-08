library(bigreadr)
library(dplyr)
library(glue)

proj_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(proj_path)
cov_use <- "cov"

out_path <- "03_output/01_scDRS/03_integrate/"
scdrs_score <- fread2(glue("{out_path}{cov_use}/CD.full_score.gz")) %>% 
  tibble::column_to_rownames(., "V1")
sc_meta_CD4T <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta.rds")
sc_meta_CD8T <- readRDS("03_output/03_clustering/CD8T/WNN_ADT_RNA/sc_meta.rds")
sc_meta_cd103trm <- readRDS("03_output/03_clustering/CD4_CD103_TRM/WNN_RNA_ATAC/sc_meta.rds")

## CD4
scdrs_score_CD4T <- scdrs_score[rownames(sc_meta_CD4T),]
system(glue("mkdir -p {out_path}CD4T/all_sub/{cov_use}/"))
fwrite2(scdrs_score_CD4T, file = glue("{out_path}CD4T/all_sub/{cov_use}/CD.full_score"), 
        row.names = T, sep = "\t")
system(glue("gzip {out_path}CD4T/all_sub/{cov_use}/CD.full_score"))

## CD8
scdrs_score_CD8T <- scdrs_score[rownames(sc_meta_CD8T),]
system(glue("mkdir -p {out_path}CD8T/all_sub/{cov_use}/"))
fwrite2(scdrs_score_CD8T, file = glue("{out_path}CD8T/all_sub/{cov_use}/CD.full_score"), 
        row.names = T, sep = "\t")
system(glue("gzip {out_path}CD8T/all_sub/{cov_use}/CD.full_score"))


## CD4_CD103_TRM
scdrs_score_CD103trm <- scdrs_score[rownames(sc_meta_cd103trm),]
system(glue("mkdir -p {out_path}CD4_CD103_TRM/all_sub/{cov_use}/"))
fwrite2(scdrs_score_CD103trm, file = glue("{out_path}CD4_CD103_TRM/all_sub/{cov_use}/CD.full_score"), 
        row.names = T, sep = "\t")
system(glue("gzip {out_path}CD4_CD103_TRM/all_sub/{cov_use}/CD.full_score"))

