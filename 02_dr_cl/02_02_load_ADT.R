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
meta_use_clinic_add <- readRDS("03_output/02_clean/meta_use_clinic_add.rds")
use_sample <- unique(meta_use_clinic_add$Batch)
scADT_obj_raw <- load.seurat(glue("02_clean/{use_sample}_adt_dsb.rds"),
                             type = "ADT",
                             merge = T)
scADT_obj_raw@meta.data <- meta_use_clinic_add[colnames(scADT_obj_raw), ]
saveRDS(scADT_obj_raw, file = "03_output/02_clean/scADT_obj_raw.rds")
##
sc_meta_raw <- readRDS("03_output/03_clustering/raw/sc_meta_raw.rds")
scADT_obj <- subset(scADT_obj_raw,
                    cells = rownames(sc_meta_raw)[sc_meta_raw$ann_raw == "Immune"])
saveRDS(scADT_obj, file = "03_output/02_clean/scADT_obj_immune.rds")

