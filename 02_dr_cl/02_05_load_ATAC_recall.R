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

##### load ATAC data #####
meta_use_clinic_add <- readRDS("03_output/02_clean/meta_use_clinic_add.rds")
use_sample <- unique(meta_use_clinic_add$Batch)
scATAC_obj_raw <- load.seurat(glue("02_clean/{use_sample}_atac.rds"), 
                              type = "peaks",
                              merge = T)
scATAC_obj_raw@meta.data <- meta_use_clinic_add[colnames(scATAC_obj_raw), ]
saveRDS(scATAC_obj_raw, file = "03_output/02_clean/scATAC_obj_raw.rds")
##
sc_meta_raw <- readRDS("03_output/03_clustering/raw/sc_meta_raw.rds")
scATAC_obj <- subset(scATAC_obj_raw, 
                     cells = rownames(sc_meta_raw)[sc_meta_raw$ann_raw == "Immune"])
saveRDS(scATAC_obj, file = "03_output/02_clean/scATAC_obj_immune.rds")

##### peak recalling on all immune cells #####
scATAC_obj_recall <- peak.recall(scATAC_obj = scATAC_obj,
                                 ref_peaks = NULL,
                                 group_by = NULL,
                                 macs2_path = "/ix1/wchen/xiangyu/conda_env/macs2/bin/macs2",
                                 assay_use = "peaks",
                                 new_assay = "peaks",
                                 combine_peaks = F,
                                 tmp_path = "/ix3/wchen/xiy231/test/03_CD_DOGMA/tmp/",
                                 out_prefix = "03_output/03_clustering/recall_comb")
##
plan("multisession", workers = 10)
scATAC_act <- GeneActivity(scATAC_obj_recall,
                           assay = "peaks",
                           extend.upstream = 2000,
                           extend.downstream = 0,
                           process_n = 2000,
                           verbose = T)
plan("sequential")
saveRDS(scATAC_act, file = "03_output/03_clustering/scATAC_act_mat_immune.rds")

##### peak recalling group by major cell types (run after level2 annotation frozen) #####
scATAC_obj_recall <- readRDS("03_output/03_clustering/recall_comb_scATAC_obj.rds")
scATAC_obj_recall@meta.data <- readRDS("03_output/03_clustering/all/WNN_ADT_RNA/sc_meta_ann_level2_refine.rds")

plan("multisession", workers = 10)
scATAC_obj_recall <- peak.recall(scATAC_obj = scATAC_obj_recall,
                                 ref_peaks = NULL,
                                 group_by = "ann_level2_refine",
                                 macs2_path = "/ix1/wchen/xiangyu/conda_env/macs2/bin/macs2",
                                 assay_use = "peaks",
                                 new_assay = "peaks",
                                 combine_peaks = T,
                                 tmp_path = "/ix3/wchen/xiy231/test/03_CD_DOGMA/tmp/",
                                 out_prefix = "03_output/03_clustering/recall_comb_ann_level2_refine")
plan("sequential")
##
scATAC_obj_recall[["chromvar"]] <- NULL
scATAC_obj_recall <- anno.motif(scATAC_obj = scATAC_obj_recall,
                                assay_name = "peaks",
                                species_use = "human")
saveRDS(scATAC_obj_recall[["chromvar"]], file = "03_output/03_clustering/recall_comb_ann_level2_refine_scchromvar_assay.rds")
saveRDS(scATAC_obj_recall, file = "03_output/03_clustering/recall_comb_ann_level2_refine_addmotif.rds")

##### output consensus_regions.bed for SCENIC+ #####
recall_comb_ann_peaks <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_peaks.rds")
consensus_regions <- as.data.frame(recall_comb_ann_peaks)[,c("seqnames", "start", "end")]
scenic_path <- "03_output/07_SCENIC/"
system(glue("mkdir -p {scenic_path}"))
write.table(consensus_regions, file = glue("{scenic_path}consensus_regions.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)

##### annotation from ChIPseeker #####
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(glue)
library(dplyr)
library(stringr)
##
recall_comb_ann_peaks <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_peaks.rds")
win_size <- 3000
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peakAnno <- annotatePeak(recall_comb_ann_peaks, 
                         TxDb = txdb, 
                         tssRegion = c(-win_size, win_size), 
                         annoDb = "org.Hs.eg.db")
saveRDS(peakAnno, file = glue("03_output/03_clustering/peakAnno_recall_comb_ann_level2_refine_peaks_win{win_size}.rds"))
##
peakAnno_df <- peakAnno@anno %>% as.data.frame()
peakAnno_df$ann <- str_split_i(peakAnno_df$annotation, "\\(", 1) %>%
  gsub(" $", "", .)
peakAnno_df$term <- paste0(peakAnno_df$seqnames, "-", peakAnno_df$start, "-", peakAnno_df$end)
saveRDS(peakAnno_df, file = glue("03_output/03_clustering/peakAnno_df_recall_comb_ann_level2_refine_peaks_win{win_size}.rds"))

