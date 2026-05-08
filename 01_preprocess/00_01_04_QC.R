#
library(bigreadr)
library(stringr)
library(optparse)
suppressMessages(library(glue))
suppressMessages(library(dplyr))

## Input parameters
args_list = list(
  make_option("--in_path_arc", type="character", default=NULL,
              help="INPUT: raw filepath of cellranger-arc", metavar="character"),
  make_option("--arc_fragfile", type="character", default=NULL,
              help="INPUT: alternative filepath of fragments", metavar="character"),
  make_option("--in_path_adt", type="character", default=NULL,
              help="INPUT: raw filepath of ADT data", metavar="character"),
  make_option("--in_path_hto", type="character", default=NULL,
              help="INPUT: raw filepath of HTO data", metavar="character"),
  make_option("--hto_valid_file", type="character", default=NULL,
              help="INPUT: file contain valid HTO", metavar="character"),
  make_option("--clean_path", type="character", default=NULL,
              help="INPUT: clean output path", metavar="character"),
  make_option("--doublets", type="character", default="NONE",
              help="INPUT: doublets file", metavar="character"),
  make_option("--threads", type="integer", default=1,
              help="INPUT: number of threads", metavar="integer"),
  make_option("--idx", type="character", default=NULL,
              help="INPUT: Smaple ID", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## set parameters
in_path_arc <- opt$in_path_arc
arc_fragfile <- opt$arc_fragfile
in_path_adt <- opt$in_path_adt
in_path_hto <- opt$in_path_hto
hto_valid_file <- opt$hto_valid_file
clean_path <- opt$clean_path
doublets <- opt$doublets
threads <- opt$threads
idx <- opt$idx
##
seed_use <- 20250528

# idx <- "Duerr_20230124_DOGMAseq-1"
# in_path_arc <- glue("/ix1/rduerr/shared/rduerr_wchen/raw/Cell_Ranger_output/cellranger-arc_count_2.0.2/Helmsley_2nd_6_months/{idx}_arc/")
# in_path_adt <- glue("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/01_raw/{idx}/ADT/featurecounts/")
# in_path_hto <- glue("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/01_raw/{idx}/HTO/featurecounts/")
# hto_valid_file <- glue("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/01_raw/{idx}/HTO/Feature_HTOs_valid.csv")
# arc_fragfile <- glue("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/01_raw/{idx}/modified_fragments.tsv.gz")
# clean_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/02_clean/"

ref_path <- "/ix1/wchen/xiangyu/Ref_data/"
code_path <- "/ix1/wchen/xiangyu/Ref_code/single_cell/"
source(glue("{code_path}00_FUNCTION/process/QC_FUN.R"))

## set files
arc_h5file <- glue("{in_path_arc}/filtered_feature_bc_matrix.h5")
arc_metafile <- glue("{in_path_arc}per_barcode_metrics.csv")
if (is.null(arc_fragfile)) {
  arc_fragfile <- glue("{in_path_arc}atac_fragments.tsv.gz")
}
grange_file <- glue("{ref_path}genome/grange_Anno/EnsDb.Hsapiens.v86.rds")

## qc on each data type
metadata_all <- read.csv(arc_metafile, header = T, row.names = 1)
metadata <- subset(metadata_all, 
                   is_cell == 1,
                   select = c("gex_barcode", "atac_barcode"))
if (!file.exists(doublets)) {
  run_dbl <- T
  metadata$pass_COM <- T
  message("Would run doublets detection as no results provided!")
} else {
  multiplet <- fread2(doublets)
  run_dbl <- F
  metadata$pass_COM <- !multiplet$multiplet_classification
  message("Use doublets detection results provided!")
  
}

##### Demultiplex and QC on HTO data #####
valid_hto <- read.csv(hto_valid_file, header = F)[, 1, drop = T]
valid_hto <- gsub("\\_", "-", valid_hto)
qc_out_hto <- qc.fun.HTO(hto_path = in_path_hto,
                         meta_barcode = metadata$gex_barcode,
                         max_ncounts = 5E3,
                         valid_hto = valid_hto,
                         idx = idx,
                         out_path = clean_path,
                         save_rds = T,
                         return_seurat = T)
message("Demultiplexing and QC done for HTO data!")

##### qc on ADT data #####
qc_out_adt <- qc.fun.ADT(adt_path = in_path_adt,
                         meta_barcode = metadata$gex_barcode,
                         min_features = 0,
                         max_ncounts = 2E4,
                         min_ncounts = 200,
                         max_percent_ctrl = 2,
                         max_ncounts_ctrl = Inf,
                         idx = idx,
                         out_path = clean_path,
                         save_rds = T,
                         return_seurat = T)
message("QC done for ADT data!")
##### qc on gene expression data #####
qc_out_gex <- qc.fun.RNA(h5_file = arc_h5file,
                         ARC = T,
                         run_dbl = F,
                         dbr_sd = 1,
                         min_features = 200,
                         min_ncounts = 1000,
                         max_ncounts = Inf,
                         max_mtpercent = 20,
                         idx = idx,
                         n_core = threads,
                         out_path = clean_path,
                         save_rds = T,
                         return_seurat = T)
message("QC done for GEX data!")
##### QC on ATAC data #####
qc_out_atac <- qc.fun.ATAC(h5_arc_file = arc_h5file,
                           frag_file = arc_fragfile,
                           grange_file = grange_file,
                           run_dbl = F,
                           dbr_sd = 1,
                           seed_db = seed_use,
                           min_features = 0,
                           min_ncounts = 1E3,
                           max_ncounts = Inf,
                           max_artef_rt = Inf,
                           min_tss = 4,
                           max_nucleosome = 4,
                           idx = idx,
                           n_core = threads,
                           out_path = clean_path,
                           save_rds = T,
                           return_seurat = T)
message("QC done for ATAC data!")
##### combine output ####
metadata$Batch_barcode <- paste0(idx, "_", metadata$gex_barcode)
qc_out_hto <- readRDS(glue("{clean_path}/{idx}_hto_tmp.rds"))
qc_out_adt <- readRDS(glue("{clean_path}/{idx}_adt_tmp.rds"))
qc_out_gex <- readRDS(glue("{clean_path}/{idx}_gex_tmp.rds"))
qc_out_atac <- readRDS(glue("{clean_path}/{idx}_atac_tmp.rds"))
sample_info <- fread2(glue("{clean_path}{idx}_sample_df.txt"))
# add qc matrix
metadata <- cbind(metadata,
                  Batch = idx,
                  qc_out_hto@meta.data[,c("nCount_HTO", "nFeature_HTO", "HTO_margin", "HTO_classification.global", "pass_HTO", "hash.ID")],
                  qc_out_adt@meta.data[,c("nCount_ADT", "nFeature_ADT", "percent.ctrl", "sum.ctrl", "pass_ADT")],
                  qc_out_gex@meta.data[,c("nFeature_RNA", "nCount_RNA", "percent.mt", "pass_RNA")],
                  qc_out_atac@meta.data[,c('nCount_peaks', "nFeature_peaks", 'nucleosome_signal', 'TSS.enrichment', 'blacklist_ratio', "pass_ATAC")])
metadata$pass <- metadata$pass_COM &
  metadata$pass_HTO & metadata$pass_ADT &
  metadata$pass_RNA & metadata$pass_ATAC
n_excl <- sum(!metadata$pass)
message(glue("A total of {n_excl} out of {nrow(metadata)} cells were excluded from this sample!"))
# add sample info
metadata <- cbind(metadata,
                  sample_info[match(metadata$hash.ID, sample_info$HTO_id),
                              c("ID", "Sample", "Sample_ID", "Sample_exp", "Sample_ID_exp",
                                "Condition", "sequence", "Section", "Date")])
## output
saveRDS(metadata, file = glue("{clean_path}/{idx}_meta.rds"))
saveRDS(subset(qc_out_hto, cells = metadata$Batch_barcode[which(metadata$pass)]),
        file = glue("{clean_path}/{idx}_hto.rds"))
saveRDS(subset(qc_out_adt, cells = metadata$Batch_barcode[which(metadata$pass)]),
        file = glue("{clean_path}/{idx}_adt.rds"))
saveRDS(subset(qc_out_gex, cells = metadata$Batch_barcode[which(metadata$pass)]),
        file = glue("{clean_path}/{idx}_gex.rds"))
saveRDS(subset(qc_out_atac, cells = metadata$Batch_barcode[which(metadata$pass)]),
        file = glue("{clean_path}/{idx}_atac.rds"))

## pre-evaluation for ADT denoising
## set files
arc_metafile <- glue("{in_path_arc}per_barcode_metrics.csv")
metadata_all <- read.csv(arc_metafile, header = T, row.names = 1)
metadata <- readRDS(glue("{clean_path}/{idx}_meta.rds"))
blank_barcode <- metadata_all$gex_barcode[which(metadata_all$is_cell == 0 &
                                                  metadata_all$excluded_reason == 0)]
#
denoise_pre_adt <- denoise.pre.ADT(adt_path = in_path_adt,
                                   clean_barcode = metadata$gex_barcode[which(metadata$pass)],
                                   blank_barcode = blank_barcode,
                                   idx = idx,
                                   out_path = clean_path)
write.table(unlist(denoise_pre_adt) %>% t,
            file = glue("{clean_path}/{idx}_blank_pick.txt"),
            sep = "\t", col.names = T, row.names = F, quote = F)
