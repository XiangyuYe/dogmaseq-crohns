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
  make_option("--in_path_adt", type="character", default=NULL,
              help="INPUT: raw filepath of ADT data", metavar="character"),
  make_option("--clean_path", type="character", default=NULL,
              help="INPUT: clean output path", metavar="character"),
  make_option("--low_prot", type="numeric", default=NULL,
              help="INPUT: low CI for prot", metavar="numeric"),
  make_option("--hi_prot", type="numeric", default=NULL,
              help="INPUT: high CI for prot", metavar="numeric"),
  make_option("--idx", type="character", default=NULL,
              help="INPUT: Smaple ID", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

## set parameters
in_path_arc <- opt$in_path_arc
in_path_adt <- opt$in_path_adt
clean_path <- opt$clean_path
low_prot <- opt$low_prot
hi_prot <- opt$hi_prot
idx <- opt$idx
##
seed_use <- 20250528

# idx <- "Duerr_20230124_DOGMAseq-1"
# in_path_adt <- glue("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/01_raw/{idx}/ADT/featurecounts/")
# clean_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/02_clean/"

code_path <- "/ix1/wchen/xiangyu/Ref_code/single_cell/"
source(glue("{code_path}FUNCTION/process/QC_FUN.R"))

## set files
arc_metafile <- glue("{in_path_arc}per_barcode_metrics.csv")
metadata_all <- read.csv(arc_metafile, header = T, row.names = 1)
metadata <- readRDS(glue("{clean_path}/{idx}_meta.rds"))
## denoise on QCed ADT data
blank_barcode <- metadata_all$gex_barcode[which(metadata_all$is_cell == 0 &
                                                  metadata_all$excluded_reason == 0)]
denoise_pre_adt <- denoise.fun.ADT(adt_path = in_path_adt,
                                   clean_barcode = metadata$gex_barcode[which(metadata$pass)],
                                   blank_barcode = blank_barcode,
                                   idx = idx,
                                   low_prot = low_prot, 
                                   hi_prot = hi_prot,
                                   out_path = clean_path,
                                   seed_use = seed_use,
                                   save_rds = T)
message("Denoise done for ADT data!")

## add to 
qc_out_adt <- readRDS(glue("{clean_path}/{idx}_adt.rds"))
qc_out_adt[["ADT"]]$data <- readRDS(glue("{clean_path}/{idx}_adt_dsbmat.rds"))
saveRDS(qc_out_adt, file = glue("{clean_path}/{idx}_adt_dsb.rds"))
