### CoMM
rm(list=ls())
gc()
library(optparse)
library(bigreadr)

## Input parameters
args_list = list(
  make_option("--megama_out", type="character", default=NULL,
              help="INPUT: megama out file path", metavar="character"),
  make_option("--grch", type="numeric", default=NULL,
              help="INPUT: GRCh version", metavar="numeric"),
  make_option("--pheno", type="character", default=NULL,
              help="INPUT: phenotype", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

megama_out <- opt$megama_out
grch <- opt$grch
pheno <- opt$pheno

# megama_out <- ""
# grch <- 37
# pheno <- 
#
ref_file <- ifelse(grch == 38,
                   "/ix1/wchen/xiangyu/Ref_data/genome/homo_sapiens/NCBI38.gene.loc",
                   "/ix1/wchen/xiangyu/Ref_data/genome/homo_sapiens/NCBI37.3.gene.loc")
ref_gene <- fread2(ref_file)
magma_df <- fread2(megama_out)
magma_df$GENE <- ref_gene$V6[match(magma_df$GENE, ref_gene$V1)]
zscore_df <- magma_df[,c("GENE", "ZSTAT")]
colnames(zscore_df)[2] <- pheno
fwrite2(zscore_df, file = paste0(dirname(megama_out), "/zscore_", pheno, ".tsv"), sep = "\t")

