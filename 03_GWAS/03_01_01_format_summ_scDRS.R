#
library(bigreadr)
library(dplyr)
library(stringr)
library(glue)
library(rtracklayer)
library(bigsnpr)
library(Signac)

#
ref_path <- "/ix1/wchen/xiangyu/Ref_data/"
summ_raw_path <- glue("{ref_path}GWAS/Raw/")
clean_path <- glue("{ref_path}GWAS/GRCh37/")

#
hm3_snp_list <- fread2(glue("{ref_path}ldsc/eur_w_ld_chr/w_hm3.snplist"))[, 1, drop = T]
hg19ToHg38_chain <- import.chain("/ix1/wchen/xiangyu/Ref_data/genome/homo_sapiens/hg19ToHg38.over.chain")
hg38ToHg19_chain <- import.chain("/ix1/wchen/xiangyu/Ref_data/genome/homo_sapiens/hg38ToHg19.over.chain")
comm_base <- c("T", "C", "G", "A")
ref_bim <- fread2("/ix1/wchen/xiangyu/Ref_data/genome/1000GP/EUR/merge.bim")[,-3]
colnames(ref_bim) <- c("chr", "rsid", "pos", "a0", "a1")

#########################
##### CD: 40266 GRCh37 #####
# Allele2 = effect_allele
pheno <- "CD"
n <- 40266
summ_CD <- fread2(glue("{summ_raw_path}cd_build37_40266_20161107.txt.gz"))
summ_CD <- cbind(summ_CD, 
                 (str_split(summ_CD$MarkerName, "_", simplify = T)[,1]) %>%
                   str_split(., "\\:", simplify = T) %>%
                   as.data.frame())
summ_CD <- summ_CD[,c("V1", "V2", "Allele1", "Allele2", "Effect", "StdErr", "P.value")]
summ_CD$V1 <- as.integer(summ_CD$V1)
summ_CD$V2 <- as.integer(summ_CD$V2)
colnames(summ_CD) <- c("chr", "pos", "a0", "a1", "beta", "se", "p")
summ_CD$a0 <- toupper(summ_CD$a0)
summ_CD$a1 <- toupper(summ_CD$a1)
summ_CD <- subset(summ_CD, a0 %in% comm_base & 
                    a1 %in% comm_base &
                    se > 0)
summ_CD$z <- summ_CD$beta / summ_CD$se
##
summ_info <- snp_match(sumstats = summ_CD, 
                       info_snp = ref_bim, 
                       strand_flip = TRUE)

summ_info <- summ_info[,c("chr", "pos", "rsid", "a0", "a1", "beta", "se", "z", "p")]
summ_info$n <- N
fwrite2(summ_info,
        file = glue("{clean_path}{pheno}.txt"),
        sep = "\t")
# ##
# # summ_info <- fread2(glue("{clean_path}{pheno}.txt"))
# gr <- StringToGRanges(paste0("chr",
#                              summ_info$chr, "-",
#                              summ_info$pos, "-",
#                              summ_info$pos + 1))
# gr38_df <- liftOver(gr, 
#                     hg19ToHg38_chain) %>%
#   unlist() %>%
#   GRangesToString(., sep = c("-", "-")) %>%
#   str_split(., "\\-", simplify = T) %>%
#   as.data.frame()
# gr38_df[, 1] <- gsub("chr", "", gr38_df[,1]) %>% as.integer()
# gr38_df[, 2] <- as.integer(gr38_df[, 2])
# summ_info_38 <- cbind(chr = gr38_df[,1],
#                       pos = gr38_df[,2],
#                       summ_info[,-c(1:2)])
# summ_info_38 <- summ_info_38[!is.na(summ_info_38$chr),]
# 
# fwrite2(summ_info_38,
#         file = glue("{ref_path}GWAS/GRCh38/{pheno}.txt"),
#         sep = "\t")
##
summ_info <- fread2(file = glue("{clean_path}{pheno}.txt"))
summ_ldsc <- subset(summ_info, 
                    select = c("rsid", "pos", "a1", "a0", "z", "n"))
colnames(summ_ldsc) <- c("SNP", "POS", "A1", "A2", "Z", "N")
fwrite2(summ_ldsc,
        file = glue("{clean_path}/{pheno}.ldsc"),
        sep = "\t")
system(glue("gzip -f {clean_path}/{pheno}.ldsc"))

##### IBD: 59957 GRCh37 #####
# Allele2 = effect_allele
pheno <- "IBD"
N <- 59957
summ_IBD <- fread2(glue("{summ_raw_path}ibd_build37_59957_20161107.txt.gz"))
summ_IBD <- cbind(summ_IBD, 
                  (str_split(summ_IBD$MarkerName, "_", simplify = T)[,1]) %>%
                    str_split(., "\\:", simplify = T) %>%
                    as.data.frame())
summ_IBD <- summ_IBD[,c("V1", "V2", "Allele1", "Allele2", "Effect", "StdErr", "P.value")]
summ_IBD$V1 <- as.integer(summ_IBD$V1)
summ_IBD$V2 <- as.integer(summ_IBD$V2)
colnames(summ_IBD) <- c("chr", "pos", "a0", "a1", "beta", "se", "p")
summ_IBD$a0 <- toupper(summ_IBD$a0)
summ_IBD$a1 <- toupper(summ_IBD$a1)
summ_IBD <- subset(summ_IBD, a0 %in% comm_base & 
                     a1 %in% comm_base &
                     se > 0)
summ_IBD$z <- summ_IBD$beta / summ_IBD$se
##
summ_info <- snp_match(sumstats = summ_IBD, 
                       info_snp = ref_bim, 
                       strand_flip = TRUE)

summ_info <- summ_info[,c("chr", "pos", "rsid", "a0", "a1", "beta", "se", "z", "p")]
summ_info$n <- N
fwrite2(summ_info,
        file = glue("{clean_path}{pheno}.txt"),
        sep = "\t")
##
summ_ldsc <- subset(summ_info, 
                    select = c("rsid", "pos", "a1", "a0", "z", "n"))
colnames(summ_ldsc) <- c("SNP", "POS", "A1", "A2", "Z", "N")
fwrite2(summ_ldsc,
        file = glue("{clean_path}/{pheno}.ldsc"),
        sep = "\t")
system(glue("gzip -f {clean_path}/{pheno}.ldsc"))

##### UC: 45975 GRCh37 #####
# Allele2 = effect_allele
pheno <- "UC"
N <- 45975
summ_UC <- fread2(glue("{summ_raw_path}uc_build37_45975_20161107.txt.gz"))
summ_UC <- cbind(summ_UC, 
                 (str_split(summ_UC$MarkerName, "_", simplify = T)[,1]) %>%
                   str_split(., "\\:", simplify = T) %>%
                   as.data.frame())
summ_UC <- summ_UC[,c("V1", "V2", "Allele1", "Allele2", "Effect", "StdErr", "P.value")]
summ_UC$V1 <- as.integer(summ_UC$V1)
summ_UC$V2 <- as.integer(summ_UC$V2)
colnames(summ_UC) <- c("chr", "pos", "a0", "a1", "beta", "se", "p")
summ_UC$a0 <- toupper(summ_UC$a0)
summ_UC$a1 <- toupper(summ_UC$a1)
summ_UC <- subset(summ_UC, a0 %in% comm_base & 
                    a1 %in% comm_base &
                    se > 0)
summ_UC$z <- summ_UC$beta / summ_UC$se
##
summ_info <- snp_match(sumstats = summ_UC, 
                       info_snp = ref_bim, 
                       strand_flip = TRUE)

summ_info <- summ_info[,c("chr", "pos", "rsid", "a0", "a1", "beta", "se", "z", "p")]
summ_info$n <- N
fwrite2(summ_info,
        file = glue("{clean_path}{pheno}.txt"),
        sep = "\t")
##
summ_ldsc <- subset(summ_info, 
                    select = c("rsid", "pos", "a1", "a0", "z", "n"))
colnames(summ_ldsc) <- c("SNP", "POS", "A1", "A2", "Z", "N")
fwrite2(summ_ldsc,
        file = glue("{clean_path}/{pheno}.ldsc"),
        sep = "\t")
system(glue("gzip -f {clean_path}/{pheno}.ldsc"))



