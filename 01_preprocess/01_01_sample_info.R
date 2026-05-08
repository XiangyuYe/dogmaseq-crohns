library(dplyr)
library(glue)

fq_path <- "/ix1/rduerr/shared/rduerr_wchen/raw/Single_cell_fastq/"
cr_path <- "/ix1/rduerr/shared/rduerr_wchen/raw/Cell_Ranger_output/cellranger-arc_count_2.0.2/"
#
fq_batch <- list.files(fq_path, full.names = F, pattern = "Helmsley")
cr_batch <- list.files(cr_path, full.names = F, pattern = "Helmsley")
#
fq_df <- lapply(fq_batch, function(batchx){
  data.frame(id = list.files(glue("{fq_path}/{batchx}/"), full.names = F) %>%
               gsub("_fastq", "", .),
             batch_fq = gsub("_fastq", "", batchx))
})%>% Reduce("rbind", .) %>% as.data.frame()
#
cr_df <- lapply(cr_batch, function(batchx){
  data.frame(id = list.files(glue("{cr_path}/{batchx}/"), full.names = F) %>%
               gsub("_arc", "", .),
             batch_cr = batchx)
})%>% Reduce("rbind", .)  %>% as.data.frame()

sample_df <- merge(fq_df, cr_df, by = "id")
write.table(sample_df, file = "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/Helmsley_sample.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(sample_df, file = "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/Helmsley_sample.txt",
            col.names = F, row.names = F, quote = F, sep = "\t")

###
library(openxlsx)
library(dplyr)
library(stringr)
library(bigreadr)
library(glue)
##
proj_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
sample_id_df <- read.table(file = glue("{proj_path}/Helmsley_sample.txt"),
                         header = F, sep = "\t")
all_idx <- sample_id_df[,1]
experi_file <- list.files("/ix1/rduerr/shared/rduerr_wchen/Experiment_summary_feature_ref_cite_seq_index/",
                          pattern = "Single_Cell_Multimodal_Omics_Experiments", 
                          full.names = T)[1]
experi_info <- openxlsx::read.xlsx(experi_file, sheet = 1)
feature_ref <- read.csv("/ix1/rduerr/shared/rduerr_wchen/Experiment_summary_feature_ref_cite_seq_index/feature_ref_TotalSeq-A-Human-Universal-Cocktail-V1_HTOs.csv")
#
lapply(all_idx, function(idx){
  
  # idx <- "Duerr_20241010_DOGMAseq"
  sample_df <- experi_info[experi_info$Experiment_Name == idx, 
                           c("Experiment_Name", "Human_Subject_Identifier", "Human_Subject_Phenotype", "Sample.name", 
                             "Enriched_Cell_Population(s)", "Cell_Culture_Stimulation_Condition", "Barcode")]
  colnames(sample_df) <- c("ID", "Sample", "Group", "Sample_ID", "Sample_type", "Condition", "sequence")
  sample_df$HTO_id <- feature_ref$id[match(sample_df$sequence, feature_ref$sequence)] %>%
    gsub("\\_", "\\-", .)
  sample_df$Condition <- ifelse(grepl("N_N", sample_df$Sample_ID), "NN",
                                ifelse(grepl("N_U", sample_df$Sample_ID), "NU", "II"))
  # LC     RC     RM     SC     TC     # TI 
  sample_df$Section <- ifelse(grepl("_TI$", sample_df$Sample_ID), "TI", "Colon")
  sample_df$Date <- str_split_i(sample_df$ID, "_", 2)
  sample_df$Sample_exp <- paste0(sample_df$Sample, 
                                 "_",
                                 sample_df$Date)
  sample_df$Sample_ID_exp <- paste0(sample_df$Sample_ID, 
                                 "_",
                                 sample_df$Date)
  write.table(sample_df, file = glue("{proj_path}02_clean/{idx}_sample_df.txt"),
              sep = "\t", row.names = F, col.names = T, quote = F)
  Feature_HTOs <- sample_df[, c("HTO_id", "sequence")]
  Feature_HTOs$HTO_id <- gsub("\\-", "\\_", Feature_HTOs$HTO_id)
  system(glue("mkdir -p {proj_path}01_raw/{idx}/HTO/"))
  fwrite2(Feature_HTOs, file = glue("{proj_path}01_raw/{idx}/HTO/Feature_HTOs_valid.csv"),
          col.names = F, row.names = F, quote = F, sep = ",")
})


