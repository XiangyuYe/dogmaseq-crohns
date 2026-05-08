## module load r/4.5.0
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(glue)
library(bigreadr)
library(cols4all)

## set path
proj_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(proj_path)
ref_path <- "/ix1/wchen/xiangyu/Ref_data/NicheNet/"
data_path <- "03_output/05_Interaction/multinichenetr/"

## load in LR-target reference
organism = "human"
ligand_target_matrix <- readRDS(glue("{ref_path}Human/ligand_target_matrix.rds"))

## set parameters
# set cell type of interest
celltype_id <- "ann_level5_final"
message(ct)
de_method <- "DESeq2"
message(de_method)

# 0.2 set cell type of interest
ct_rs <- c("CD4_CD103_TRM_Th17", "CD4T")
ct_s <- ct_rs[1]
ct_r <- ct_rs[2]

message(glue("{ct_s} to {ct_r}."))
celltype_id_rs <- c("ann_level5_final", "ann_level4_final")

#
ct_level <- readRDS("03_output/03_clustering/ct_levels_refine.rds")
all_ct_s <- ct_level[[celltype_id_rs[1]]]
senders_oi <- all_ct_s[grep("CD4_CD103_TRM_Th17$|CD4_CD103_TRM_Th17_", all_ct_s)]
all_ct_r <- ct_level[[celltype_id_rs[2]]]
receivers_oi <- all_ct_r[grep("Treg", all_ct_r)]

# plot parameters
source("code/FUNCTION/multinichenet/make_sample_lr_prod_activity_plots_modify.R")
top_lr_circle <- 50 # show top 50 LR pairs in circle plot
top_lr_lpa <- 30 # show top 50 LR pairs in circle plot
cor_thresh <- 0.5
p_thresh <- 0.05

cd4cd103trm_cols <- c(
  "#FF7B36", "#FE0000",
  "#FDD501", "#DDCB82",
  "#4691CC", "#A6D6FA",
  "#082C8C", "#112E5A", 
  "#EC4B88"
)
cd4_cols <- c("#00897B", "#06C3F2", "#016FAD", "#004F78", "#00A0E9",
              "#FFEB3B", "#FFB300", "#F57C33", "#FF5722", "#D32F2F",  "#A14141",  "#5D4037",
              "#F8AAC0", "#D81B60", "#F06292")
all_cor <- c(cd4cd103trm_cols, cd4_cols[13:15])
names(all_cor) <- c(senders_oi, receivers_oi)

# for (type in c("TI", "Colon")) {
lapply(c("TI", "Colon"), function(type){
  
  # type <- "Colon"
  message(type)
  
  # for (contrast_name in c("II_vs_NN", "NU_vs_NN")) {
  lapply(c("II_vs_NN", "NU_vs_NN"), function(contrast_name){
    
    # set contrast
    # contrast_name <- c("II_vs_NN", "NU_vs_NN")[2]
    conditions_keep <- strsplit(contrast_name, "_vs_")[[1]]
    contrasts_oi <- c(paste(conditions_keep, collapse = "-"),
                      paste(rev(conditions_keep), collapse = "-"))
    contrast_tbl <- tibble(contrast = contrasts_oi, 
                           group = conditions_keep)
    ##### 0. load data #####
    out_path <- glue("{data_path}/{type}/{ct_s}/{de_method}/{celltype_id}/")
    message(glue("Loading data: {out_path}multinichenet_output_{contrast_name}_{ct_r}.rds"))
    multinichenet_output <- readRDS(glue("{out_path}/multinichenet_output_{contrast_name}_{ct_r}.rds"))
    
    ##### 1. ChordDiagram circos plots #####
    # title for circos plot
    all_g <- c("NN", "NU", "II")
    if (type == "TI") {
      title_all_g <- glue("{all_g} in Ileum")
    } else {
      title_all_g <- glue("{all_g} in Colon")
    }
    names(title_all_g) <- all_g
    ## 1.1 from B to CD4T
    prior_tbl_oi_inner <- get.prior.tbl(multinichenet_output = multinichenet_output,
                                        n_top_lr = top_lr_circle,
                                        senders_oi = senders_oi,
                                        receivers_oi = receivers_oi)
    # plot
    pdf(glue("{out_path}circos_{ct_s}_to_{ct_r}_{contrast_name}_top{top_lr_circle}.pdf"),
        height = 6, width = 6)
    make_circos_group_comparison_modify(prioritized_tbl_oi = prior_tbl_oi_inner, 
                                        colors_sender = all_cor[prior_tbl_oi_inner$sender %>% unique()], 
                                        colors_receiver = all_cor[prior_tbl_oi_inner$receiver %>% unique()],
                                        title_goi = title_all_g)
    dev.off()
    ##### 2. filter LR-target df #####
    lr_target_cor_f <- lr.target.filter(multinichenet_output = multinichenet_output,
                                        contrast_tbl = contrast_tbl,
                                        cor_thresh = cor_thresh,
                                        p_thresh = p_thresh)
    ##### 3. plot from each CD4CD103TRM17 to CD4Treg #####
    ## 3.1 up
    lapply(senders_oi, function(senders_oix){
      
      senders_oix <- senders_oi[1]
      # format plot
      lrt_plot_upx <- mnn.lrt.plot(multinichenet_output = multinichenet_output,
                                   contrast_tbl = contrast_tbl,
                                   lr_target_cor_f = lr_target_cor_f,
                                   lr_prod_activity = T,
                                   ligand_target = T,
                                   lr_target_cor = T,
                                   n_top_lpa = top_lr_lpa,
                                   receivers_oi = receivers_oi,
                                   senders_oi = senders_oix,
                                   groups_oi = conditions_keep[1])
      # output plot
      plt_suffix_3_1 <- glue("{senders_oix}_to_{ct_r}_{contrast_name}_top{top_lr_lpa}_up")
      mnn.lrt.plot.out(lr_prod_plot_list = lrt_plot_upx$lr_prod_plot_list,
                       ligand_target_plot_list = lrt_plot_upx$ligand_target_plot_list,
                       lr_target_cor_plot_list = lrt_plot_upx$lr_target_cor_plot_list,
                       out_path = out_path,
                       plt_suffix = plt_suffix_3_1)
      #
      return(plt_suffix_3_1)
      
    })
    ## 3.2 down
    lapply(senders_oi, function(senders_oix){
      
      # senders_oix <- senders_oi[1]
      # format plot
      lrt_plot_downx <- mnn.lrt.plot(multinichenet_output = multinichenet_output,
                                     contrast_tbl = contrast_tbl,
                                     lr_target_cor_f = lr_target_cor_f,
                                     lr_prod_activity = T,
                                     ligand_target = T,
                                     lr_target_cor = T,
                                     n_top_lpa = top_lr_lpa,
                                     receivers_oi = receivers_oi,
                                     senders_oi = senders_oix,
                                     groups_oi = conditions_keep[2])
      # output plot
      plt_suffix_3_2 <- glue("{senders_oix}_to_{ct_r}_{contrast_name}_top{top_lr_lpa}_down")
      mnn.lrt.plot.out(lr_prod_plot_list = lrt_plot_downx$lr_prod_plot_list,
                       ligand_target_plot_list = lrt_plot_downx$ligand_target_plot_list,
                       lr_target_cor_plot_list = lrt_plot_downx$lr_target_cor_plot_list,
                       out_path = out_path,
                       plt_suffix = plt_suffix_3_2)
      #
      return(plt_suffix_3_2)
      
    })
    
    ##### 4. plot from all CD4CD103TRM17 to each CD4Treg #####
    ## 4.1 up
    lapply(receivers_oi, function(receivers_oix){
      
      # receivers_oix <- receivers_oi[1]
      # format plot
      lrt_plot_upx <- mnn.lrt.plot(multinichenet_output = multinichenet_output,
                                   contrast_tbl = contrast_tbl,
                                   lr_target_cor_f = lr_target_cor_f,
                                   lr_prod_activity = T,
                                   ligand_target = T,
                                   lr_target_cor = T,
                                   n_top_lpa = top_lr_lpa,
                                   receivers_oi = receivers_oix,
                                   senders_oi = senders_oi,
                                   groups_oi = conditions_keep[1])
      # output plot
      plt_suffix_4_1 <- glue("{ct_s}_to_{receivers_oix}_{contrast_name}_top{top_lr_lpa}_up")
      mnn.lrt.plot.out(lr_prod_plot_list = lrt_plot_upx$lr_prod_plot_list,
                       ligand_target_plot_list = lrt_plot_upx$ligand_target_plot_list,
                       lr_target_cor_plot_list = lrt_plot_upx$lr_target_cor_plot_list,
                       out_path = out_path,
                       plt_suffix = plt_suffix_4_1)
      #
      return(plt_suffix_4_1)
      
    })
    
    ## 4.2 down
    lapply(receivers_oi, function(receivers_oix){
      
      # receivers_oix <- receivers_oi[1]
      lrt_plot_downx <- mnn.lrt.plot(multinichenet_output = multinichenet_output,
                                     contrast_tbl = contrast_tbl,
                                     lr_target_cor_f = lr_target_cor_f,
                                     lr_prod_activity = T,
                                     ligand_target = T,
                                     lr_target_cor = T,
                                     n_top_lpa = top_lr_lpa,
                                     receivers_oi = receivers_oix,
                                     senders_oi = senders_oi,
                                     groups_oi = conditions_keep[2])
      # output plot
      plt_suffix_4_2 <- glue("{ct_s}_to_{receivers_oix}_{contrast_name}_top{top_lr_lpa}_down")
      mnn.lrt.plot.out(lr_prod_plot_list = lrt_plot_downx$lr_prod_plot_list,
                       ligand_target_plot_list = lrt_plot_downx$ligand_target_plot_list,
                       lr_target_cor_plot_list = lrt_plot_downx$lr_target_cor_plot_list,
                       out_path = out_path,
                       plt_suffix = plt_suffix_4_2) 
      #
      return(plt_suffix_4_2)
      
    })
    
  })
  #
  return(0)  
  
})
