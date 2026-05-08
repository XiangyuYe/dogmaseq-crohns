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
celltype_id <- "ann_level4_final"
ct = "B"
message(ct)
de_method <- "DESeq2"
message(de_method)

# set receiver and sender
ct_level <- readRDS("03_output/03_clustering/ct_levels_refine.rds")
all_ct <- ct_level[[celltype_id]]
receivers_oi <- all_ct[grep("^CD4", all_ct)]
senders_oi <- all_ct[grep("^B_|^ABC|^GC|^Plasma_", all_ct)]
receivers_oi2 <- senders_oi
senders_oi2 <- receivers_oi
ct1 <- "B"
ct2 <- "CD4T"

# plot parameters
source("code/FUNCTION/multinichenet/make_sample_lr_prod_activity_plots_modify.R")
top_lr_circle <- 50 # show top 50 LR pairs in circle plot
top_lr_lpa <- 30 # show top 50 LR pairs in circle plot
cor_thresh <- 0.5
p_thresh <- 0.05

all_cor <- c("#098476", "#06C3F2", "#016FAD", "#004F78", "#002D6D",
             "#FFFF00", "#FED037", "#E65731", "#F57C33", "#904A30",
             "#C02020", "#E61737", "#EC6E73", "#F8AAC0")
names(all_cor) <- receivers_oi


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
    out_path <- glue("{data_path}/{type}/{ct}/{de_method}/{celltype_id}/")
    message(glue("Loading data: {out_path}multinichenet_output_{contrast_name}.rds"))
    multinichenet_output <- readRDS(glue("{out_path}/multinichenet_output_{contrast_name}.rds"))
    
    ##### 1. ChordDiagram circos plots #####
    ## 1.1 from B to CD4T
    prior_tbl_oi_inner <- get.prior.tbl(multinichenet_output = multinichenet_output,
                                        n_top_lr = top_lr_circle,
                                        senders_oi = senders_oi,
                                        receivers_oi = receivers_oi)
    # plot
    pdf(glue("{out_path}circos_{ct1}_to_{ct2}_{contrast_name}_top{top_lr_circle}.pdf"),
        height = 6, width = 6)
    make_circos_group_comparison_modify(prioritized_tbl_oi = prior_tbl_oi_inner, 
                                        colors_sender = all_cor[prior_tbl_oi_inner$sender %>% unique()], 
                                        colors_receiver = all_cor[prior_tbl_oi_inner$receiver %>% unique()])
    dev.off()
    ## 1.2 from CD4T to B
    prior_tbl_oi_outer <- get.prior.tbl(multinichenet_output = multinichenet_output,
                                        n_top_lr = 50,
                                        senders_oi = senders_oi2,
                                        receivers_oi = receivers_oi2)
    #
    pdf(glue("{out_path}circos_{ct2}_to_{ct1}_{contrast_name}_top{top_lr_circle}.pdf"),
        height = 6, width = 6)
    make_circos_group_comparison_modify(prioritized_tbl_oi = prior_tbl_oi_outer,
                                        colors_sender = all_cor[prior_tbl_oi_outer$sender %>% unique()],
                                        colors_receiver = all_cor[prior_tbl_oi_outer$receiver %>% unique()])
    dev.off()
    ##### 2. filter LR-target df #####
    lr_target_cor_f <- lr.target.filter(multinichenet_output = multinichenet_output,
                                        contrast_tbl = contrast_tbl,
                                        cor_thresh = cor_thresh,
                                        p_thresh = p_thresh)
    ##### 3. plot from each B to CD4 #####
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
      plt_suffix_3_1 <- glue("{senders_oix}_to_{ct2}_{contrast_name}_top{top_lr_lpa}_up")
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
      plt_suffix_3_2 <- glue("{senders_oix}_to_{ct2}_{contrast_name}_top{top_lr_lpa}_down")
      mnn.lrt.plot.out(lr_prod_plot_list = lrt_plot_downx$lr_prod_plot_list,
                       ligand_target_plot_list = lrt_plot_downx$ligand_target_plot_list,
                       lr_target_cor_plot_list = lrt_plot_downx$lr_target_cor_plot_list,
                       out_path = out_path,
                       plt_suffix = plt_suffix_3_2)
      #
      return(plt_suffix_3_2)
      
    })
    
    ##### 4. plot from all B to each CD4 #####
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
      plt_suffix_4_1 <- glue("{ct1}_to_{receivers_oix}_{contrast_name}_top{top_lr_lpa}_up")
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
      plt_suffix_4_2 <- glue("{ct1}_to_{receivers_oix}_{contrast_name}_top{top_lr_lpa}_down")
      mnn.lrt.plot.out(lr_prod_plot_list = lrt_plot_downx$lr_prod_plot_list,
                       ligand_target_plot_list = lrt_plot_downx$ligand_target_plot_list,
                       lr_target_cor_plot_list = lrt_plot_downx$lr_target_cor_plot_list,
                       out_path = out_path,
                       plt_suffix = plt_suffix_4_2) 
      #
      return(plt_suffix_4_2)
      
    })
    
    ##### 5. plot from all CD4 to each B #####
    ## 5.1 up
    lapply(receivers_oi2, function(receivers_oix){
      
      # receivers_oix <- receivers_oi2[1]
      # format plot
      lrt_plot_upx <- mnn.lrt.plot(multinichenet_output = multinichenet_output,
                                   contrast_tbl = contrast_tbl,
                                   lr_target_cor_f = lr_target_cor_f,
                                   lr_prod_activity = T,
                                   ligand_target = T,
                                   lr_target_cor = T,
                                   n_top_lpa = top_lr_lpa,
                                   receivers_oi = receivers_oix,
                                   senders_oi = senders_oi2,
                                   groups_oi = conditions_keep[1])
      # output plot
      plt_suffix_5_1 <- glue("{ct2}_to_{receivers_oix}_{contrast_name}_top{top_lr_lpa}_up")
      mnn.lrt.plot.out(lr_prod_plot_list = lrt_plot_upx$lr_prod_plot_list,
                       ligand_target_plot_list = lrt_plot_upx$ligand_target_plot_list,
                       lr_target_cor_plot_list = lrt_plot_upx$lr_target_cor_plot_list,
                       out_path = out_path,
                       plt_suffix = plt_suffix_5_1)
      #
      return(plt_suffix_5_1)
      
    })
    
    ## 5.2 down
    lapply(receivers_oi2, function(receivers_oix){
      
      # receivers_oix <- receivers_oi2[1]
      lrt_plot_downx <- mnn.lrt.plot(multinichenet_output = multinichenet_output,
                                     contrast_tbl = contrast_tbl,
                                     lr_target_cor_f = lr_target_cor_f,
                                     lr_prod_activity = T,
                                     ligand_target = T,
                                     lr_target_cor = T,
                                     n_top_lpa = top_lr_lpa,
                                     receivers_oi = receivers_oix,
                                     senders_oi = senders_oi2,
                                     groups_oi = conditions_keep[2])
      # output plot
      plt_suffix_5_2 <- glue("{ct2}_to_{receivers_oix}_{contrast_name}_top{top_lr_lpa}_down")
      mnn.lrt.plot.out(lr_prod_plot_list = lrt_plot_downx$lr_prod_plot_list,
                       ligand_target_plot_list = lrt_plot_downx$ligand_target_plot_list,
                       lr_target_cor_plot_list = lrt_plot_downx$lr_target_cor_plot_list,
                       out_path = out_path,
                       plt_suffix = plt_suffix_5_2)
      #
      return(plt_suffix_5_2)
      
    })
    #
    return(0)  
  })
  #
  return(0)  
  
})
