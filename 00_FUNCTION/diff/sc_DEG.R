# Libraries
library(MAST)
library(SingleCellExperiment)
library(ggpubr)
library(qs)
library(optparse)
library(data.table)
library(dittoSeq)
library(EnhancedVolcano)

###### MAST-LMM ######
mast.lmm <- function(seurat_obj = NULL,
                     do_norm = T,
                     group_col = NULL,
                     cov_col = NULL,
                     block_col = NULL,
                     group_levels = NULL,
                     use_assay = "RNA",
                     do_para = T
){
  ## trans to SCA object
  DefaultAssay(seurat_obj) <- use_assay
  if (do_norm) {
    seurat_obj <- NormalizeData(seurat_obj)
  }
  sca <- as.SingleCellExperiment(seurat_obj,
                                 assay = use_assay) %>% SceToSingleCellAssay
  g1 <- group_levels[2] # Case
  g2 <- group_levels[1] # Control
  ## set regresion formula
  if (!is.null(block_col)) {
    rd_term <- glue("(1|{block_col})")
  } else {
    rd_term <- NULL
  }
  formula_zlm <- paste0("~", paste(c(group_col, cov_col, rd_term), collapse = " + ")) %>% as.formula()
  message(glue("Fit: {formula_zlm}"))
  ## Run MAST
  zlmCond <- zlm(formula_zlm, 
                 sca, 
                 method = 'glmer', 
                 ebayes = F, 
                 fitArgsD = list(nAGQ = 0),
                 parallel = do_para)
  # saveRDS(zlmCond, file = glue("03_output/04_Diff/DEG/MAST/zlmCond_{g1}_vs_{g2}_MAST.rds"))
  message("Done zlm")
  ## Run likelihood ratio test for the condition coefficient.
  targ_contra <- paste0(group_col, g1)
  res_sum <- summary(zlmCond, 
                     doLRT = targ_contra,
                     parallel = do_para)$datatable %>% as.data.frame()
  saveRDS(res_sum, file = glue("03_output/04_Diff/DEG/MAST/res_sum_{g1}_vs_{g2}_MAST.rds"))
  ##
  res <- merge(res_sum[res_sum$contrast == targ_contra & res_sum$component == "H",
                       c("primerid", "Pr(>Chisq)")], #hurdle P values
               res_sum[res_sum$contrast == targ_contra & res_sum$component == "logFC", 
                       c("primerid", "coef", "z")], 
               by='primerid') #logFC coefficients
  res$padj <- p.adjust(res$`Pr(>Chisq)`, "BH")
  res$log2FC <- log2(exp(res$coef))
  
  res_df <- res[,c("primerid", "log2FC", "z", "Pr(>Chisq)", "padj")]
  colnames(res_df) <- c("Term", "log2FC", "stat", "pvalue", "padj")
  res_df <- res_df[!is.na(res_df$log2FC) & !is.nan(res_df$log2FC),]
  return(res_df)
  
}

##
sc.deg <- function(seurat_obj = NULL,
                   use_assay = NULL,
                   ident_col = NULL,
                   group_col = NULL,
                   cov_col = NULL,
                   sample_col = NULL,
                   block_col = NULL,
                   do_norm = T,
                   min_percent = 0.01,
                   paired_only = F,
                   test_use = "MAST-lmm",
                   do_para = T){
  ##
  if (any(!c(ident_col, group_col, cov_col, block_col, sample_col) %in% 
          colnames(seurat_obj@meta.data))) {
    stop("Not all required variables in the metadata of Seurat object!")
  }
  ##
  if (is.null(use_assay)) {
    use_assay <- DefaultAssay(seurat_obj)
  }
  # set cell type levels
  if (!is.null(levels(seurat_obj@meta.data[[ident_col]]))) {
    all_ct <- levels(seurat_obj@meta.data[[ident_col]])
  } else {
    all_ct <- unique(seurat_obj@meta.data[[ident_col]]) %>%sort()
  }
  #
  if (!is.null(levels(seurat_obj@meta.data[[group_col]]))) {
    group_levels <- levels(seurat_obj@meta.data[[group_col]])
  } else {
    group_levels <- unique(seurat_obj@meta.data[[group_col]]) %>%sort()
  }
  #
  if (!grepl("MAST", test_use)) {
    stop("Please use MAST for DEG identification!")
  }
  ## DEG by clusters
  deg_by_ct <- lapply(all_ct, function(ctx){
    print(ctx)
    ## set subset cells
    cell_use <- colnames(seurat_obj)[seurat_obj@meta.data[[ident_col]] == ctx]
    seurat_obj_ct <- subset(seurat_obj, cells = cell_use)
    if (paired_only) {
      ## keep paired samples only
      message("Only retained paired samples!")
      bulk_meta <- seurat_obj_ct@meta.data[!duplicated(seurat_obj_ct@meta.data[[sample_col]]), 
                                          c(sample_col, block_col, group_col, cov_col)]
      sample_retain <- bulk_meta[[sample_col]][duplicated(bulk_meta[[block_col]])]
      cell_use <- colnames(seurat_obj_ct)[seurat_obj_ct@meta.data[[sample_col]] %in% sample_retain]
      seurat_obj_ct <- subset(seurat_obj_ct, cells = cell_use)
      
    } else {
      message("Use all samples!")
    }
    ## filtering genes
    prop_exp_group <- lapply(unique(seurat_obj_ct@meta.data[[group_col]]), function(gx){
      
      cell_gx <- colnames(seurat_obj_ct)[seurat_obj_ct@meta.data[[group_col]] == gx]
      seurat_objx <- subset(seurat_obj_ct, cells = cell_gx)
      prop_exp <- rowSums(seurat_objx[[use_assay]]$counts > 0) / ncol(seurat_objx)
      
    }) %>% Reduce("cbind", .) %>% as.data.frame()
    ft_use <- Features(seurat_obj_ct)[rowMeans(prop_exp_group > min_percent) == 1]
    seurat_obj_ct[[use_assay]] <- subset(seurat_obj_ct[[use_assay]], features = ft_use)
    #
    if (min(table(seurat_obj_ct@meta.data[[group_col]])) > 3 &
        nrow(seurat_obj_ct) > 1) {
      ## Run DEG analysis
      res <- switch(test_use,
                    "MAST-lmm" = {
                      mast.lmm(seurat_obj = seurat_obj_ct,
                               group_col = group_col,
                               cov_col = cov_col,
                               block_col = block_col,
                               group_levels = group_levels,
                               do_norm = do_norm,
                               use_assay = use_assay,
                               do_para = do_para)
                    })
      res$cluster <- ctx
    } else {
      res <- data.frame(Term = NA, 
                        log2FC = NA, 
                        stat = NA, 
                        pvalue = NA, 
                        padj = NA,
                        cluster = ctx)
    }
    
    return(res)
    
  })
  # deg_by_ct <- lapply(all_ct, function(ctx){
  #   readRDS(glue("03_output/04_Diff/DAR/small_intestine/pseudo_bulk_N_U_vs_N_N_{ctx}.rds"))
  # })
  names(deg_by_ct) <- all_ct
  return(deg_by_ct)
}


