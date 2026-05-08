##
library(dplyr)
library(Seurat)


###### pseudoBulk ######
pseudo.pair.diff <- function(seurat_obj = NULL,
                             assay_use = NULL,
                             feature_sel = NULL,
                             ident_col = NULL,
                             group_col = NULL,
                             sample_col = NULL,
                             block_col = NULL,
                             test_use = "wilcox"){
  ##
  if (any(!c(ident_col, group_col, block_col, sample_col) %in% 
          colnames(seurat_obj@meta.data))) {
    stop("Not all required variables in the metadata of Seurat object!")
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
  if (is.null(assay_use)) {
    assay_use <- DefaultAssay(seurat_obj)
    message("No assay provide. Use default assay!")
  }
  if (!assay_use %in% names(seurat_obj@assays)) {
    stop("Assay provided not fund!")
  }
  #
  if (is.null(feature_sel)) {
    feature_sel <- rownames(seurat_obj)
  } else {
    feature_sel <- intersect(feature_sel, rownames(seurat_obj))
  }
  if (length(feature_sel) == 0) {
    stop("No features provided found!")
  }
  ## format meta data
  seurat_obj@meta.data[[sample_col]] <- gsub("\\_", "\\.", seurat_obj@meta.data[[sample_col]]) %>%
    make.names()
  bulk_meta <- seurat_obj@meta.data[!duplicated(seurat_obj@meta.data[[sample_col]]), 
                                   c(sample_col, block_col, group_col)]
  rownames(bulk_meta) <- bulk_meta[[sample_col]]
  
  ## DEG by clusters
  deg_by_ct <- lapply(all_ct, function(ctx){
    print(ctx)
    ## set subset matrix and meta data
    sel_cells <- colnames(seurat_obj)[seurat_obj@meta.data[[ident_col]] == ctx]
    seurat_objx <- subset(seurat_obj, cells = sel_cells)
    pseudo_dfx <- AverageExpression(seurat_objx, 
                                    assays = assay_use, 
                                    features = feature_sel,
                                    group.by = sample_col, 
                                    return.seurat = F, 
                                    slot = "data")[[assay_use]] %>%
      as.data.frame()
    #
    bulk_metax <- bulk_meta[colnames(pseudo_dfx),]
    sample_tab <- table(bulk_metax[[block_col]])
    pair_sample <- names(sample_tab)[sample_tab == 2]
    
    res <- data.frame(Term = NA,
                      median_diff = NA,
                      pvalue = NA,
                      padj = NA,
                      cluster = ctx)
    if (length((pair_sample) >= 3)) {
      ## split by group and sort by paired samples
      bulk_meta_listx <- lapply(group_levels, function(lvx){
        
        bulk_metaxx <- bulk_metax[bulk_metax[[group_col]] == lvx,]
        bulk_metaxx[match(pair_sample, bulk_metaxx[[block_col]]),]
        
      })
      pseudo_mat_list <- lapply(bulk_meta_listx, function(bulk_metaxx){
        
        pseudo_dfx[,bulk_metaxx[[sample_col]]]
        
      })
      names(pseudo_mat_list) <- names(bulk_meta_listx) <- group_levels
      ## diff = case - control
      avg_diff <- apply(pseudo_mat_list[[group_levels[2]]], 1, mean) - 
        apply(pseudo_mat_list[[group_levels[1]]], 1, mean)
      
      pvalue <- lapply(rownames(pseudo_dfx), function(ftx){
        
        xx <- pseudo_mat_list[[1]][ftx,] %>% as.numeric()
        yy  <- pseudo_mat_list[[2]][ftx,] %>% as.numeric()
        return(wilcox.test(xx, yy, paired = T, alternative = "two.sided")$p.value)
        
      }) %>% unlist
      padj <- p.adjust(pvalue, method = "BH")
      res <- data.frame(Term = rownames(pseudo_dfx),
                        avg_diff = avg_diff,
                        pvalue = pvalue,
                        padj = padj,
                        cluster = ctx)
      
    }
    return(res)
  })
  names(deg_by_ct) <- all_ct
  return(deg_by_ct)
}


