##
library(miloR)
library(SingleCellExperiment)
library(Seurat)
library(scater)
library(scran)
library(bigreadr)
library(dplyr)
library(glue)
library(patchwork)
library(BiocParallel)

##
milor.seurat <- function(seurat_obj = NULL,
                         n_dims = 20,
                         k_nn = 60,
                         group_col = NULL,
                         cov_col = NULL,
                         block_col = NULL,
                         sample_col = NULL,
                         ann_col = NULL,
                         use_reduc = "harmony_SCT",
                         reduc_mat = NULL,
                         plot_reduc = "UMAP",
                         seed_use = 20250528,
                         n_core = 1){
  ## check inputs
  if (any(!c(group_col, cov_col, block_col, sample_col) %in% 
          colnames(seurat_obj@meta.data))) {
    stop("Not all required variables in the metadata of Seurat object!")
  }
  #
  if (!is.null(levels(seurat_obj@meta.data[[group_col]]))) {
    group_levels <- levels(seurat_obj@meta.data[[group_col]])
  } else {
    group_levels <- unique(seurat_obj@meta.data[[group_col]]) %>%sort()
  }
  ## create milo object
  message("Create milo object")
  traj_milo <- as.SingleCellExperiment(seurat_obj) %>%
    Milo()
  #
  traj_milo[[group_col]] <- factor(traj_milo[[group_col]],
                                   levels = group_levels)
  traj_milo[[sample_col]] <- make.names(traj_milo[[sample_col]])
  #
  reducedDim(traj_milo, "UMAP") <- Embeddings(seurat_obj, plot_reduc)
  if (!is.null(reduc_mat)) {
    reducedDim(traj_milo, "ext") <- as.matrix(reduc_mat)
  }
  if (is.null(use_reduc)) {
    use_reduc <- "ext"
    message(glue("Use provided reduction matrix to build graph and make neighbours"))
  } else {
    reducedDim(traj_milo, use_reduc) <- Embeddings(seurat_obj, use_reduc)
  }
  ##
  if(n_dims <= ncol(reducedDim(traj_milo, use_reduc))){
    message(glue("Use {n_dims} dims to build graph and make neighbours"))
  } else {
    n_dims <- ncol(reducedDim(traj_milo, use_reduc))
    message(glue("Max {n_dims} dims present! Use all dims to build graph and make neighbours"))
  }
  # Construct KNN graph
  set.seed(seed_use)
  traj_milo <- buildGraph(traj_milo, 
                          k = k_nn, 
                          d = n_dims,
                          reduced.dim = use_reduc,
                          BPPARAM = MulticoreParam(workers = n_core))
  # Defining representative neighbourhoods
  traj_milo <- makeNhoods(traj_milo, 
                          prop = 0.05, 
                          k = k_nn, 
                          d = n_dims, 
                          reduced_dims = use_reduc,
                          refined = T,
                          refinement_scheme = "graph")
  # Counting cells in neighbourhoods
  traj_milo <- countCells(traj_milo, 
                          meta.data = data.frame(colData(traj_milo)), 
                          samples = sample_col)
  # Differential abundance testing
  traj_design <- data.frame(colData(traj_milo))[, c(group_col, cov_col, block_col, sample_col)]
  traj_design <- distinct(traj_design)
  rownames(traj_design) <- traj_design[[sample_col]]
  # set regression formula
  if (!is.null(block_col)) {
    if (sum(duplicated(traj_design[[block_col]])) >= 40) {
      message(glue("Treat {block_col} as random effect since enough blocks present"))
      rd_term <- glue("(1|{block_col})")
    }
  } else {
    message(glue("Use fix effects only as block not specifify or no enough blocks present"))
    rd_term <- NULL
  }
  formula_milo <- paste0("~", paste(c(cov_col, group_col, rd_term), collapse = " + ")) %>% as.formula()
  
  ## Reorder rownames to match columns of nhoodCounts(milo)
  traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop = F]
  da_results <- testNhoods(traj_milo, 
                           design = formula_milo, 
                           design.df = traj_design, 
                           reduced.dim = use_reduc,
                           fdr.weighting = "graph-overlap",
                           glmm.solver = "Fisher",
                           REML = T,
                           norm.method = "TMM", 
                           BPPARAM = MulticoreParam(workers = n_core))
  ##
  traj_milo <- buildNhoodGraph(traj_milo, overlap = 25)
  plt_milo <- plotNhoodGraphDA(traj_milo, 
                               da_results, 
                               alpha = 0.1) +
    labs(fill = "Log2FC") +
    scale_fill_gradient2(high = scales::muted("red"),
                         mid = "white",
                         low = scales::muted("blue"),
                         midpoint = 0) + 
    plot_layout(guides="auto")
  ##
  if (!is.null(ann_col)) {
    #
    da_results <- annotateNhoods(traj_milo, da_results, coldata_col = ann_col)
    if (!is.null(levels(seurat_obj@meta.data[[ann_col]]))) {
      ann_levels <- intersect(levels(seurat_obj@meta.data[[ann_col]]), unique(da_results[[ann_col]]))
    } else {
      ann_levels <- unique(da_results[[ann_col]]) %>% as.character() %>% sort
    }
    # da_results[[ann_col]] <- ifelse(da_results[[glue("{ann_col}_fraction")]] < 0.7, 
    #                                 "Mixed", 
    #                                 da_results[[ann_col]])
    if (sum(da_results[[ann_col]] == "Mixed") > 0) {
      ann_levels <- c(ann_levels, "Mixed")
    }
    #
    da_results[[ann_col]] <- factor(da_results[[ann_col]], levels = rev(ann_levels))
    plt_milo2 <- plotDAbeeswarm(da_results, group.by = ann_col) +
      ylab("Log2 Fold Change") + 
      xlab("") + 
      theme(axis.title.y = element_text(size = 12, face = "bold"))
    if (sum(da_results$SpatialFDR < 0.05) > 0) {
      plt_milo2 <- plt_milo2 + 
        scale_color_gradient2(high = scales::muted("red"),
                              mid = "white",
                              low = scales::muted("blue"),
                              midpoint = 0)
    } else {
      plt_milo2 <- plt_milo2 + 
        scale_color_manual(values = "white")
    }
  } else {
    plt_milo2 <- NA
  }
  ##
  return(list("plt_milo" = plt_milo,
              "plt_milo2" = plt_milo2,
              "da_results" = da_results))
}

