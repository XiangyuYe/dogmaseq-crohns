#
library(monocle)
library(pheatmap)
library(glue)
library(dplyr)
library(glue)
library(ggplot2)

###### monocle2.fun ######
monocle2.fun <- function(seurat_obj = NULL,
                         order_genes = NULL,
                         cluster_col = NULL,
                         group_col = NULL,
                         n_sample = 20000,
                         sed_use = 20250323,
                         out_path = NULL){
  ## 0. set parameters
  # 0.1 set marker gene per cluster for ordering alternatively (use deg from seurat wilcox test)
  if (is.null(order_genes)) {
    message("Perform DEGs as no genes provided for ordering!")
    Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col]]
    seurat_obj <- NormalizeData(seurat_obj)
    sc_markers <- FindAllMarkers(seurat_obj, 
                                 assay = "RNA", 
                                 only.pos = T, 
                                 return.thresh = 1E-3,
                                 test.use = "wilcox")
    sc_markers <- subset(sc_markers, sc_markers$p_val_adj < 0.05 &
                           sc_markers$avg_log2FC > 1)
    top_markers <- sc_markers %>%
      group_by(cluster) %>%
      dplyr::top_n(n = 200, wt = -p_val_adj) %>%
      dplyr::top_n(n = 200, wt = avg_log2FC) %>%
      ungroup()
    order_genes <- unique(top_markers$gene)
  }
  # 0.2 sample cells to ensure normal running of monocle2
  if (ncol(seurat_obj) > n_sample) {
    message(glue("Sample {n_sample} cells to ensure normal running of monocle2!"))
    set.seed(sed_use)
    sel_cells <- sample(colnames(seurat_obj), n_sample)
    seurat_obj <- subset(seurat_obj, cells = sel_cells)
  }
  if (length(seurat_obj@reductions) == 0) {
    seurat_obj <- NormalizeData(seurat_obj) %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA()
  }
  ## 1. Run monocle2
  # 1.1 creat cds object
  cds <- as.CellDataSet(seurat_obj, 
                        assay  = "RNA")
  # 1.2 Estimate size factor
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  # 1.3 set genes for Ordering
  cds <- monocle::setOrderingFilter(cds, order_genes)
  # 1.4 dimension reduciton
  cds <- monocle::reduceDimension(cds, method = 'DDRTree')
  # 1.5 ordering cells
  cds <- monocle::orderCells(cds)
  reduc_ddr <- cds@reducedDimS %>% t %>% as.data.frame()
  colnames(reduc_ddr) <- c("X", "Y")
  traj_df <- cbind(data.frame(Pseudotime = cds$Pseudotime,
                              cluster = cds[[cluster_col]],
                              Group = cds[[group_col]]),
                   reduc_ddr) %>%
    tibble::rownames_to_column(., "cells")
  saveRDS(cds, file = glue("{out_path}/cds_monocle2.rds"))
  saveRDS(traj_df, glue("{out_path}/Pseudotime_monocle2.rds"))
  
  message("Monocle2 complete!")
  return(cds)
}

###### monocle.plot ######
monocle.plot <- function(monocle_obj = NULL,
                         ct_cols = NULL,
                         cluster_col = cluster_col){
  
  # 2.1 set theme for plot
  my_theme <- theme(legend.position = "right",
                    legend.title = element_text(face = "bold"),
                    panel.spacing.x = unit(0, "cm"), 
                    panel.spacing.y = unit(0, "cm"), 
                    strip.placement = "outside", 
                    strip.background = element_blank(),
                    strip.text.x = element_text(size = 12, 
                                                face = "bold", 
                                                color = "black", 
                                                angle = 0),
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    panel.grid = element_blank())
  # 2.2 plot for pseudotime
  traj_plt1 <- monocle::plot_cell_trajectory(monocle_obj, 
                                             color_by = "Pseudotime",
                                             cell_size = 1,
                                             show_branch_points = F) +
    scale_color_viridis_c()+
    theme_bw() + 
    my_theme
  # 2.2 colored by clusters
  traj_plt2 <- monocle::plot_cell_trajectory(monocle_obj, 
                                             color_by = cluster_col,
                                             cell_size = 1,
                                             show_branch_points = F) +
    scale_color_manual(name = "Cell types",
                       values = ct_cols)+
    guides(color = guide_legend(title = "Cell Type", 
                                byrow = F, 
                                ncol = 1,
                                override.aes = list(size = 3)))+
    theme_bw() + 
    my_theme
  
  return(list("traj_plt1" = traj_plt1,
              "traj_plt2" = traj_plt2))
}

###### make_add_annotation_col ######
make_add_annotation_col <- function(cds_subset, cols, bins = 100) {
  stopifnot("Pseudotime" %in% colnames(pData(cds_subset)))
  pd <- pData(cds_subset)
  #
  brks <- seq(min(pd$Pseudotime, na.rm = TRUE),
              max(pd$Pseudotime, na.rm = TRUE),
              length.out = bins + 1)
  # 
  bin_id <- cut(pd$Pseudotime, breaks = brks, include.lowest = TRUE, labels = FALSE)
  
  # 
  Mode <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    tab <- table(x)
    names(tab)[which.max(tab)]
  }
  
  ann_list <- lapply(cols, function(cl) {
    v <- pd[[cl]]
    if (is.numeric(v)) {
      tapply(v, bin_id, function(z) mean(z, na.rm = TRUE))
    } else {
      tapply(v, bin_id, Mode)
    }
  })
  ann <- as.data.frame(ann_list, stringsAsFactors = FALSE)
  colnames(ann) <- cols
  rownames(ann) <- seq_len(bins)     
  return(ann)
}

##### plot_pseudotime_heatmap_modify #####
plot_pseudotime_heatmap_modify <- function(cds_subset, 
                                           cluster_rows = TRUE, 
                                           hclust_method = "ward.D2", 
                                           num_clusters = 6, 
                                           hmcols = NULL, 
                                           add_annotation_row = NULL, 
                                           add_annotation_col = NULL, 
                                           annotation_colors = NULL,
                                           show_rownames = FALSE, 
                                           use_gene_short_name = TRUE, 
                                           norm_method = c("log", "vstExprs"), 
                                           scale_max = 3, 
                                           scale_min = -3, 
                                           trend_formula = "~sm.ns(Pseudotime, df=3)", 
                                           return_heatmap = FALSE, 
                                           cores = 1){
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                         max(pData(cds_subset)$Pseudotime), length.out = 100))
  m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, 
                       relative_expr = T, new_data = newdata)
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && 
      is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    m = vstExprs(cds_subset, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log10(m + pseudocount)
  }
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- monocle:::blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T, 
                 cluster_cols = FALSE, 
                 cluster_rows = cluster_rows, 
                 show_rownames = F, 
                 show_colnames = F, 
                 clustering_distance_rows = row_dist, 
                 clustering_method = hclust_method, 
                 cutree_rows = num_clusters, 
                 silent = TRUE, 
                 filename = NA, 
                 breaks = bks, 
                 border_color = NA, 
                 color = hmcols)
  if (cluster_rows) {
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                         num_clusters)))
  }
  else {
    annotation_row <- NULL
  }
  if (!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
    ])
    colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
  }
  if (!is.null(add_annotation_col)) {
    if (nrow(add_annotation_col) != 100) {
      stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
    }
    annotation_col <- add_annotation_col
  }
  else {
    annotation_col <- NA
  }
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 
                                                      "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 
                                                       "gene_short_name"])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    if (!is.null(annotation_row)) 
      row_ann_labels <- row.names(annotation_row)
  }
  row.names(heatmap_matrix) <- feature_label
  if (!is.null(annotation_row)) 
    row.names(annotation_row) <- row_ann_labels
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  ph_res <- pheatmap(heatmap_matrix[, ], 
                     useRaster = T, 
                     cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, 
                     show_rownames = show_rownames, 
                     show_colnames = F, 
                     clustering_distance_rows = row_dist, 
                     clustering_method = hclust_method, 
                     cutree_rows = num_clusters, 
                     annotation_row = annotation_row, 
                     annotation_col = annotation_col, 
                     annotation_colors = annotation_colors,
                     treeheight_row = 20, 
                     breaks = bks, 
                     fontsize = 6, 
                     color = hmcols, 
                     border_color = NA, 
                     silent = TRUE, 
                     filename = NA)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap) {
    return(ph_res)
  }
}

