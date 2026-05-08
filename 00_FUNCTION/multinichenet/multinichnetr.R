## module load r/4.5.0
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(nichenetr)
library(multinichenetr)
library(glue)
library(bigreadr)
##### 1. create sce object from seurat object #####
create.from.tenX <- function(tenX_path = NULL,
                             meta_path = NULL,
                             sample_id = NULL,
                             group_id = NULL,
                             celltype_id = NULL,
                             covariates = NA,
                             batch = NA){
  #
  meta_data <- fread2(meta_path)
  count_mat <- Read10X(tenX_path)
  col_data <- meta_data[match(colnames(count_mat), meta_data$V1),
                        c(sample_id, group_id, celltype_id)]
  for (nx in 1:ncol(col_data)) {
    col_data[,nx] <- make.names(col_data[,nx])
  }
  sce <- SingleCellExperiment(assays = list(counts = count_mat),
                              colData = DataFrame(col_data))
  colData(sce)$celltype_id <- make.names(meta_data[[celltype_id]])
  if (!is.na(covariates)) {
    colData(sce)$covariates <- make.names(meta_data[[covariates]])
  }
  if (!is.na(batches)) {
    colData(sce)$batches <- make.names(meta_data[[batches]])
  }
  return(sce)  
}
##### 2. format forward and reverse DEGs #####
format.deg <- function(defile_list = NULL,
                       gene_col = "Term",
                       ct_col = "cluster",
                       fc_col = "log2FC",
                       stat_col = "stat",
                       p_col = "pvalue",
                       padj = "padj",
                       contrasts_oi = NULL){
  #
  celltype_de_fw <- lapply(defile_list, function(x){
    readRDS(x) %>% Reduce("rbind", .)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  # format
  celltype_de_fw <- celltype_de_fw[, c(gene_col, ct_col, fc_col, stat_col, 
                                       p_col, padj, padj)]
  celltype_de_fw$contrast <- contrasts_oi[1]
  colnames(celltype_de_fw) <- c("gene", "cluster_id", "logFC", "F", 
                                "p_val", "p_adj.loc", "p_adj", "contrast")
  celltype_de_fw <- celltype_de_fw[!is.na(celltype_de_fw$p_val) &
                                     celltype_de_fw$cluster_id %in% c(receivers_oi, senders_oi),]
  celltype_de_fw$p_adj.loc[is.na(celltype_de_fw$p_adj.loc) ] <- max(celltype_de_fw$p_adj.loc, na.rm = T)
  celltype_de_fw$p_adj[is.na(celltype_de_fw$p_adj) ] <- max(celltype_de_fw$p_adj, na.rm = T)
  
  # combine reverse and forward
  celltype_de_rev <- celltype_de_fw
  celltype_de_rev$logFC <- -celltype_de_rev$logFC
  celltype_de_rev$contrast <- contrasts_oi[2]
  celltype_de <- rbind(celltype_de_fw, celltype_de_rev)
  
  #
  return(celltype_de)
}

##### 3. Run Multinichenet #####
run.multinichenet <- function(sce = NULL,
                              sample_id = NULL,
                              group_id = NULL,
                              celltype_id = NULL,
                              covariates = NA,
                              batches = NA,
                              senders_oi = NULL,
                              receivers_oi = NULL,
                              contrast_tbl = NULL,
                              celltype_de = NULL,
                              ligand_target_matrix = NULL,
                              lr_network = NULL,
                              weighted_networks = NULL,
                              min_cells = 5,
                              min_sample_prop = 0.5,
                              fraction_cutoff = 0.05,
                              logFC_threshold = 0.5,
                              p_val_threshold = 0.05,
                              p_val_adj = T,
                              top_n_target = 250,
                              n_cores = 1){
  ##### 1. Cell-type filtering #####
  message("1. Run Cell-type filtering......")
  abundance_info <- get_abundance_info(
    sce = sce, 
    sample_id = sample_id, 
    group_id = group_id, 
    celltype_id = celltype_id, 
    min_cells = min_cells, 
    senders_oi = senders_oi, 
    receivers_oi = receivers_oi, 
    batches = batches
  )
  ##
  sample_group_celltype_df <- abundance_info$abundance_data %>% 
    filter(n > min_cells) %>% 
    ungroup() %>% 
    distinct(sample_id, group_id) %>% 
    cross_join(
      abundance_info$abundance_data %>% 
        ungroup() %>% 
        distinct(celltype_id)
    ) %>% 
    arrange(sample_id)
  
  abundance_df <- sample_group_celltype_df %>% 
    left_join(abundance_info$abundance_data %>% ungroup())
  
  abundance_df$n[is.na(abundance_df$n)] <- 0
  abundance_df$keep[is.na(abundance_df$keep)] <- FALSE
  abundance_df_summarized <- abundance_df %>% 
    mutate(keep = as.logical(keep)) %>% 
    group_by(group_id, celltype_id) %>% 
    summarise(samples_present = sum((keep)))
  #
  celltypes_absent_one_condition <- abundance_df_summarized %>% 
    filter(samples_present == 0) %>% 
    pull(celltype_id) %>% 
    unique() 
  celltypes_present_one_condition <- abundance_df_summarized %>% 
    filter(samples_present >= 2) %>% 
    pull(celltype_id) %>% 
    unique() 
  condition_specific_celltypes <- intersect(celltypes_absent_one_condition, 
                                            celltypes_present_one_condition)
  #
  total_nr_conditions <- colData(sce)[,group_id] %>% 
    unique() %>% length() 
  
  absent_celltypes <- abundance_df_summarized %>% 
    filter(samples_present < 2) %>% 
    group_by(celltype_id) %>% 
    dplyr::count() %>% 
    filter(n == total_nr_conditions) %>% 
    pull(celltype_id)
  
  print("condition-specific celltypes:")
  print(condition_specific_celltypes)
  print("absent celltypes:")
  print(absent_celltypes)
  
  analyse_condition_specific_celltypes = FALSE
  if(analyse_condition_specific_celltypes){
    senders_oi <- setdiff(senders_oi,, absent_celltypes)
    receivers_oi <- setdiff(receivers_oi, absent_celltypes)
  } else {
    senders_oi <- setdiff(senders_oi, 
                          union(absent_celltypes, condition_specific_celltypes))
    receivers_oi <- setdiff(receivers_oi,
                            union(absent_celltypes, condition_specific_celltypes))
  }
  
  sce <- sce[, colData(sce)[,celltype_id] %in% 
               c(senders_oi, receivers_oi)]
  ##### 2.Gene filtering: determine #####
  message("2. Run Gene filtering......")
  frq_list <- get_frac_exprs(
    sce = sce, 
    sample_id = sample_id, 
    celltype_id = celltype_id, 
    group_id = group_id, 
    batches = batches, 
    min_cells = min_cells, 
    fraction_cutoff = fraction_cutoff, 
    min_sample_prop = min_sample_prop)
  genes_oi <- frq_list$expressed_df %>% 
    filter(expressed == TRUE) %>% 
    pull(gene) %>% 
    unique() 
  sce <- sce[genes_oi, ]
  
  ##### 3. Pseudobulk expression calculation #####
  message("3. calculate Pseudobulk expression......")
  abundance_expression_info <- process_abundance_expression_info(
    sce = sce, 
    sample_id = sample_id, 
    group_id = group_id, 
    celltype_id = celltype_id, 
    min_cells = min_cells, 
    senders_oi = senders_oi, 
    receivers_oi = receivers_oi, 
    lr_network = lr_network, 
    batches = batches, 
    frq_list = frq_list, 
    abundance_info = abundance_info)
  
  ##### 4. Differential expression analysis #####
  message("4. Get DEGs......")
  if (is.null(celltype_de)) {
    
    DE_info <- get_DE_info(
      sce = sce, 
      sample_id = sample_id, 
      group_id = group_id, 
      celltype_id = celltype_id, 
      batches = batches, 
      covariates = covariates, 
      contrasts_oi = contrasts_oi, 
      min_cells = min_cells, 
      expressed_df = frq_list$expressed_df)
    
    ##
    empirical_pval = T
    if(empirical_pval){
      DE_info_emp <- get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
      celltype_de <- DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
        rename(p_val = p_emp, p_adj = p_adj_emp)
    } else {
      celltype_de <- DE_info$celltype_de$de_output_tidy
    } 
    
  }
  
  ##
  sender_receiver_de <- combine_sender_receiver_de(
    sender_de = celltype_de,
    receiver_de = celltype_de,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network
  )
  ##### 5. Ligand activity prediction #####
  message("5. Run Ligand activity prediction......")
  ligand_activities_targets_DEgenes <- get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = T, 
    n.cores = n_cores
  )
  ##### 6. Prioritization CCC #####
  message("6. Run CCC Prioritization......")
  ligand_activity_down <- FALSE
  sender_receiver_tbl <- sender_receiver_de %>% distinct(sender, receiver)
  metadata_combined <- colData(sce) %>% tibble::as_tibble()
  #
  if(!is.na(batches)){
    grouping_tbl <- metadata_combined[,c(sample_id, group_id, batches)] %>% 
      tibble::as_tibble() %>% distinct()
    colnames(grouping_tbl) <- c("sample", "group", batches)
  } else {
    grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>% 
      tibble::as_tibble() %>% distinct()
    colnames(grouping_tbl) <- c("sample", "group")
  }
  #
  prioritization_tables <- generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    scenario = "regular", # all prioritization criteria will be weighted equally
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender,
    ligand_activity_down = ligand_activity_down
  )
  
  ##### 7. Cross-samples expression correlation #####
  message("7. Run Cross-samples expression correlation......")
  lr_target_prior_cor <- lr_target_prior_cor_inference(
    receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
    abundance_expression_info = abundance_expression_info, 
    celltype_de = celltype_de, 
    grouping_tbl = grouping_tbl, 
    prioritization_tables = prioritization_tables, 
    ligand_target_matrix = ligand_target_matrix, 
    logFC_threshold = logFC_threshold, 
    p_val_threshold = p_val_threshold, 
    p_val_adj = p_val_adj
  )
  
  ##### 8. Save all the output of MultiNicheNet #####
  message("8. Format Output......")
  multinichenet_output <- list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 
  multinichenet_output <- make_lite_output(multinichenet_output)
  return(multinichenet_output)
  
}


