##
suppressMessages(library(dplyr))
suppressMessages(library(glue))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(variancePartition))
suppressMessages(library(usdm))
suppressMessages(library(BiocParallel))

###### process covariates #####
process.cov <- function(meta_df = NULL,
                        vif_thre = 5,
                        cov_col = NULL){
  #
  meta_df <- apply(meta_df, 2, function(xx){
    
    if (is.character(xx) | is.factor(xx)) {
      xx <- as.factor(xx) %>% as.integer()
    }
    return(xx)
    
  }) %>% as.data.frame()
  #
  var_cov <- lapply(cov_col, function(x){
    
    var(meta_df[[x]])
    
  }) %>% unlist
  cov_col <- cov_col[var_cov > 0]
  #
  if (length(cov_col) > 1) {
    excl_var <- vifstep(meta_df[, cov_col, drop = F], th = 3)@excluded
    cov_col <- setdiff(cov_col, excl_var)
  }
  
  return(cov_col)
  
}

###### DESeq2 ######
DESeq2.fun <- function(exp_mat = NULL,
                       meta_data = NULL,
                       cov_col = NULL,
                       group_col = NULL,
                       group_levels = NULL,
                       deseq_test = "Wald",
                       fit_Type = "parametric",
                       n_core = 4){
  
  ## set DESeq object
  design_formu <- as.formula(paste0("~ ",
                                    paste(c(cov_col, group_col), collapse = " + ")))
  print(design_formu)
  dds <- DESeqDataSetFromMatrix(countData = exp_mat,
                                colData = meta_data,
                                design = design_formu) 
  dds[[group_col]] <- relevel(dds[[group_col]], ref = group_levels[1]) 
  if (deseq_test == "LRT") {
    if (length(cov_col) == 0) {
      reduc_term <- as.formula("~1")
    } else {
      reduc_term <- as.formula(paste0("~ ",
                                      paste(c(cov_col), collapse = " + ")))
    }
  } else {
    reduc_term <- NULL
  }
  ## Run DESeq2
  suppressMessages(dds <- DESeq(dds, 
                                test = deseq_test,
                                fitType = fit_Type,
                                reduced = reduc_term,
                                sfType = "poscounts",
                                betaPrior = F, 
                                parallel = T, 
                                BPPARAM = MulticoreParam(n_core)))
  
  # dds <- estimateSizeFactors(dds,
  #                            type = "poscounts",
  #                            quiet = F)
  # dds <- estimateDispersionsGeneEst(dds, 
  #                                   type = "parametric")
  # dispersions(dds) <- mcols(dds)$dispGeneEst
  # dds <- nbinomWaldTest(dds,
  #                       betaPrior = F,
  #                       useQR = TRUE,
  #                       minmu = 0.5)
  ## summary out
  res <- results(dds, 
                 name = paste0(group_col, "_",
                               group_levels[2], "_vs_", group_levels[1]),
                 independentFiltering = T) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(., "Term")
  res_df <- res[,c("Term", "log2FoldChange", "stat", "pvalue", "padj")]
  colnames(res_df) <- c("Term", "log2FC", "stat", "pvalue", "padj")
  res_df$Contrast <- paste0(group_col, "_",
                            group_levels[2], "_vs_", group_levels[1])
  return(res_df)
  
}

###### limma ######
limma.fun <- function(exp_mat = NULL,
                      meta_data = NULL,
                      cov_col = NULL,
                      block_col = NULL,
                      do_norm = T,
                      group_col = NULL,
                      group_levels = NULL){
  
  ## create DGE object
  if (length(unique(meta_data[[group_col]])) != 2 |
      length(group_levels) != 2) {
    message("Warning: Group should only have TWO levels!")
    res_df <- data.frame("Term" = NA, 
                         "log2FC" = NA, 
                         "stat" = NA, 
                         "pvalue" = NA, 
                         "padj" = NA)
  } else {
    ## design matrix
    meta_data$Group <- meta_data[[group_col]] %>%
      factor(., levels = group_levels, labels = c("Control", "Case"))
    design_formu <- as.formula(paste0("~ 0 +",
                                      paste(c("Group", cov_col), collapse = " + ")))
    design <- model.matrix(design_formu, data = meta_data)
    rownames(design) <- rownames(meta_data)
    print(design_formu)
    
    ## contrast matrix
    contrast_matrix <- makeContrasts(
      "contra_g" = GroupCase - GroupControl,
      levels = design)
    
    ## sample cells for duplicateCorrelation
    if(!is.null(block_col)){
      
      use_block <- meta_data[[block_col]]
      if(sum(duplicated(use_block)) >= 2){
        if (do_norm) {
          lcpm <- edgeR::cpm(exp_mat, log = T)
        } else {
          lcpm <- exp_mat
        }
        cor <- duplicateCorrelation(lcpm, 
                                    design, 
                                    block = use_block)
        cor_block <- cor$consensus.correlation
        print(cor_block)
        
      } else {
        # with less than two pairs
        use_block <- cor_block <- NULL
      }
      
    } else {
      use_block <- cor_block <- NULL
    }
    
    ## fit lmm
    if (do_norm) {
      dge_list <- DGEList(counts = exp_mat, samples = meta_data)
      dge_list <- calcNormFactors(dge_list, method = "TMM")
      v <- voom(dge_list, 
                design = design, 
                block = use_block, 
                correlation = cor_block)
    } else {
      v <- exp_mat
    }
    
    fit_lmm <- lmFit(object = v, 
                     design = design, 
                     block = use_block, 
                     correlation = cor_block)
    fit_lmm2 <- contrasts.fit(fit_lmm, contrast_matrix) 
    fit_lmm2 <- eBayes(fit_lmm2)
    ## summary out
    res <- topTable(fit_lmm2, coef = 1, number = Inf, adjust = "fdr") %>%
      as.data.frame() %>%
      tibble::rownames_to_column(., "Term")
    res_df <- res[,c("Term", "logFC", "t", "P.Value", "adj.P.Val")]
    colnames(res_df) <- c("Term", "log2FC", "stat", "pvalue", "padj")
  }
  res_df$Contrast <- paste0(group_col, "_",
                            group_levels[2], "_vs_", group_levels[1])
  return(res_df)
}

###### limma-dream ######
limma.dream.fun <- function(exp_mat = NULL,
                            meta_data = NULL,
                            cov_col = NULL,
                            block_col = NULL,
                            group_col = NULL,
                            group_levels = NULL,
                            n_core = 4){
  
  ## create DGE object
  if (length(unique(meta_data[[group_col]])) != 2 |
      length(group_levels) != 2) {
    message("Warning: Group should only have TWO levels!")
    res_df <- data.frame("Term" = NA, 
                         "log2FC" = NA, 
                         "stat" = NA, 
                         "pvalue" = NA, 
                         "padj" = NA)
  } else {
    meta_data$Group <- meta_data[[group_col]] %>%
      factor(., levels = group_levels, labels = c("Control", "Case"))
    ## set regresion formula
    if (!is.null(block_col)) {
      rd_term <- glue("(1|{block_col})")
    } else {
      rd_term <- NULL
    }
    formula_dream <- paste0("~ 0 + ", paste(c("Group", cov_col, rd_term), collapse = " + ")) %>% as.formula()
    message(formula_dream)
    ## contrast matrix
    L <- makeContrastsDream(formula = formula_dream, 
                            data = meta_data,
                            contrasts = c(contra_g = "GroupCase - GroupControl")
    )
    ## fit lmm
    if (do_norm) {
      dge_list <- DGEList(counts = exp_mat, samples = meta_data)
      dge_list <- calcNormFactors(dge_list, method = "TMM")
    } else {
      dge_list <- exp_mat
    }
    use_param <- SnowParam(n_core, 
                           "SOCK", 
                           progressbar = T)
    vobjDream <- voomWithDreamWeights(counts = dge_list, 
                                      formula = formula_dream, 
                                      data = meta_data,
                                      plot = F,
                                      save.plot = F,
                                      rescaleWeightsAfter = F,
                                      scaledByLib = F,
                                      priorWeightsAsCounts = F,
                                      normalize.method = "none",
                                      BPPARAM = use_param)
    fitmm <- dream(exprObj = vobjDream,
                   formula = formula_dream,
                   data = meta_data,
                   L = L,
                   ddf = "adaptive",
                   useWeights = T,
                   hideErrorsInBackend = F,
                   REML = T,
                   BPPARAM = use_param)
    fitmm <- eBayes(fitmm)
    ## summary out
    res <- topTable(fitmm, coef = "contra_g", number = Inf, adjust = "fdr") %>%
      as.data.frame() %>%
      tibble::rownames_to_column(., "Term")
    res_df <- res[,c("Term", "logFC", "t", "P.Value", "adj.P.Val")]
    colnames(res_df) <- c("Term", "log2FC", "stat", "pvalue", "padj")
  }
  res_df$Contrast <- paste0(group_col, "_",
                            group_levels[2], "_vs_", group_levels[1])
  return(res_df)
}


###### edgeR ######
edgeR.fun <- function(exp_mat = NULL,
                      meta_data = NULL,
                      edger_test = "QLF",
                      cov_col = NULL,
                      group_col = NULL,
                      group_levels = NULL){
  
  ## create DGE object
  if (length(unique(meta_data[[group_col]])) != 2 |
      length(group_levels) != 2) {
    message("Warning: Group should only have TWO levels!")
    res_df <- data.frame("Term" = NA, 
                         "log2FC" = NA, 
                         "stat" = NA, 
                         "pvalue" = NA, 
                         "padj" = NA)
  } else {
    meta_data$Group <- meta_data[[group_col]] %>%
      factor(., levels = group_levels, labels = c("Control", "Case"))
    ##
    dge_list <- DGEList(counts = exp_mat, samples = meta_data)
    dge_list$samples$lib.size <- colSums(dge_list$counts)
    dge_list <- normLibSizes(dge_list, method = "TMM")
    ## design matrix
    design_formu <- as.formula(paste0("~ 0 +",
                                      paste(c(cov_col, "Group"), collapse = " + ")))
    design <- model.matrix(design_formu, data = meta_data)
    dimnames(design) <- list(rownames(meta_data),
                             make.names(colnames(design)))
    ##
    dge_list <- estimateDisp(dge_list, design, robust = T)
    res_df <- switch(edger_test,
                     "QLF" = {
                       fit <- glmQLFit(dge_list, design)
                       test <- glmQLFTest(fit, coef = ncol(fit$design))
                       res <- topTags(test, n = Inf, sort.by = "PValue")$table %>%
                         as.data.frame() %>%
                         tibble::rownames_to_column(., "Term")
                       res_df <- res[,c("Term", "logFC", "F", "PValue", "FDR")]
                     },
                     "LRT" = {
                       fit <- glmFit(dge_list, design = design)
                       test <- glmLRT(fit, coef = ncol(fit$design))
                       res <- topTags(test, n = Inf, sort.by = "PValue")$table %>%
                         as.data.frame() %>%
                         tibble::rownames_to_column(., "Term")
                       res_df <- res[,c("Term", "logFC", "LR", "PValue", "FDR")]
                     })
    ## summary out
    colnames(res_df) <- c("Term", "log2FC", "stat", "pvalue", "padj")
  }
  res_df$Contrast <- paste0(group_col, "_",
                            group_levels[2], "_vs_", group_levels[1])
  return(res_df)
  
}

###### pseudoBulk ######
pseudo.bulk.deg <- function(seurat_obj = NULL,
                            use_assay = "RNA",
                            ident_col = NULL,
                            group_col = NULL,
                            cov_col = NULL,
                            agg_strategy = "sum",
                            vif_thre = 5,
                            sample_col = NULL,
                            block_col = NULL,
                            paired_only = T,
                            do_norm = T,
                            min_cell_per_sample = 10,
                            min_percent = 0.1,
                            filter_ByExpr = F,
                            test_use = "DESeq2-Wald",
                            n_core = 4){
  ##
  if (any(!c(ident_col, group_col, cov_col, block_col, sample_col) %in% 
          colnames(seurat_obj@meta.data))) {
    stop("Not all required variables in the metadata of Seurat object!")
  }
  # set cell type levels
  if (!is.null(levels(seurat_obj@meta.data[[ident_col]]))) {
    all_ct <- seurat_obj@meta.data[[ident_col]] %>% droplevels() %>% levels()
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
  if (!grepl("DESeq2|limma|edgeR", test_use)) {
    stop("Please use DEseq2, limma, or edgeR for DEG identification!")
  }
  ## format meta data
  seurat_obj@meta.data[[sample_col]] <- gsub("\\_", "\\.", seurat_obj@meta.data[[sample_col]]) %>%
    make.names()
  bulk_meta <- seurat_obj@meta.data[!duplicated(seurat_obj@meta.data[[sample_col]]), 
                                    c(sample_col, block_col, group_col, cov_col)]
  rownames(bulk_meta) <- bulk_meta[[sample_col]]
  
  ## DEG by clusters
  deg_by_ct <- lapply(all_ct, function(ctx){
    print(ctx)
    ## set subset matrix and meta data
    seurat_objx <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[ident_col]] == ctx])
    ## removed samples with less than min_cell_per_sample cells
    sample_countx <- table(seurat_objx@meta.data[[sample_col]])
    sample_retain <- names(sample_countx)[sample_countx >= min_cell_per_sample]
    ## at least 3+3 samples needed
    if (length(sample_retain) >= 10) {
      
      seurat_objx <- subset(seurat_objx, cells = colnames(seurat_objx)[seurat_objx@meta.data[[sample_col]] %in% sample_retain])
      #
      if (agg_strategy %in% c("sum", "SUM", "Sum")) {
        pseudo_matx <- AggregateExpression(seurat_objx, 
                                           assays = use_assay, 
                                           group.by = sample_col, 
                                           return.seurat = F, 
                                           slot = "counts")[[use_assay]] %>%
          as.matrix()
      } else {
        pseudo_matx <- AverageExpression(seurat_objx, 
                                         assays = use_assay, 
                                         group.by = sample_col, 
                                         return.seurat = F, 
                                         slot = "data")[[use_assay]] %>%
          as.data.frame()
      }
      bulk_metax <- bulk_meta[colnames(pseudo_matx),]
      
      ## filtering when enabling block/paired design
      if (!is.null(block_col)) {
        if (grepl("DESeq2|edgeR", test_use)) {
          
          sample_pair <- bulk_metax[[block_col]][duplicated(bulk_metax[[block_col]])]
          
          if (length(sample_pair) >= 3) {
            
            ## keep paired samples only for DESeq2 and edgeR
            message(glue("Only retained paired samples as {test_use} does not support LMM!"))
            bulk_metax <- bulk_metax[bulk_metax[[block_col]] %in% sample_pair,]
            pseudo_matx <- pseudo_matx[,rownames(bulk_metax)]
            # then treat block as covariates only
            cov_col <- block_col
            
          } else {
            message(glue("Less than Three sample pairs present! Retained unpaired samples to adopt CR design!"))
            max_groupx <- table(bulk_metax[[group_col]]) %>% sort() %>% tail(1) %>% names()
            dup_idx <- which(bulk_metax[[block_col]] %in% sample_pair & bulk_metax[[group_col]] == max_groupx)
            if (length(dup_idx) > 0) {
              bulk_metax <- bulk_metax[-dup_idx,]
              pseudo_matx <- pseudo_matx[,rownames(bulk_metax)]
            }
          }
          
        } else {
          
          if (paired_only) {
            ## keep paired samples only
            message("Only retained paired samples!")
            sample_pair <- bulk_metax[[block_col]][duplicated(bulk_metax[[block_col]])]
            bulk_metax <- bulk_metax[bulk_metax[[block_col]] %in% sample_pair,]
            pseudo_matx <- pseudo_matx[,rownames(bulk_metax)]
            
            
          } else {
            message("Use all samples!")
          }
          
        }
        ## keep all samples for limma-lmm
        bulk_metax[[block_col]] <- as.factor(bulk_metax[[block_col]])
      }
      
      ## remove redundant cov
      cov_col <- process.cov(meta_df = bulk_metax,
                             vif_thre = vif_thre,
                             cov_col = cov_col)
      
      ## filtering genes (take the library size in to account)
      var_pseudo_matx <- apply(pseudo_matx, 1, function(x){var(x)}) %>% unlist
      var_exprs <- which(var_pseudo_matx > 0)
      if (filter_ByExpr) {
        ## filter by cpm
        if (ncol(pseudo_matx) * min_percent < 3) {
          min_percent_use <- 3/ncol(pseudo_matx)
        } else {
          min_percent_use <- min_percent
        }
        med_lib <- colSums(pseudo_matx) %>% median()
        min_cpm <- 10
        min_count <- ifelse(ceiling(min_cpm * med_lib / 1E6) > 10, 10,
                            ceiling(min_cpm * med_lib / 1E6))
        min_total_count <- ifelse(3 * min_count < 10, 10,
                                  3 * min_count)
        keep_exprs <- filterByExpr(pseudo_matx,
                                   group = bulk_metax[[group_col]],
                                   min.count = min_count,
                                   min.total.count = min_total_count,
                                   min.prop = min_percent)
        keep_exprs <- which(keep_exprs)
      } else {
        ## filter by counts (proportion with expression)
        keep_exprs <- split(t(pseudo_matx) %>% as.data.frame, f = bulk_metax[[group_col]]) %>%
          lapply(., function(x){
            
            if (min_percent > 0 & nrow(x) * min_percent < 3) {
              min_percent_use <- 3/nrow(x)
            } else {
              min_percent_use <- min_percent
            }
            
            which(colSums(x >= 1) >= nrow(x) * min_percent_use)
          }) %>% unlist %>% unique
        # keep_exprs <- rowSums(pseudo_matx >= 5) >= ncol(pseudo_matx) * min_percent
      }
      pseudo_matx <- pseudo_matx[intersect(var_exprs, keep_exprs),]
      message(glue("{ncol(pseudo_matx)} samples and {nrow(pseudo_matx)} features were included!"))
      #
      if (min(table(bulk_metax[[group_col]])) >= 3 &
          nrow(pseudo_matx) > 1) {
        message(glue("{ncol(pseudo_matx)} samples and {nrow(pseudo_matx)} features were included!"))
        ## Run DEG analysis
        res <- tryCatch(
          expr = {
            res_try <- switch(test_use,
                   "DESeq2-LRT" = {
                     DESeq2.fun(exp_mat = pseudo_matx,
                                meta_data = bulk_metax,
                                cov_col = cov_col,
                                group_col = group_col,
                                group_levels = group_levels,
                                deseq_test = "LRT",
                                fit_Type = "glmGamPoi",
                                n_core = n_core)
                   },
                   "DESeq2-Wald" = {
                     DESeq2.fun(exp_mat = pseudo_matx,
                                meta_data = bulk_metax,
                                cov_col = cov_col,
                                group_col = group_col,
                                group_levels = group_levels,
                                deseq_test = "Wald",
                                fit_Type = "parametric",
                                n_core = n_core)
                   },
                   "limma" = {
                     limma.fun(exp_mat = pseudo_matx,
                               meta_data = bulk_metax,
                               cov_col = cov_col,
                               group_col = group_col,
                               block_col = block_col,
                               do_norm = do_norm,
                               group_levels = group_levels)
                   },
                   "limma-dream" = {
                     limma.dream.fun(exp_mat = pseudo_matx,
                                     meta_data = bulk_metax,
                                     cov_col = cov_col,
                                     block_col = block_col,
                                     group_col = group_col,
                                     group_levels = group_levels,
                                     n_core = n_core)
                   },
                   "edgeR-LRT" = {
                     edgeR.fun(exp_mat = pseudo_matx,
                               meta_data = bulk_metax,
                               edger_test = "LRT",
                               cov_col = cov_col,
                               group_col = group_col,
                               group_levels = group_levels)
                   },
                   "edgeR-QLF" = {
                     edgeR.fun(exp_mat = pseudo_matx,
                               meta_data = bulk_metax,
                               edger_test = "QLF",
                               cov_col = cov_col,
                               group_col = group_col,
                               group_levels = group_levels)
                   })
            res_try$cluster <- ctx
            res_try$cov_use <- paste(cov_col, collapse = ",")
            return(res_try)
          },
          error = function(e) {
            data.frame(Term = NA, 
                       log2FC = NA, 
                       stat = NA, 
                       pvalue = NA, 
                       padj = NA,
                       Contrast = paste0(group_col, "_",
                                         group_levels[2], "_vs_", group_levels[1]),
                       cluster = ctx,
                       cov_use = NA)
          }
        )

        
      } else {
        res <- data.frame(Term = NA, 
                          log2FC = NA, 
                          stat = NA, 
                          pvalue = NA, 
                          padj = NA,
                          Contrast = paste0(group_col, "_",
                                            group_levels[2], "_vs_", group_levels[1]),
                          cluster = ctx,
                          cov_use = NA)
        
      }
      
    } else {
      res <- data.frame(Term = NA, 
                        log2FC = NA, 
                        stat = NA, 
                        pvalue = NA, 
                        padj = NA,
                        Contrast = paste0(group_col, "_",
                                          group_levels[2], "_vs_", group_levels[1]),
                        cluster = ctx,
                        cov_use = NA)
      
    }
    return(res)
  })
  # deg_by_ct <- lapply(all_ct, function(ctx){
  #   readRDS(glue("03_output/04_Diff/DAR/small_intestine/pseudo_bulk_N_U_vs_N_N_{ctx}.rds"))
  # })
  names(deg_by_ct) <- all_ct
  return(deg_by_ct)
}

