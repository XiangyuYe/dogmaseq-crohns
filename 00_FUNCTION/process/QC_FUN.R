## functions for quality control
library(bigreadr)
library(stringr)
library(SoupX)
library(BiocParallel)
library(dsb)
suppressMessages(library(dplyr))
suppressMessages(library(glue))
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(scDblFinder))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Matrix))
suppressMessages(library(mclust))

##### Demultiplex and QC on HTO data #####
qc.fun.HTO <- function(hto_path = NULL,
                       meta_barcode = NULL,
                       max_ncounts = 2E4,
                       valid_hto = NULL,
                       idx = NULL,
                       add_idx = T,
                       out_path = NULL,
                       save_rds = T,
                       return_seurat = F){
  ## load data and create object
  mtx_file <- glue('{hto_path}/featurecounts.mtx')
  feature_file <- glue('{hto_path}/featurecounts.genes.txt')
  barcode_file <- glue('{hto_path}/featurecounts.barcodes.txt')
  
  if (all(file.exists(mtx_file, feature_file, barcode_file))) {
    #
    hto_mat_raw <- readMM(mtx_file) %>% t
    dimnames(hto_mat_raw) <- list(fread2(feature_file, header = FALSE)[[1]],
                                  fread2(barcode_file, header = FALSE)[[1]])
    colnames(hto_mat_raw) <- paste0(colnames(hto_mat_raw), "-1")
    hto_mat_raw <- hto_mat_raw[, colnames(hto_mat_raw) %in% meta_barcode]
    
    scHTO_obj <- CreateAssayObject(hto_mat_raw, 
                                   min.cells = 0, 
                                   min.features = 0) %>%
      CreateSeuratObject(counts = .,
                         assay = "HTO", 
                         project = idx)
    if(!is.null(valid_hto)){
      scHTO_obj <- subset(scHTO_obj, features = valid_hto)
    }
    n_forclus <- ifelse(ncol(scHTO_obj) > 1000, 1000, ncol(scHTO_obj) * 0.1)
    scHTO_obj <- NormalizeData(scHTO_obj, 
                               assay = "HTO", 
                               nsamples = n_forclus,
                               normalization.method = "CLR")
    keep_idx <- apply(scHTO_obj[["HTO"]]$data, 1, var)
    scHTO_obj <- subset(scHTO_obj, features = names(keep_idx)[keep_idx > 0.1])
    scHTO_obj <-  HTODemux(scHTO_obj, 
                           assay = "HTO", 
                           positive.quantile = 0.99)
    if(!is.null(valid_hto)){
      scHTO_obj$HTO_classification.global[scHTO_obj$HTO_classification.global == "Singlet" &
                                            !scHTO_obj$hash.ID %in% valid_hto] <- "Mismatched"
    }
    Idents(scHTO_obj) <- scHTO_obj$orig.ident
    
    # Visualize QC metrics as a violin plot
    pre_violin <- VlnPlot(scHTO_obj, 
                          features = c("nCount_HTO", "nFeature_HTO", "HTO_margin"), 
                          ncol = 3, 
                          raster = T)
    ##
    scHTO_obj$pass_HTO <- scHTO_obj$nCount_HTO < max_ncounts & 
      scHTO_obj$HTO_classification.global == "Singlet"
    post_violin <- VlnPlot(subset(scHTO_obj, pass_HTO),
                           features = c("nCount_HTO", "nFeature_HTO", "HTO_margin"),
                           ncol = 3,
                           raster = T)
    
    tiff(glue("{out_path}/{idx}_qc_hto.tiff"), 
         height = 15, width = 15, units = "in", compression = "lzw", res = 300)
    (pre_violin/post_violin) %>% print()
    dev.off()
    pdf(glue("{out_path}/{idx}_heat_hto.pdf"))
    (HTOHeatmap(scHTO_obj, assay = "HTO") + scale_fill_viridis_c()) %>% print()
    dev.off()
    
    if(add_idx){
      colnames(scHTO_obj) <- paste0(idx, "_", colnames(scHTO_obj))
    }
    if (save_rds) {
      saveRDS(scHTO_obj, file = glue("{out_path}/{idx}_hto_tmp.rds"))
    }
    if (return_seurat) {
      return(scHTO_obj)
    } else {
      return(colnames(scHTO_obj))
    }
    
  } else {
    
    return(NULL)
    
  }
  
}

##### QC on ADT data  #####
qc.fun.ADT <- function(adt_path = NULL,
                       meta_barcode = NULL,
                       min_features = 80,
                       max_ncounts_ctrl = 50,
                       max_percent_ctrl = 2,
                       max_ncounts = 2E4,
                       min_ncounts = 1E3,
                       idx = NULL,
                       add_idx = T,
                       out_path = NULL,
                       save_rds = T,
                       return_seurat = F){
  ## load data and create object
  mtx_file <- glue('{adt_path}/featurecounts.mtx')
  feature_file <- glue('{adt_path}/featurecounts.genes.txt')
  barcode_file <- glue('{adt_path}/featurecounts.barcodes.txt')
  
  if (all(file.exists(mtx_file, feature_file, barcode_file))) {
    
    adt_mat_raw <- readMM(mtx_file) %>% t
    dimnames(adt_mat_raw) <- list(fread2(feature_file, header = FALSE)[[1]] %>%
                                    str_split_i(., "-A[0-1]", 1),
                                  fread2(barcode_file, header = FALSE)[[1]] %>%
                                    paste0(., "-1"))
    ##     
    scADT_obj <- CreateAssayObject(adt_mat_raw[, colnames(adt_mat_raw) %in% meta_barcode], 
                                   min.features = 0, 
                                   min.cells = 0) %>%
      CreateSeuratObject(counts = .,
                         assay = "ADT", 
                         project = idx)
    scADT_obj[["percent.ctrl"]] <- PercentageFeatureSet(scADT_obj, pattern = "Ctrl$")
    scADT_obj[["sum.ctrl"]] <- scADT_obj[["percent.ctrl"]] / 100 * scADT_obj$nCount_ADT
    # Visualize QC metrics as a violin plot
    pre_violin <- VlnPlot(scADT_obj, 
                          features = c("nCount_ADT", "nFeature_ADT", "sum.ctrl"),
                          ncol = 3, 
                          raster = T)
    ## 
    scADT_obj$pass_ADT <- scADT_obj$nFeature_ADT > min_features &
      scADT_obj$nCount_ADT < max_ncounts &
      scADT_obj$nCount_ADT > min_ncounts
    if (!is.null(max_percent_ctrl)) {
      scADT_obj$pass_ADT <- scADT_obj$pass_ADT &
        scADT_obj$percent.ctrl < max_percent_ctrl
    }
    if (!is.null(max_ncounts_ctrl)) {
      scADT_obj$pass_ADT <- scADT_obj$pass_ADT &
        scADT_obj$sum.ctrl < max_ncounts_ctrl
    }
    scADT_obj$pass_ADT[is.na(scADT_obj$pass_ADT)] <- FALSE
    ##
    post_violin <- VlnPlot(subset(scADT_obj, pass_ADT),
                           features = c("nCount_ADT", "nFeature_ADT", "sum.ctrl"),
                           ncol = 3,
                           raster = T)
    
    tiff(glue("{out_path}/{idx}_qc_adt.tiff"), 
         height = 15, width = 15, units = "in", compression = "lzw", res = 300)
    (pre_violin/post_violin) %>% print()
    dev.off()
    
    #
    if(add_idx){
      colnames(scADT_obj) <- paste0(idx, "_", colnames(scADT_obj))
    }
    
    if (save_rds) {
      saveRDS(scADT_obj, file = glue("{out_path}/{idx}_adt_tmp.rds"))
    }
    if (return_seurat) {
      return(scADT_obj)
    } else {
      return(colnames(scADT_obj))
    }
    
  } else {
    
    return(NULL)
    
  }
  
}
## denoise on ADT data
denoise.pre.ADT <- function(adt_path = NULL,
                            clean_barcode = NULL,
                            blank_barcode = NULL,
                            idx = NULL,
                            out_path = NULL){
  ## load data and create object
  mtx_file <- glue('{adt_path}/featurecounts.mtx')
  feature_file <- glue('{adt_path}/featurecounts.genes.txt')
  barcode_file <- glue('{adt_path}/featurecounts.barcodes.txt')
  
  if (all(file.exists(mtx_file, feature_file, barcode_file))) {
    
    adt_mat_raw <- readMM(mtx_file) %>% t
    dimnames(adt_mat_raw) <- list(fread2(feature_file, header = FALSE)[[1]] %>%
                                    str_split_i(., "-A[0-1]", 1),
                                  fread2(barcode_file, header = FALSE)[[1]] %>%
                                    paste0(., "-1"))
    ## denoise
    pre_dsb_df <- data.frame(prot.size = log10(colSums(adt_mat_raw)),
                             group = ifelse(colnames(adt_mat_raw) %in% blank_barcode,
                                            "bk", 
                                            ifelse(colnames(adt_mat_raw) %in% clean_barcode, 
                                                   "cell", NA)))
    pre_dsb_df <- pre_dsb_df[pre_dsb_df$prot.size >= 0 & !is.na(pre_dsb_df$group),]
    blank_pick <- pick.blank(pre_dsb_df, 
                             log_col = "prot.size", 
                             group_col = "group")
    pre_dsb_plt <- ggplot(pre_dsb_df) + 
      geom_histogram(aes(x = prot.size, fill = group), bins = 100) + 
      geom_vline(xintercept = blank_pick$suggested_interval$lo, 
                 linetype = "dashed", colour = "blue") +
      geom_vline(xintercept = blank_pick$suggested_interval$hi, 
                 linetype = "dashed", colour = "blue") +
      geom_vline(xintercept = blank_pick$peak_stats$t_cut, 
                 linetype = "dotted", colour = "red") +
      scale_x_continuous(minor_breaks = seq(1, 3, 0.2)) + 
      theme_bw()
    ##
    ## plot
    pdf(glue("{out_path}/{idx}_pre_dsb.pdf"), 
        height = 10, width = 10)
    print(pre_dsb_plt)
    dev.off()
    
    return(blank_pick) 
  } else {
    return("NULL")
  }
  
}
##
pick.blank <- function(df, 
                       log_col = "prot.size.log10", 
                       group_col = "group",
                       min_bk = 1000) {
  x_bk <- df[[log_col]][df[[group_col]] == "bk"]
  # fitting GMM
  gm <- Mclust(x_bk, G = 2, verbose = FALSE)
  means <- gm$parameters$mean
  sds   <- sqrt(gm$parameters$variance$sigmasq)
  low_comp <- which.min(means)
  mu_low <- means[low_comp]
  sd_low <- sds[low_comp]
  
  # junction
  post <- gm$z[, low_comp]
  t_cut <- median(x_bk[abs(post - 0.5) == min(abs(post - 0.5))])
  
  # quantile
  low_idx <- which(post > 0.5)
  q1 <- quantile(x_bk[low_idx], 0.25, na.rm = TRUE)
  q95<- quantile(x_bk[low_idx], 0.95, na.rm = TRUE)
  
  # blank range candidate
  lo <- max(q1, mu_low - 2 * sd_low)
  hi <- min(q95, mu_low + 2 * sd_low)
  
  # lower than cell peak
  hi_adj <- min(hi, t_cut * 0.9)
  
  return(list("suggested_interval" = list("lo" = lo, "hi" = hi_adj),
              "peak_stats" = list("mu_low" = mu_low, "sd_low" = sd_low, "t_cut" = t_cut)))
}
## denoise on ADT data
denoise.fun.ADT <- function(adt_path = NULL,
                            clean_barcode = NULL,
                            blank_barcode = NULL,
                            idx = NULL,
                            add_idx = T,
                            low_prot = 1, 
                            hi_prot = 2.5,
                            out_path = NULL,
                            seed_use = 42,
                            save_rds = T){
  ## load data and create object
  mtx_file <- glue('{adt_path}/featurecounts.mtx')
  feature_file <- glue('{adt_path}/featurecounts.genes.txt')
  barcode_file <- glue('{adt_path}/featurecounts.barcodes.txt')
  
  if (all(file.exists(mtx_file, feature_file, barcode_file))) {
    
    adt_mat_raw <- readMM(mtx_file) %>% t
    dimnames(adt_mat_raw) <- list(fread2(feature_file, header = FALSE)[[1]] %>%
                                    str_split_i(., "-A[0-1]", 1),
                                  fread2(barcode_file, header = FALSE)[[1]] %>%
                                    paste0(., "-1"))
    ## denoise
    # set parameters
    isotype_ctrl <- rownames(adt_mat_raw)[grep('Ctrl', rownames(adt_mat_raw))]
    # qc on empty droplets
    blank_idx_qc <- which(colnames(adt_mat_raw) %in% blank_barcode &
                            colSums(adt_mat_raw) > 10^low_prot &
                            colSums(adt_mat_raw) < 10^hi_prot)
    
    # run
    message(glue("Run denoise on ADT data using {length(blank_idx_qc)} empty droplets!"))
    set.seed(seed_use)
    adt_dsb <- DSBNormalizeProtein(cell_protein_matrix = adt_mat_raw[, clean_barcode],
                                   empty_drop_matrix = adt_mat_raw[, blank_idx_qc],
                                   denoise.counts = T,
                                   use.isotype.control = T,
                                   isotype.control.name.vec = isotype_ctrl,
                                   return.stats = T)
    # create new ADT object
    scADT_dsb <- adt_dsb$dsb_normalized_matrix
    if(add_idx){
      colnames(scADT_dsb) <- paste0(idx, "_", colnames(scADT_dsb))
    }
    saveRDS(adt_dsb$protein_stats, file = paste0(out_path, "/", idx, "_adt_dsbout.rds"))
    
    #
    if (save_rds) {
      saveRDS(scADT_dsb, file = glue("{out_path}/{idx}_adt_dsbmat.rds"))
    }
    return(scADT_dsb)
    
  } else {
    
    return(NULL)
    
  }
  
}

##### QC function for scRNA-seq ######
qc.fun.RNA <- function(h5_file = NULL,
                       ARC = T,
                       run_dbl = T,
                       dbr_sd = 0.015,
                       seed_db = 20240719,
                       min_features = 200,
                       min_ncounts = 1000,
                       max_ncounts = 20000,
                       max_mtpercent = 10,
                       idx = NULL,
                       add_idx = T,
                       n_core = 1,
                       out_path = NULL,
                       save_rds = T,
                       return_seurat = F){
  ## create object
  if (ARC) {
    scRNA_obj <- CreateSeuratObject(counts = Read10X_h5(h5_file)[["Gene Expression"]],
                                    assay = "RNA", 
                                    min.cells = 0, 
                                    min.features = 0,
                                    project = idx)
  } else {
    scRNA_obj <- CreateSeuratObject(counts = Read10X_h5(h5_file),
                                    assay = "RNA", 
                                    min.cells = 0, 
                                    min.features = 0,
                                    project = idx)
  }
  
  ## remove doublets
  if (run_dbl) {
    
    message("Run scDblFinder on GEX data")
    pre_sce <- as.SingleCellExperiment(scRNA_obj)
    set.seed(seed_db)
    pre_sce <- scDblFinder(pre_sce, 
                           verbose = T, 
                           dbr.sd = dbr_sd,
                           BPPARAM = MulticoreParam(n_core))
    tmp <- gc()
    scRNA_obj@meta.data$scDblFinder_RNA <- pre_sce$scDblFinder.class[match(colnames(scRNA_obj),
                                                                           colnames(pre_sce))]
    # scRNA_obj <- subset(scRNA_obj, cells = colnames(pre_sce)[pre_sce$scDblFinder.class == "singlet"])
    
  }
  ## add prop of UMI from mt genes
  scRNA_obj[["percent.mt"]] <- PercentageFeatureSet(scRNA_obj, pattern = "^MT-")
  # Visualize QC metrics as a violin plot
  pre_violin <- VlnPlot(scRNA_obj, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                        slot = "counts", 
                        ncol = 3, 
                        raster = T)
  
  ## remove cells with too low expression or too high mt genes expression
  scRNA_obj$pass_RNA <- scRNA_obj$nFeature_RNA > min_features &
    scRNA_obj$nCount_RNA > min_ncounts &
    scRNA_obj$nCount_RNA < max_ncounts &
    scRNA_obj$percent.mt < max_mtpercent
  scRNA_obj$pass_RNA[is.na(scRNA_obj$pass_RNA)] <- FALSE
  
  # Visualize QC metrics as a violin plot
  post_violin <- VlnPlot(subset(scRNA_obj, pass_RNA),
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                         slot = "counts",
                         ncol = 3,
                         raster = T)
  tiff(glue("{out_path}/{idx}_qc_gex.tiff"), 
       height = 15, width = 15, units = "in", compression = "lzw", res = 300)
  (pre_violin/post_violin) %>% print()
  dev.off()
  
  if(add_idx){
    colnames(scRNA_obj) <- paste0(idx, "_", colnames(scRNA_obj))
  }
  
  if (save_rds) {
    saveRDS(scRNA_obj, file = glue("{out_path}/{idx}_gex_tmp.rds"))
  }
  if (return_seurat) {
    return(scRNA_obj)
  } else {
    return(colnames(scRNA_obj))
  }
  
}

##### QC function for scATAC-seq ######
qc.fun.ATAC <- function(h5_arc_file = NULL,
                        ARC = T,
                        frag_file = NULL,
                        grange_file = NULL,
                        run_dbl = T,
                        dbr_sd = 1,
                        seed_db = 20240719,
                        min_features = 500,
                        min_ncounts = 1000,
                        max_ncounts = 1E5,
                        max_artef_rt = 0.01,
                        min_tss = 4,
                        max_nucleosome = 4,
                        idx = NULL,
                        add_idx = T,
                        n_core = 1,
                        out_path = NULL,
                        save_rds = T,
                        return_seurat = F){
  ## create object
  grange_anno <- readRDS(grange_file)
  if (ARC) {
    scATAC_obj <- CreateChromatinAssay(counts = Read10X_h5(h5_arc_file)[["Peaks"]],
                                       sep = c(":", "-"),
                                       annotation = grange_anno,
                                       min.cells = 0,
                                       min.features = 0) %>%
      CreateSeuratObject(counts = .,
                         assay = "peaks", 
                         project = idx)
  } else {
    scATAC_obj <- CreateChromatinAssay(counts = Read10X_h5(h5_arc_file),
                                       sep = c(":", "-"),
                                       annotation = grange_anno,
                                       min.cells = 0,
                                       min.features = 0) %>%
      CreateSeuratObject(counts = .,
                         assay = "peaks", 
                         project = idx)
  }
  
  ## add prefix and fragment file
  if(add_idx){
    new_cell_names <- paste0(idx, "_", colnames(scATAC_obj))
    names(new_cell_names) <- new_cell_names
    scATAC_obj <- RenameCells(scATAC_obj, add.cell.id = idx)
  }
  ## add prefix and fragment file
  if (!is.null(frag_file)) {
    frag_obj <- CreateFragmentObject(
      path = frag_file,
      cells = colnames(scATAC_obj),
      validate.fragments = T
    )
    Fragments(scATAC_obj[["peaks"]]) <- frag_obj
  }
  ## remove doublets
  if (run_dbl) {
    
    # 1. use scDblFinder
    message("Run scDblFinder on ATAC data")
    pre_sce <- as.SingleCellExperiment(scATAC_obj)
    set.seed(seed_db)
    pre_sce <- scDblFinder(pre_sce, 
                           aggregateFeatures = T, 
                           # dbr.sd = dbr_sd,
                           processing = "normFeatures",
                           BPPARAM = MulticoreParam(n_core),
                           verbose = T)
    tmp <- gc()
    table(pre_sce$scDblFinder.class)
    
    # 2. use amulet
    message("Run Amulet on ATAC data")
    # set excluded regions
    repeats <- GRanges("chr6", IRanges(1000,2000))
    chr_list <- readRDS("/ix1/wchen/xiangyu/Ref_data/genome/grange_Anno/chr_list.rds")
    otherChroms <- GRanges(c(chr_list$exclude_chr, chr_list$exclude_chr_format),
                           IRanges(1L, width = 10^8))
    fragfile <- Fragments(scATAC_obj)[[1]]@path
    res <- amulet(fragfile,
                  regionsToExclude = c(repeats, otherChroms),
                  BPPARAM = MulticoreParam(n_core))
    tmp <- gc()
    
    # 3. combine two methods
    pre_sce$scDblFinder.p <- 1 - pre_sce$scDblFinder.score
    pre_sce$amulet.p <- res[colnames(pre_sce), "p.value"]
    pre_sce$amulet.q <- res[colnames(pre_sce), "q.value"]
    pre_sce$combined.p <- apply(colData(pre_sce)[, c("scDblFinder.p", "amulet.p")],
                                1,
                                function(x){
                                  x[x < 0.001] <- 0.001
                                  aggregation::fisher(x)
                                })
    pre_sce$combined.class <- ifelse(pre_sce$combined.p > 0.05,
                                     "singlet", "doublet")
    scATAC_obj@meta.data$scDblFinder_RNA <- pre_sce$scDblFinder.class[match(colnames(scATAC_obj),
                                                                            colnames(pre_sce))]
    # scATAC_obj <- subset(scATAC_obj, cells = colnames(pre_sce)[pre_sce$combined.class == "singlet"])
    
  }
  ##
  peaks_keep <- seqnames(granges(scATAC_obj)) %in% standardChromosomes(granges(scATAC_obj))
  scATAC_obj <- scATAC_obj[as.vector(peaks_keep), ]
  
  ## add qc metrics
  # compute nucleosome signal score per cell
  scATAC_obj <- NucleosomeSignal(object = scATAC_obj)
  # compute TSS enrichment score per cell
  scATAC_obj <- TSSEnrichment(object = scATAC_obj)
  # add blacklist ratio
  scATAC_obj$blacklist_ratio <- FractionCountsInRegion(scATAC_obj, 
                                                       assay = 'peaks',
                                                       regions = blacklist_hg38_unified)
  # Visualize QC metrics as a violin plot
  pre_violin <- VlnPlot(scATAC_obj,
                        features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio'),
                        pt.size = 0.1,
                        ncol = 4, 
                        raster = T)
  # apply quality control
  scATAC_obj$pass_ATAC <- scATAC_obj$nFeature_peaks > min_features &
    scATAC_obj$nCount_peaks > min_ncounts &
    scATAC_obj$nCount_peaks < max_ncounts &
    scATAC_obj$nucleosome_signal < max_nucleosome &
    scATAC_obj$TSS.enrichment > min_tss &
    scATAC_obj$blacklist_ratio < max_artef_rt
  scATAC_obj$pass_ATAC[is.na(scATAC_obj$pass_ATAC)] <- FALSE
  
  # Visualize QC metrics as a violin plot
  post_violin <- VlnPlot(subset(scATAC_obj, pass_ATAC),
                         features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio'),
                         pt.size = 0.1,
                         ncol = 4,
                         raster = T)
  
  tiff(glue("{out_path}/{idx}_qc_atac.tiff"), 
       height = 15, width = 20, units = "in", compression = "lzw", res = 300)
  (pre_violin/post_violin) %>% print()
  dev.off()
  
  if (save_rds) {
    saveRDS(scATAC_obj, file = glue("{out_path}/{idx}_atac_tmp.rds"))
  }
  if (return_seurat) {
    return(scATAC_obj)
  } else {
    return(colnames(scATAC_obj))
  }
}

