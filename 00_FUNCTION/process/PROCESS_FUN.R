#
library(Seurat)
library(Signac)
library(ArchR)
library(bigreadr)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)
library(harmony)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(patchwork)
##### Process functions #####
##
load.seurat <- function(file_list,
                        type = NULL,
                        common_peaks = T,
                        join_RNA = T,
                        merge = T){
  
  sc_list <- lapply(file_list, function(x){
    sc_objx <- readRDS(x)
  })
  #
  if (type %in% c("ATAC", "scATAC") & common_peaks) {
    combined_peaks <- lapply(sc_list, function(sc_objx){
      granges(sc_objx)
    }) %>% Reduce("c", .) %>% GenomicRanges::reduce()
    #
    peakwidths <- width(combined_peaks)
    combined_peaks <- combined_peaks[peakwidths  < 10000 & peakwidths > 20]
    
    sc_list <- lapply(sc_list, function(sc_objx){
      
      frag_objx <- Fragments(sc_objx)
      sc_objx <- FeatureMatrix(
        fragments = frag_objx,
        features = combined_peaks,
        cells = colnames(sc_objx)) %>%
        CreateChromatinAssay(., 
                             fragments = frag_objx) %>%
        CreateSeuratObject(., 
                           assay = "peaks", 
                           meta.data = sc_objx@meta.data)
      
    })
    
  }
  #
  if (merge) {
    sc_obj <- merge(x = sc_list[[1]], 
                    y = sc_list[2:length(sc_list)],
                    merge.data = TRUE)
    ## join RNA layers
    if (type %in% c("RNA", "scRNA") & join_RNA) {
      DefaultAssay(sc_obj) <- "RNA"
      sc_obj <- JoinLayers(sc_obj)
    }
    return(sc_obj)
    
  } else {
    return(sc_list)
  }
}
##
SCT.harmony <- function(scRNA_obj = NULL,
                        batch_col = "orig.ident",
                        n_higvar = 3000,
                        exclude_gene = NULL,
                        seed_use = 42){
  ## SCT
  scRNA_list <- SplitObject(scRNA_obj, split.by = batch_col) %>%
    lapply(., 
           FUN = SCTransform, 
           method = "glmGamPoi", 
           min_cells = 1,
           return.only.var.genes = F)
  ## select high variable genes
  var_sel_sct <- SelectIntegrationFeatures(object.list = scRNA_list, 
                                           nfeatures = n_higvar + length(exclude_gene)) %>%
    setdiff(., exclude_gene)
  var_sel_sct <- var_sel_sct[1:n_higvar]
  ## merge and join object
  scRNA_obj <- merge(x = scRNA_list[[1]], 
                     y = scRNA_list[2:length(scRNA_list)], 
                     merge.data = T)
  scRNA_obj <- JoinLayers(scRNA_obj, assay = "RNA")
  DefaultAssay(scRNA_obj) <- "SCT"
  VariableFeatures(scRNA_obj) <- var_sel_sct
  scRNA_obj[["SCT"]] <- subset(scRNA_obj[["SCT"]], features = var_sel_sct)
  
  ##
  scRNA_obj <- RunPCA(scRNA_obj, 
                      npcs = 50, 
                      assay = "SCT",
                      verbose = F, 
                      ndims.print = 1:5, 
                      nfeatures.print = 10) %>%
    RunHarmony(.,
               group.by.vars = batch_col, 
               reduction.use = "pca",
               reduction.save = "harmony_SCT",
               assay.use = "SCT",
               project.dim = F)
  
  return(scRNA_obj)
}
# ##
# dr.cl.RNA <- function(scRNA_obj = NULL,
#                       var_genes = NULL,
#                       use_assay = "RNA",
#                       batch_var = "orig.ident",
#                       vars_reg = NULL,
#                       run_harmony = T,
#                       var_exp = 0.8,
#                       res = 0.6,
#                       n_iter = 200,
#                       run_umap = T,
#                       n_neig = 30L,
#                       n_epochs = 300,
#                       neg_rate = 10L,
#                       min_dist = 0.3,
#                       seed_use = 20240816){
#   
#   DefaultAssay(scRNA_obj) <- use_assay
#   scRNA_obj <- NormalizeData(scRNA_obj) 
#   if (!is.null(var_genes)) {
#     VariableFeatures(scRNA_obj) <- var_genes
#   } else {
#     scRNA_obj <- FindVariableFeatures(scRNA_obj)
#   }
#   scRNA_obj <- ScaleData(scRNA_obj, 
#                          vars.to.regress = vars_reg) %>%
#     RunPCA(., npcs = 50)
#   ##
#   if (run_harmony) {
#     use_dr <- "harmony"
#     scRNA_obj <- RunHarmony(scRNA_obj,
#                             group.by.vars = batch_var, 
#                             reduction.use = "pca",
#                             reduction.save = use_dr,
#                             assay.use = use_assay,
#                             project.dim = F)
#   } else {
#     use_dr <- "pca"
#   }
#   eig_val <- (scRNA_obj@reductions[[use_dr]]@stdev)^2
#   var_explained <- eig_val / sum(eig_val)
#   pc_clust_gex <- min(which(cumsum(var_explained) >= var_exp))
#   
#   ##
#   scRNA_obj <- FindNeighbors(scRNA_obj, 
#                              reduction = use_dr,
#                              n.trees = 200,
#                              k.param = 20, 
#                              dims = 1:pc_clust_gex)
#   scRNA_obj <- FindClusters(scRNA_obj, 
#                             n.iter = n_iter,
#                             resolution = res,
#                             random.seed = seed_use)
#   
#   if (run_umap) {
#     scRNA_obj <- RunUMAP(scRNA_obj, 
#                          dims = 1:pc_clust_gex, 
#                          reduction = use_dr,
#                          n.neighbors = n_neig,
#                          umap.method = "uwot",
#                          n.epochs = n_epochs,
#                          negative.sample.rate = neg_rate,
#                          min.dist = min_dist,
#                          seed.use = seed_use)
#   }
#   print(pc_clust_gex)
#   return(scRNA_obj)
# }
# ##
# drcl.harmony <- function(seurat_obj = NULL,
#                          re_dr = T,
#                          n_var_ft = 2000,
#                          exclude_ft = NULL,
#                          use_assay = NULL,
#                          batch_var = NULL,
#                          var_reg = NULL,
#                          pc_clust = 20,
#                          res = 1,
#                          run_UMAP = T,
#                          seed_use = 123L){
#   ##
#   if(re_dr | !any(grepl(glue("harmony_{use_assay}"), names(seurat_obj@reductions)))){
#     #
#     n_var_ft <- min(nrow(seurat_obj[[use_assay]]$counts), n_var_ft)
#     seurat_obj <- NormalizeData(seurat_obj) %>%
#       FindVariableFeatures(nfeatures = n_var_ft + length(exclude_ft))
#     var_sel <- VariableFeatures(seurat_obj, assay = use_assay) %>%
#       setdiff(., exclude_ft)
#     var_sel <- var_sel[1:n_var_ft]
#     VariableFeatures(seurat_obj) <- var_sel
#     #
#     seurat_obj <- ScaleData(seurat_obj, vars.to.regress = var_reg) %>%
#       RunPCA(., 
#              assay = use_assay,
#              reduction.name = glue("pca_{use_assay}"),
#              verbose = T) %>%
#       RunHarmony(.,
#                  group.by.vars = batch_var, 
#                  reduction.use = glue("pca_{use_assay}"),
#                  reduction.save = glue("harmony_{use_assay}"),
#                  assay.use = use_assay,
#                  project.dim = F)
#   }
#   # clustering
#   seurat_obj <- FindNeighbors(seurat_obj, 
#                               reduction = glue("harmony_{use_assay}"),
#                               n.trees = 200,
#                               k.param = 20, 
#                               dims = 1: pc_clust)
#   seurat_obj <- FindClusters(seurat_obj, 
#                              n.iter = 300,
#                              resolution = res,
#                              random.seed = seed_use)
#   if (run_UMAP) {
#     seurat_obj <- RunUMAP(seurat_obj, 
#                           reduction = glue("harmony_{use_assay}"), 
#                           dims = 1: pc_clust, 
#                           reduction.name = glue("umap.harmony_{use_assay}"),
#                           return.model = T, # for project
#                           random.seed = seed_use)
#     
#   }
#   return(seurat_obj)
# }

##
dr.cl.ATAC <- function(scATAC_obj = NULL,
                       cutoff_q = "q5",
                       use_assay = "peaks",
                       run_harmony = T,
                       batch_var = "orig.ident",
                       var_exp = 0.8,
                       res = 0.6,
                       run_umap = T,
                       seed_use = 20240816){
  ##
  DefaultAssay(scATAC_obj) <- use_assay
  scATAC_obj <- RunTFIDF(scATAC_obj,
                         method = 1) %>%
    FindTopFeatures(min.cutoff = cutoff_q) %>%
    RunSVD()
  
  if (run_harmony) {
    use_dr <- "harmony_lsi"
    scATAC_obj <- RunHarmony(scATAC_obj,
                             group.by.vars = batch_var, 
                             reduction.use = "lsi",
                             dims.use = 2:50,
                             reduction.save = use_dr,
                             assay.use = use_assay,
                             project.dim = F)
    eig_val <- (scATAC_obj@reductions[[use_dr]]@stdev)^2
    var_explained <- eig_val / sum(eig_val)
    pc_clust_hlsi <- min(which(cumsum(var_explained) >= var_exp))
  } else {
    use_dr <- "lsi"
    eig_val <- (scATAC_obj@reductions[[use_dr]]@stdev)^2
    var_explained <- eig_val / sum(eig_val[-1])
    pc_clust_hlsi <- min(which(cumsum(var_explained[-1]) >= var_exp)) + 1
  }
  ##
  scATAC_obj <- FindNeighbors(scATAC_obj, 
                              reduction = use_dr,
                              n.trees = 200,
                              k.param = 20, 
                              dims = 1:pc_clust_hlsi)
  scATAC_obj <- FindClusters(scATAC_obj, 
                             n.iter = 200,
                             resolution = res,
                             random.seed = seed_use)
  
  if (run_umap) {
    scATAC_obj <- RunUMAP(scATAC_obj, 
                          dims = 1:pc_clust_hlsi, 
                          reduction = use_dr,
                          n.neighbors = 30L,
                          umap.method = "uwot",
                          n.epochs = 300,
                          negative.sample.rate = 10L,
                          min.dist = 0.3,
                          seed.use = seed_use)
  } 
  print(pc_clust_hlsi)
  return(scATAC_obj)  
}
##
peak.recall <- function(scATAC_obj = NULL,
                        ref_peaks = NULL,
                        group_by = "seurat_clusters",
                        macs2_path = "/ihome/wchen/zhongli/.conda/envs/MACS2/bin/macs2",
                        assay_use = "peaks",
                        new_assay = "recall_peaks",
                        combine_peaks = T,
                        tmp_path = NULL,
                        out_prefix = NULL){
  
  tim1 <-Sys.time()
  # call peaks using MACS2
  DefaultAssay(scATAC_obj) <- assay_use
  if (is.null(ref_peaks)) {
    if(is.null(tmp_path)){
      tmp_path <- tempdir()
    } else {
      system(glue("mkdir -p {tmp_path}"))
    }
    peaks <- CallPeaks(scATAC_obj, 
                       assay = assay_use,
                       group.by = group_by,
                       combine_peaks = T,
                       macs2.path = macs2_path,
                       outdir = tmp_path)
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, 
                                     pruning.mode = "coarse")
    peaks <- subsetByOverlaps(peaks, 
                              ranges = blacklist_hg38_unified, 
                              invert = TRUE)
    # combine recalled peaks with original peaks
    if (combine_peaks) {
      peaks <- GenomicRanges::reduce(c(granges(scATAC_obj), peaks))
      peakwidths <- width(peaks)
      peaks <- peaks[peakwidths  < 10000 & peakwidths > 20]
    }
  } else {
    peaks <- ref_peaks
  }

  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(scATAC_obj),
    features = peaks,
    cells = colnames(scATAC_obj)
  )
  # create a new assay using the MACS2 peak set and add it to the Seurat object
  scATAC_obj[[new_assay]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = Fragments(scATAC_obj),
    annotation = Annotation(scATAC_obj)
  )
  if (!is.null(out_prefix)) {
    saveRDS(peaks, file = paste0(out_prefix, "_peaks.rds"))
    saveRDS(macs2_counts, file = paste0(out_prefix, "_macs2_counts.rds"))
    saveRDS(scATAC_obj, file = paste0(out_prefix, "_scATAC_obj.rds"))
  }
  tim2 <- Sys.time()
  message(tim2 - tim1)
  return(scATAC_obj)
}
##
anno.motif <- function(scATAC_obj = NULL,
                       assay_name = "peaks",
                       species_use = "human"){
  ##
  DefaultAssay(scATAC_obj) <- assay_name
  if (species_use == "human") {
    opt_use <- list(species = 9606, 
                    all_versions = FALSE)
    genome_use <- BSgenome.Hsapiens.UCSC.hg38
    
  } else {
    opt_use <- list(collection = "CORE", 
                    tax_group = "vertebrates", 
                    species = "Mus musculus", 
                    all_versions = FALSE)
    genome_use <- BSgenome.Mmusculus.UCSC.mm10
  }
  
  pwm_set <- getMatrixSet(x = JASPAR2020, 
                          opts = opt_use)
  # add motif information
  scATAC_obj <- AddMotifs(object = scATAC_obj,
                          genome = genome_use,
                          pfm = pwm_set)
  # Note that this step can take 30-60 minutes 
  scATAC_obj <- RunChromVAR(object = scATAC_obj,
                            genome = genome_use)
  return(scATAC_obj)
}

##
dr.cl.wnn <- function(sc_obj = NULL,
                      redc_list = NULL,
                      dim_list = NULL,
                      k_nn = 20,
                      prune_SNN = 1/15,
                      n_iter = 300,
                      res = 0.6,
                      cl_method = 1,
                      run_umap = T,
                      n_neig = 30L,
                      n_epochs = 300,
                      neg_rate = 10L,
                      min_dist = 0.3,
                      sprd = 0.5,
                      seed_use = 42){
  #
  sc_obj <- FindMultiModalNeighbors(sc_obj,
                                    reduction.list = redc_list, 
                                    k.nn = k_nn,
                                    prune.SNN = prune_SNN,
                                    dims.list = dim_list)
  sc_obj <- FindClusters(sc_obj, 
                         n.iter = n_iter,
                         graph.name = "wsnn", 
                         algorithm = cl_method,
                         resolution = res, 
                         verbose = T,
                         random.seed = seed_use)
  if (run_umap) {
    sc_obj <- RunUMAP(sc_obj, 
                      nn.name = "weighted.nn", 
                      reduction.name = "wnn.umap", 
                      reduction.key = "wnnUMAP_",
                      n.neighbors = n_neig,
                      umap.method = "uwot",
                      n.epochs = n_epochs,
                      negative.sample.rate = neg_rate,
                      min.dist = min_dist,
                      spread = sprd,
                      seed.use = seed_use)
  }
  return(sc_obj)
}
##
scsub.renorm <- function(sc_obj = NULL,
                         do_ADT = T,
                         do_ATAC = T,
                         q_atac = "q25",
                         do_RNA = T,
                         do_SCT = F,
                         do_harmony = T,
                         batch_col = "orig.ident"){
  #
  if (do_ADT) {
    DefaultAssay(sc_obj) <- "ADT"
    sc_obj <- ScaleData(sc_obj) %>% 
      RunPCA(assay = "ADT")
    if (do_harmony) {
      sc_obj <- RunHarmony(sc_obj,
                           group.by.vars = batch_col, 
                           reduction.use = "pca",
                           reduction.save = "harmony_adt",
                           assay.use = "ADT",
                           project.dim = F)
    }
    
  }
  if (do_ATAC) {
    assay_atac <- intersect(names(sc_obj@assays), c("ATAC", "peaks", "peak"))[1]
    DefaultAssay(sc_obj) <- assay_atac
    sc_obj <- dr.cl.ATAC(scATAC_obj = sc_obj,
                         cutoff_q = q_atac,
                         use_assay = assay_atac,
                         run_harmony = do_harmony,
                         batch_var = batch_col,
                         run_umap = F,
                         seed_use = seed_use)
  }
  if (do_SCT) {
    DefaultAssay(sc_obj) <- "SCT"
    sc_obj <- RunPCA(sc_obj, assay = "SCT") 
    if (do_harmony) {
      sc_obj <- RunHarmony(sc_obj,
                           group.by.vars = batch_col, 
                           reduction.use = "pca",
                           reduction.save = "harmony_SCT",
                           assay.use = "SCT",
                           project.dim = F)
    }
  }
  if (do_RNA) {
    DefaultAssay(sc_obj) <- "RNA"
    sc_obj <- RunPCA(sc_obj, assay = "SCT") #%>%
    if (do_harmony) {
      sc_obj <- RunHarmony(sc_obj,
                           group.by.vars = batch_col, 
                           reduction.use = "pca",
                           reduction.save = "harmony_SCT",
                           assay.use = "RNA",
                           project.dim = F)
    }
  }
  return(sc_obj)
}

##### Plot functions #####
plt.fun <- function(seurat_obj = NULL,
                    marker_list = NULL,
                    use_assay = NULL,
                    group_by = "seurat_clusters",
                    reduc_use = "umap.harmony",
                    maker_ft = NULL,
                    plt_prefix = "./test"){
  #
  p1 <- DimPlot(seurat_obj,
                reduction = reduc_use,
                group.by = group_by, 
                raster = F, 
                label = T)
  tiff(glue("{plt_prefix}_UMAP_raw.tiff"), 
       height = 6, width = 7, units = "in", res = 300, compression  = "lzw")
  print(p1)
  dev.off()
  
  #
  dot_mk <- DotPlot(seurat_obj,
                    assay = use_assay, 
                    col.min = -2, 
                    col.max= 2,
                    group.by = group_by,
                    features = marker_list) +
    theme(axis.text.x = element_text(angle = 90))
  tiff(glue("{plt_prefix}_dot_mk.tiff"), 
       height = 6, width = 15, units = "in", res = 300, compression = "lzw")
  print(dot_mk)
  dev.off()
  #
  wd_plt <- sqrt(length(maker_ft)) %>% ceiling()
  ht_plt <- (length(maker_ft)/wd_plt) %>% ceiling() 
  dt_mk <- FeaturePlot(seurat_obj, 
                       reduction = reduc_use,
                       features = maker_ft, 
                       ncol = wd_plt,
                       raster = F, 
                       min.cutoff = "q1", 
                       max.cutoff = "q99")
  tiff(glue("{plt_prefix}_ft_mk.tiff"), 
       height = ht_plt *4, width = wd_plt * 4 + 2, units = "in", res = 300, compression = "lzw")
  print(dt_mk)
  dev.off()
  #
  return("Plot saved!")
}


dot.minmax <- function(seurat_obj = NULL,
                       assay = "RNA",
                       features = NULL,
                       col_min = 0,
                       col_max = 1,
                       group = "seurat_clusters",
                       col_use = c("lightgrey", "blue")){
  
  dot_seurat <- Seurat::DotPlot(object = seurat_obj, 
                                assay = assay, 
                                group.by = group,
                                features = features, 
                                scale = F)
  dot_data <- dot_seurat$data
  dot_data <- dot_data %>%
    dplyr::group_by(features.plot) %>%
    dplyr::mutate(avg.exp.scaled = (avg.exp - min(avg.exp)) / (max(avg.exp) - min(avg.exp))) %>%
    dplyr::ungroup()
  # dot_data <- split(dot_data, dot_data$features.plot) %>%
  #   lapply(., function(datax){
  #     datax$avg.exp.scaled <- (datax$avg.exp - min(datax$avg.exp)) / (max(datax$avg.exp) - min(datax$avg.exp))
  #     return(datax)
  #   }) %>% Reduce("rbind", .) %>% as.data.frame()
  dot_data$avg.exp.scaled[dot_data$avg.exp.scaled < col_min] <- 0
  dot_data$avg.exp.scaled[dot_data$avg.exp.scaled > col_max] <- 1
  ##
  dot_plt <- ggplot(dot_data, aes(x = features.plot, y = id)) +
    geom_point(aes(fill = avg.exp.scaled, size = pct.exp), 
               shape = 21, color = "grey30") +
    scale_fill_gradient(name = "Scaled Expression",
                         low = col_use[1], 
                         high = col_use[2],
                         breaks = c(0, 1),
                         labels = c("min", "max")) +
    scale_size_continuous(name = "Percent Expressed",
                          breaks = seq(20, 100, 20), 
                          range = c(0, 6))
  if (!is.null(names(features))) {
    dot_plt <- dot_plt + scale_x_discrete(labels = names(features)) 
    
  }
  dot_plt <- dot_plt + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
  
  return(dot_plt)
} 
##
peak.set.plt <- function(sc_obj = NULL,
                         assay_use = "peaks",
                         mk_list = NULL,
                         hl_list = NULL,
                         extend_kb_up = 1000,
                         extend_kb_down = 1000,
                         mk_strand = "+",
                         extend_type = c("prop", "kb"),
                         color_use = NULL){
  ##
  sc_obj[["peaks"]] <- sc_obj[[assay_use]]
  
  extend_kb_up <- rep(extend_kb_up, times = length(mk_list))[seq_along(mk_list)]
  extend_kb_down <- rep(extend_kb_down, times = length(mk_list))[seq_along(mk_list)]
  mk_strand <- rep(mk_strand, times = length(mk_list))[seq_along(mk_list)]
  ##
  plt_list <- lapply(seq_along(mk_list), function(xx){
    ##
    mkx <- mk_list[xx]
    rgx <- LookupGeneCoords(sc_obj, 
                            mkx,
                            assay = "peaks")
    mk_strandx <- mk_strand[xx]
    if (extend_type == "kb") {
      extend_kb_upx <- min(c(extend_kb_up[xx], rgx@ranges@width))
      extend_kb_downx <- min(c(extend_kb_down[xx], rgx@ranges@width))
    } else {
      extend_kb_upx <- round(extend_kb_up[xx] * rgx@ranges@width)
      extend_kb_downx <- round(extend_kb_down[xx] * rgx@ranges@width)
    }
    if (mk_strandx == "-") {
      rgx <- Extend(rgx, upstream = extend_kb_downx, downstream = extend_kb_upx)
    } else {
      rgx <- Extend(rgx, upstream = extend_kb_upx, downstream = extend_kb_downx)
    }
    
    pks_gr <- granges(sc_obj[["peaks"]])
    pk_hits <- findOverlaps(query = rgx, 
                            subject = pks_gr, 
                            ignore.strand = T)
    if (length(pk_hits) > 0) {
      overlap_peaks <- pks_gr[subjectHits(pk_hits)]
      rgx_union <- union(rgx, overlap_peaks)
    } else {
      rgx_union <- rgx
    }
    ##
    cvgx <- CoveragePlot(
      object = sc_obj, 
      assay = "peaks",
      region = rgx_union,
      annotation = F,
      region.highlight = hl_list[[mkx]],
      peaks = F,
      extend.upstream = 0,
      extend.downstream = 0
    ) + 
      ggtitle(mkx) + 
      theme(axis.text = element_blank(), 
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            strip.placement = "outside")
    if (!is.null(color_use)) {
      cvgx <- cvgx + 
        scale_fill_manual(values = rep(color_use, length = nlevels(Idents(sc_obj))))
    }
    if (mkx != mk_list[1]) {
      cvgx <- cvgx + 
        theme(strip.text.y.left = element_blank())
    } else {
      cvgx <- cvgx + 
        theme(strip.text.y.left = element_text(size = 12, color = "black"))
    }
    return(cvgx)
  })
  return(patchwork::wrap_plots(plt_list, nrow = 1))
}
##
heat.adt <- function(sc_obj = NULL,
                     clus_col = "seurat_clusters",
                     assay_use = "ADT",
                     adt_use = NULL){
  ##
  adt_plot <- as.data.frame(t(sc_obj[[assay_use]]$data))
  adt_plot[adt_plot < 0] <- 0
  adt_plot <- cbind(clusters = sc_obj@meta.data[[clus_col]],
                    adt_plot) 
  adt_plot <- adt_plot %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarize_at(.vars = adt_use, .funs = mean) %>%
    tibble::remove_rownames() %>%
    reshape2::melt()
  
  dot_adt <- ggplot(adt_plot) +
    geom_tile(aes(x = variable, y = clusters, fill = value)) +
    scale_fill_viridis_c(option = "B") +
    # scale_fill_manual(values = viridis::viridis(25, option = "B")) + 
    # scale_fill_gradient(low = "grey90", high = "darkgreen") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12, color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          axis.title = element_blank(),
          panel.grid = element_blank())
  
  # adt_plot = cbind(clusters = sc_obj@meta.data[[clus_col]],
  #                  as.data.frame(t(sc_obj[[assay_use]]$scale.data))) %>%
  #   dplyr::group_by(clusters) %>%
  #   dplyr::summarize_at(.vars = adt_use, .funs = median) %>%
  #   tibble::remove_rownames() %>%
  #   tibble::column_to_rownames("clusters")
  # # plot a heatmap of the average dsb normalized values for each cluster
  # pheatmap::pheatmap(adt_plot,
  #                    scale = "none",
  #                    color = viridis::viridis(20, option = "B"),
  #                    fontsize_row = 8, 
  #                    cluster_rows = F, 
  #                    cluster_cols = F,
  #                    border_color = NA)
  
  return(dot_adt)
}
##
feature.plt <- function(seurat_obj = NULL,
                        axis_x = NULL,
                        axis_y = NULL,
                        assay_use = "RNA",
                        layer_use = "data",
                        col_use = c("lightgrey", "#FF7000", "#FF5000","#FF3000"),
                        feature_list = NULL,
                        ncol = 4){
  # choose marker genes
  DefaultAssay(seurat_obj) <- assay_use
  norm_dat <- FetchData(seurat_obj, 
                        vars = feature_list,
                        layer = layer_use) %>%
    rownames_to_column(., var = "cell")
  # add annotation
  norm_dat$UMAP1 <- axis_x
  norm_dat$UMAP2 <- axis_y
  plot_dat <- reshape2::melt(norm_dat, 
                             id.vars = c("cell", "UMAP1", "UMAP2"), 
                             variable.name = "Gene",
                             value.name = "Expression")
  
  plt_list <- lapply(feature_list, function(genex){
    
    plot_datx <- subset(plot_dat, plot_dat$Gene == genex)
    
    pltx <- ggplot(plot_datx, aes(x = UMAP1, y = UMAP2, color = Expression)) + 
      geom_point(alpha = 1, size = 0.001) + 
      scale_colour_gradientn(colours = col_use)+
      ggtitle(genex) + 
      theme_bw()+
      theme(legend.position = "none",
            plot.title = element_text(size = 15, face = "bold", 
                                      color = "black", hjust = 0.5),
            panel.background = element_rect(fill = "white", colour = "black", size = 1),
            panel.spacing.x = unit(0, "cm"), 
            panel.spacing.y = unit(0, "cm"), 
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank())
  })
  
  plt <- wrap_plots(plt_list, ncol = ncol)
  return(plt)
}
##
vln.plt <- function(seurat_obj = NULL,
                    clus_col = NULL,
                    assay_use = NULL,
                    min_exp = 5,
                    feature_list = NULL,
                    color_use = NULL){
  # choose marker genes
  DefaultAssay(seurat_obj) <- assay_use
  norm_df <- FetchData(seurat_obj, 
                       vars = feature_list, 
                       assay = assay_use)
  # add annotation
  norm_df$celltype <- seurat_obj@meta.data[[clus_col]]
  plot_dat <- norm_df %>%
    rownames_to_column("cell") %>%
    pivot_longer(cols = -c(cell, celltype), 
                 names_to = "Marker", 
                 values_to = "Expression")
  plot_dat$Marker <- factor(plot_dat$Marker, levels = feature_list)
  # correct on dsb norm data
  if (assay_use == "ADT" & min(norm_df[, feature_list]) < 0) {
    plot_dat$Expression[plot_dat$Expression < min_exp] <- 0
  }
  # plot
  # plt <- ggplot(plot_dat, aes(y = Expression, x = celltype))+
  #   geom_violin(aes(fill = celltype), scale = "width", trim = T)+
  #   facet_grid(.~plot_dat$Gene, scales = "free_y")+
  #   scale_fill_manual(values = ct_cols)+
  #   scale_x_discrete("")+
  #   scale_y_continuous("")+
  #   theme_bw()+
  #   theme(panel.background = element_rect(fill = "white", colour = "black", size = 1),
  #         panel.spacing.y = unit(0, "cm"), 
  #         strip.placement = "outside", 
  #         strip.background = element_blank(),
  #         strip.text.y = element_text(size = 12,color = "black", angle = 0),
  #         axis.text.y.left  = element_blank(),
  #         axis.ticks.y.left = element_blank(),
  #         axis.text.x.bottom = element_text(size = 12,color = "black",
  #                                           angle = 90,hjust = 1,vjust = 1),
  #         legend.position = "none")
  plt <- ggplot(plot_dat, aes(celltype, Expression))+
    geom_violin(aes(fill = celltype), scale = "width", trim = T)+
    facet_grid(.~plot_dat$Marker, scales = "free", switch = "x")+
    scale_fill_manual(values = rep(color_use, length = nlevels(plot_dat$celltype))) +
    scale_x_discrete("", limits = rev(levels(plot_dat$celltype)))+
    scale_y_continuous("")+
    theme_bw()+
    coord_flip()+
    theme(panel.spacing.y = unit(0, "cm"), 
          strip.placement = "outside", 
          strip.background = element_blank(),
          strip.text.x = element_text(size = 12, color = "black", angle = 90),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
  
  return(plt)
  
}
##
pt.plt.fun <- function(sc_meta = NULL,
                       ann_col = NULL,
                       wd_group_col = NULL,
                       ht_group_col = NULL,
                       color_use = NULL){
  # set levels
  if (!is.factor(sc_meta[[wd_group_col]])) {
    sc_meta[[wd_group_col]] <- as.factor(sc_meta[[wd_group_col]])
  }
  if (!is.factor(sc_meta[[ht_group_col]])) {
    sc_meta[[ht_group_col]] <- as.factor(sc_meta[[ht_group_col]])
  }
  # table width * height
  sc_meta$group_comb <- paste0(sc_meta[[wd_group_col]], ":", sc_meta[[ht_group_col]])
  tab_prop <- table(sc_meta[[ann_col]], sc_meta$group_comb) %>% 
    prop.table(2) %>% as.data.frame()
  colnames(tab_prop) <- c("Cell type", "group_comb", "Proportion")
  tab_prop$Proportion <- tab_prop$Proportion * 100
  tab_prop$wd_group <- str_split_i(tab_prop$group_comb, ":", 1) %>%
    factor(levels = levels(sc_meta[[wd_group_col]]))
  tab_prop$ht_group <- str_split_i(tab_prop$group_comb, ":", 2) %>%
    factor(levels = levels(sc_meta[[ht_group_col]]))
  #
  pie_plt <- ggplot(tab_prop, aes(x = 3, y = Proportion, fill = `Cell type`)) + 
    geom_col(width = 1.5, color = NA) + 
    facet_grid(ht_group ~ wd_group, switch = "y") + 
    coord_polar(theta = "y") + 
    xlim(c(0.2, 3.8)) + 
    theme_void()+
    theme(strip.text = element_text(size = 15, face = "bold"), 
          legend.title = element_text(size = 15, face = "bold"), 
          legend.text = element_text(size = 12))
  if (!is.null(color_use)) {
    pie_plt <- pie_plt + scale_fill_manual(values = color_use)
  }
  return(pie_plt)
}
##
vln.compare.plt <- function(seurat_obj = NULL,
                            use_assay = "RNA",
                            diff_df = NULL,
                            sel_ct = NULL,
                            sel_ct_nm = NULL,
                            ct_col = NULL,
                            group_col = NULL,
                            sel_gene = NULL,
                            col_use = NULL){
  
  sel_cts <- paste(sel_ct, collapse = ",")
  if (is.null(sel_ct_nm)) {
    sel_ct_nm <- sel_cts
  }
  DefaultAssay(seurat_obj) <- use_assay
  # seurat_obj <- NormalizeData(seurat_obj)
  ## format exp df
  exp_df <- data.frame(Exp = seurat_obj[[use_assay]]$data[sel_gene,],
                       Group = seurat_obj@meta.data[[group_col]])
  exp_df$Celltype <- ifelse(seurat_obj@meta.data[[ct_col]] %in% sel_ct, sel_ct_nm, "Others") %>%
    factor(levels = c(sel_ct_nm, "Others"))
  ## extract test result
  diff_use <- subset(diff_df, gene == sel_gene &
                       cluster == sel_cts)
  diff_use$Group <- diff_use[[group_col]]
  diff_use <- diff_use[match(levels(exp_df$Group), diff_use$Group),]
  diff_use$p_adj_label <- format(diff_use$p_val_adj, scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
  n_comp <- nrow(diff_use)
  # format p value and fold change df for plot
  p_df <- data.frame(y.position = rep(max(exp_df$Exp)*0.98, n_comp),
                     x = diff_use$Group,
                     xmin = 1:n_comp - 0.2,
                     xmax = 1:n_comp + 0.2,
                     group2 = rep(sel_ct_nm, n_comp),
                     group1 = rep("Others", n_comp),
                     padj_sig = paste0("Padj = ", diff_use$p_adj_label))
  fc_df <- data.frame(y.position = rep(max(exp_df$Exp) * 1.05, n_comp),
                      x = diff_use$Group,
                      xmin = 1:n_comp - 0.2,
                      xmax = 1:n_comp + 0.2,
                      group2 = rep(sel_ct_nm, n_comp),
                      group1 = rep("Others", n_comp),
                      fc_sig = paste0("Log2FC = ",
                                      round(diff_use$avg_log2FC, 2)))
  ## plot
  vln_comp_plot <- ggplot(exp_df) + 
    geom_violin(aes(x = Group, 
                    y = Exp, 
                    fill = Celltype), scale = "width", trim = T) +
    scale_y_continuous(limits = c(0, max(exp_df$Exp) * 1.1)) + 
    scale_fill_manual(values = col_use) + 
    labs(x = "", 
         y = glue("Expression of {sel_gene}")) +
    stat_pvalue_manual(fc_df,
                       label = "fc_sig",
                       tip.length = 0,
                       size = 0,
                       label.size = 4,
                       bracket.size = 0) +
    stat_pvalue_manual(p_df,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.title = element_text(size = 15, face = "bold"), 
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12))
  
  return(vln_comp_plot)
}

##
test.marker.plt <- function(sc_obj = NULL,
                            clus_col = "seurat_clusters",
                            reduc_map = "wnn.umap",
                            adt_assay = "ADT",
                            rna_assay = "RNA",
                            motif_assay = "chromvar",
                            atac_assay = "ATAC",
                            adt_mk = NULL,
                            rna_mk = NULL,
                            motif_mk = NULL,
                            atac_mk = NULL,
                            ft_mk = NULL,
                            out_path = "./",
                            prefix = "sc_obj"){
  if (!is.null(adt_mk)) {
    dot_adt <- heat.adt(sc_obj = sc_obj,
                        clus_col = clus_col,
                        assay_use = adt_assay,
                        adt_use = adt_mk)
    tiff(glue("{out_path}{prefix}_dot_adt.tiff"),
         height = 6, width = 4, units = "in", res = 300, compression = "lzw")
    print(dot_adt)
    dev.off()
    
  }
  if (!is.null(rna_mk)) {
    dot_rna <-     DotPlot(sc_obj, 
                           assay = "RNA", 
                           group.by = clus_col,
                           features = rna_mk) + 
      theme(axis.text.x = element_text(angle = 90))
    tiff(glue("{out_path}{prefix}_dot_rna.tiff"), 
         height = 6, width = 12, units = "in", res = 300, compression = "lzw")
print(dot_rna)
    dev.off()
    
  }
  if (!is.null(motif_mk)) {
    if(!is.null(motif_mk)){
      names_motif <- names(motif_mk)
      names(motif_mk) <- NULL
    } else {
      names_motif <- motif_mk
    }
    dot_motif <- DotPlot(sc_obj, 
                         assay = "chromvar", 
                         col.min = -2,
                         col.max = 2,
                         features = motif_mk,
                         group.by = clus_col,
                         cols = c("lightgrey", "brown")) + 
      scale_x_discrete(labels = names_motif) +
      theme(axis.text.x = element_text(angle = 90))
    tiff(glue("{out_path}{prefix}_dot_motif.tiff"), 
         height = 8, width = 5, units = "in", res = 300, compression = "lzw")
    print(dot_motif)
    dev.off()
  }
  if (!is.null(atac_mk)) {
    
  }
  if (!is.null(ft_mk)) {
    ft_plt <- FeaturePlot(sc_obj, 
                          reduction = reduc_map,
                          features = ft_mk, 
                          ncol = 5, min.cutoff = "q5", max.cutoff = "q95")
    tiff(glue("{out_path}{prefix}_ft.tiff"), 
         height = 12, width = 22, units = "in", res = 300, compression = "lzw")
    print(ft_plt)
    dev.off()
    
  }
}
