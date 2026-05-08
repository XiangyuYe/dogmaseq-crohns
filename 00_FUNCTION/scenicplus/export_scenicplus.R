##
library(Seurat)
library(Signac)
library(bigreadr)
library(tidyverse)
library(stringr)
library(tibble)
library(glue)
library(GenomicRanges)
library(reticulate)
library(DropletUtils)
library(Matrix)


export.scenic <- function(seurat_obj = NULL,
                          assay_rna = "RNA",
                          assay_atac = "ATAC",
                          sel_cells = NULL,
                          outpath = NULL){
  #
  if (is.null(sel_cells)) {
    sel_cells <- colnames(seurat_obj)
  }
  seurat_obj <- subset(seurat_obj, cells = sel_cells)
  
  if (!is.null(assay_rna)) {
    scRNA_count <- seurat_obj[[assay_rna]]$counts %>%
      as.matrix() %>% as(., "CsparseMatrix")
    system(glue("mkdir -p {outpath}"))
    if (file.exists(paste0(outpath, "/gex_count/"))) {
      system(glue("rm -rf {outpath}/gex_count"))
    }
    write10xCounts(
      paste0(outpath, "/gex_count/"),
      scRNA_count,
      gene.type = "Gene Expression"
    )
  }
  ##
  fwrite2(seurat_obj@meta.data, 
          file = paste0(outpath, "/metadata.tsv"), 
          sep = "\t", , row.names = T)
  ##
  if (!is.null(assay_atac)) {
    scATAC_assay <- seurat_obj[[assay_atac]]
    peak_regions <- granges(scATAC_assay) %>% as.data.frame()
    bed_data <- peak_regions[, c("seqnames", "start", "end")]
    colnames(bed_data) <- c("chr", "start", "end")
    fwrite2(bed_data,
            file = paste0(outpath, "/consensus_regions.bed"),
            sep = "\t", , row.names = F, col.names = F)
    ##
    scATAC_count <- scATAC_assay@counts
    rownames(scATAC_count) <- paste0(bed_data[,1], ":",
                                     bed_data[,2], "-",
                                     bed_data[,3])
    ##
    writeMM(as(scATAC_count, "CsparseMatrix"),
            file = paste0(outpath, "/peaks.mtx"))
    fwrite2(data.frame(rownames(scATAC_count)),
            file = paste0(outpath, "/features.tsv"),
            sep = "\t", , row.names = F, col.names = F)
    fwrite2(data.frame(colnames(seurat_obj)),
            file = paste0(outpath, "/barcodes.tsv"),
            sep = "\t", , row.names = F, col.names = F)
    
  }
  return(0)
}