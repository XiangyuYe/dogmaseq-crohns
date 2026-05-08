## module load r/4.5.0
library(Seurat)
library(reticulate)
library(DropletUtils)
library(Matrix)
library(dplyr)
library(glue)
library(bigreadr)

## set path
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
##
export.sc <- function(sc_obj = NULL,
                      RNA_assay = "RNA",
                      ADT_assay = "ADT",
                      ATAC_assay = "peaks",
                      sel_cells = NULL,
                      outpath = NULL){
  
  sc_obj <- subset(sc_obj, cells = sel_cells)
  if (!is.null(RNA_assay)) {
    scRNA_count <- sc_obj[[RNA_assay]]$counts %>%
      as.matrix() %>% as(., "CsparseMatrix")
    system(glue("mkdir -p {outpath}"))
    if (file.exists(glue("{outpath}/RNA/"))) {
      system(glue("rm -rf {outpath}/RNA"))
    }
    write10xCounts(
      paste0(outpath, "/RNA/"),
      scRNA_count,
      gene.type = "Gene Expression"
    )
    ##
    fwrite2(sc_obj@meta.data, 
            file = paste0(outpath, "/metadata.tsv"), 
            row.names = T)
  }
  ##
  if (!is.null(ATAC_assay)) {
    scATAC_assay <- sc_obj[[ATAC_assay]]
    peak_regions <- granges(scATAC_assay) %>% as.data.frame()
    bed_data <- peak_regions[, c("seqnames", "start", "end")]
    scATAC_count <- scATAC_assay@counts
    rownames(scATAC_count) <- paste0(bed_data[,1], ":",
                                     bed_data[,2], "-",
                                     bed_data[,3])
    ##
    system(glue("mkdir -p {outpath}{ATAC_assay}/"))
    writeMM(as(scATAC_count, "CsparseMatrix"),
            file = glue("{outpath}/{ATAC_assay}/matrix.mtx"))
    fwrite2(data.frame(rownames(scATAC_count)),
            file = glue("{outpath}/{ATAC_assay}/features.tsv"),
            sep = "\t", , row.names = F, col.names = F)
    fwrite2(data.frame(colnames(sc_obj)),
            file = glue("{outpath}/{ATAC_assay}/barcodes.tsv"),
            sep = "\t", , row.names = F, col.names = F)
  }
  ##
  if (!is.null(ADT_assay)) {
    system(glue("mkdir -p {outpath}{ADT_assay}"))
    scADT_count <- sc_obj[[ADT_assay]]@counts
    writeMM(as(scADT_count, "CsparseMatrix"),
            file = glue("{outpath}/{ADT_assay}/matrix.mtx"))
    fwrite2(data.frame(rownames(scADT_count)),
            file = glue("{outpath}/{ADT_assay}/features.tsv"),
            sep = "\t", , row.names = F, col.names = F)
    fwrite2(data.frame(colnames(sc_obj)),
            file = glue("{outpath}/{ADT_assay}/barcodes.tsv"),
            sep = "\t", , row.names = F, col.names = F)
  }
  
  return(0)
}
##
clean_path <- "03_output/02_clean/"
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
sc_meta <- readRDS("03_output/03_clustering/sc_meta_ann_comb_all1205.rds")
scRNA_obj <- subset(scRNA_obj, cells = rownames(sc_meta))
scRNA_obj@meta.data <- sc_meta
## all
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta),
          outpath = glue("{clean_path}all/"))
## all CD4T
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$ann_level2_refine == "CD4T"],
          outpath = glue("{clean_path}all/CD4T/"))
## all CD8T
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$ann_level2_refine == "CD8T"],
          outpath = glue("{clean_path}all/CD8T/"))
## all B
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$ann_level2_refine == "B"],
          outpath = glue("{clean_path}all/B/"))
#########
## all TI
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "TI"],
          outpath = glue("{clean_path}TI/all/"))
## all TI CD4T
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "TI" &
                                          sc_meta$ann_level2_refine == "CD4T"],
          outpath = glue("{clean_path}TI/CD4T/"))
## all TI CD4CD103TRM17
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "TI" &
                                          sc_meta$ann_level4_final == "CD4_CD103_TRM_Th17"],
          outpath = glue("{clean_path}TI/CD4_CD103_TRM_Th17/"))
## all TI CD8T
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "TI" &
                                          sc_meta$ann_level2_refine == "CD8T"],
          outpath = glue("{clean_path}TI/CD8T/"))
## all TI B
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "TI" &
                                          sc_meta$ann_level2_refine %in% c("B", "Plasma")],
          outpath = glue("{clean_path}TI/B/"))
########
## all Colon
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "Colon"],
          outpath = glue("{clean_path}Colon/all/"))
## all Colon CD4T
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "Colon" &
                                          sc_meta$ann_level2_refine == "CD4T"],
          outpath = glue("{clean_path}Colon/CD4T/"))
## all Colon CD4CD103TRM17
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "Colon" &
                                          sc_meta$ann_level4_final == "CD4_CD103_TRM_Th17"],
          outpath = glue("{clean_path}Colon/CD4_CD103_TRM_Th17/"))
## all Colon CD8T
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "Colon" &
                                          sc_meta$ann_level2_refine == "CD8T"],
          outpath = glue("{clean_path}Colon/CD8T/"))
## all Colon B
export.sc(sc_obj = scRNA_obj,
          RNA_assay = "RNA",
          ADT_assay = NULL,
          ATAC_assay = NULL,
          sel_cells = rownames(sc_meta)[sc_meta$Section == "Colon" &
                                          sc_meta$ann_level2_refine %in% c("B", "Plasma")],
          outpath = glue("{clean_path}Colon/B/"))
