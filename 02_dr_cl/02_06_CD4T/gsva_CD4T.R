## module load r/4.5.0
library(Seurat)
library(clusterProfiler)
library(GSVA)
library(glue)
library(limma)
library(SummarizedExperiment)
library(BiocParallel)
#
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
##
ct <- "CD4T"
cov_use = c("Age", "Gender", 
            "use_Aminosalic", "use_Immunomodu", "use_MTX",
            "use_AntiTNF", "use_AntiIL23", "use_AntiIntegrin")
block_col = "Sample_exp"
sample_col = "Sample_ID_exp"
ann_col <- "ann_level4_final"
group_col <- "Condition"
subgroup_col <- "Section"

##
out_path <- "03_output/03_clustering/CD4T/WNN_ADT_RNA/"
data_path <- "03_output/03_clustering/"
sc_meta_cd4 <- readRDS(glue("{data_path}{ct}/WNN_ADT_RNA/sc_meta.rds"))
sc_meta_cd4$Age <- scale(sc_meta_cd4$Age, center = T, scale = T) %>% as.numeric()
#
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
scRNA_obj <- subset(scRNA_obj, cells = rownames(sc_meta_cd4))
DefaultAssay(scRNA_obj)<-"RNA"
## 
sc_meta <- sc_meta_cd4[, c(block_col, sample_col, ann_col, 
                           group_col, subgroup_col, cov_use)]
scRNA_obj@meta.data <- sc_meta
exp_bulk <- AggregateExpression(scRNA_obj, 
                                assays = "RNA", 
                                return.seurat = F, 
                                group.by = c("Sample_ID_exp", "ann_level4_final"))[["RNA"]]
bulk_meta <- sc_meta_cd4
bulk_meta$bulk_ID <- paste0(gsub("_", "-", bulk_meta$Sample_ID_exp),
                            "_",
                            gsub("_", "-", bulk_meta$ann_level4_final))
bulk_meta <- bulk_meta[!duplicated(bulk_meta$bulk_ID),]
rownames(bulk_meta) <- bulk_meta$bulk_ID
#
msig_path <- "/ix1/wchen/xiangyu/Ref_data/pathway/MSigDB/Hs.symbols.v2024.1/"
go_gmt <- read.gmt(glue("{msig_path}c5.go.bp.v2024.1.Hs.symbols.gmt"))
reactome_gmt <- read.gmt(glue("{msig_path}c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
kegg_gmt <- read.gmt(glue("{msig_path}c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt"))
imsigdb_gmt <- read.gmt(glue("{msig_path}c7.immunesigdb.v2024.1.Hs.symbols.gmt"))
#
path_gmt <- rbind(go_gmt, reactome_gmt, kegg_gmt, imsigdb_gmt)
path_gmt <- subset(path_gmt, gene %in% rownames(scRNA_obj[["RNA"]]))
path_list <- split(path_gmt$gene, f = path_gmt$term)

## run GSVA
exp_bulk_lcpm <- edgeR::cpm(exp_bulk, log = T, normalized.lib.sizes = T)
gsva_param <- gsvaParam(exprData = exp_bulk_lcpm,
                        geneSets = path_list, 
                        minSize = 20, 
                        maxSize = 500,
                        kcdf = "Gaussian",
                        maxDiff = T)
gsva_df <- gsva(gsva_param,
                BPPARAM = MulticoreParam(workers = 20))
saveRDS(gsva_df, file = glue("{out_path}gsva_df.rds"))
##
gsva_obj <- CreateSeuratObject(gsva_df, assay = "GSVA", meta.data = bulk_meta)
gsva_obj[["GSVA"]]$data <- gsva_obj[["GSVA"]]$counts
gsva_obj <- ScaleData(gsva_obj)
Idents(gsva_obj) <- gsva_obj$ann_level4_final
saveRDS(gsva_obj, file = glue("{out_path}gsva_obj.rds"))




