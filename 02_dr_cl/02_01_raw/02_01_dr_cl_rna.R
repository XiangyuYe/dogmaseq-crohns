## module load r/4.5.0
library(Seurat)
library(bigreadr)
library(dplyr)
library(stringr)
library(glue)
library(tibble)
library(harmony)

## set work dir
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)
code_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/code/00_FUNCTION/"
source(glue("{code_path}process/PROCESS_FUN.R"))

##### 1. merge meta data #####
## set sample id used
all_sample <- read.table(file = "Helmsley_sample.txt",
                         header = F, sep = "\t")[,1]
exclude_sample <- c("Duerr_20221104_DOGMAseq-1",
                    "Duerr_20221104_DOGMAseq-2",
                    "Duerr_20221108_DOGMAseq")
use_sample <- setdiff(all_sample, exclude_sample)
## load meta data
meta_df_all <- lapply(all_sample, function(idx){
  print(idx)
  readRDS(glue("02_clean/{idx}_meta.rds"))
}) %>% Reduce("rbind", .) %>% as.data.frame()
rownames(meta_df_all) <- meta_df_all$Batch_barcode
## exclude non-CD samples
meta_df_use <- subset(meta_df_all, Batch %in% use_sample)
#
saveRDS(meta_df_all, file = "03_output/02_clean/meta_all.rds")
saveRDS(meta_df_use, file = "03_output/02_clean/meta_use.rds")

##### 2. add clinical info #####
clinic_df <- readRDS("03_output/02_clean/clinic_df.rds")
meta_use <- readRDS("03_output/02_clean/meta_use.rds")
sel_cov <- c("Age", "Gender", "Race", "Section_specific",
             "med_strategy", "med_strategy_comb",
             "use_Aminosalic", "use_Immunomodu", "use_MTX",
             "use_AntiTNF", "use_AntiIL23", "use_AntiIntegrin")
meta_use <- cbind(meta_use[, !colnames(meta_use) %in% sel_cov],
                  clinic_df[match(meta_use$Sample_ID_exp, clinic_df$Sample_ID_exp),
                            sel_cov])
saveRDS(meta_use, file = "03_output/02_clean/meta_use_clinic_add.rds")

##### 3. load RNA data #####
meta_use_clinic_add <- readRDS("03_output/02_clean/meta_use_clinic_add.rds")
use_sample <- unique(meta_use_clinic_add$Batch)
scRNA_obj_raw <- load.seurat(glue("02_clean/{use_sample}_gex.rds"), 
                         type = "RNA",
                         merge = T)
scRNA_obj_raw@meta.data <- meta_use_clinic_add[colnames(scRNA_obj_raw), ]
saveRDS(scRNA_obj_raw, file = "03_output/02_clean/scRNA_obj_raw.rds")

## set uninterested genes (will be excluded for DRCL)
all_gene <- rownames(scRNA_obj_raw)
gene_del_label <- grepl("^MT-", all_gene) |
  grepl("^RP[1-9]", all_gene) |
  grepl("^RP[SL]", all_gene) |
  grepl("^MIR", all_gene) |
  grepl("^LINC", all_gene) |
  (grepl("^A[A-Z][0-9][0-9]", all_gene) &
     grepl("\\.", all_gene)) |
  (grepl("^B[A-Z][0-9][0-9]", all_gene) &
     grepl("\\.", all_gene)) |
  (grepl("^C[A-Z][0-9][0-9]", all_gene) &
     grepl("\\.", all_gene))
exclude_gene <- all_gene[gene_del_label]
saveRDS(exclude_gene, file = "03_output/02_clean/exclude_gene.rds")

##### 4. DRCL on all clean cells #####
## load data and run SCT+harmony (not recommended to be honest!)
exclude_gene <- readRDS("03_output/02_clean/exclude_gene.rds")
scRNA_obj_raw <- readRDS("03_output/02_clean/scRNA_obj_raw.rds")
scRNA_obj_raw <- subset(scRNA_obj_raw, Group == "Crohn’s_disease")
seed_use <- 20250528
pc_clust_gex <- 20
resx <- 0.6
scRNA_obj <- SCT.harmony(scRNA_obj = scRNA_obj_raw,
                         batch_col = "Batch",
                         n_higvar = 3000,
                         exclude_gene = exclude_gene,
                         seed_use = seed_use)
saveRDS(scRNA_obj, file = "03_output/02_clean/scRNA_obj_SCT_harmony.rds")

## clustering
scRNA_obj <- FindNeighbors(scRNA_obj, 
                           reduction = "harmony_SCT",
                           n.trees = 300,
                           k.param = 20, 
                           dims = 1:pc_clust_gex)
scRNA_obj <- FindClusters(scRNA_obj, 
                          n.iter = 300,
                          resolution = resx,
                          random.seed = seed_use)
scRNA_obj <- RunUMAP(scRNA_obj, 
                     dims = 1:pc_clust_gex, 
                     reduction = "harmony_SCT",
                     n.neighbors = 20L,
                     umap.method = "uwot",
                     n.epochs = 300,
                     negative.sample.rate = 20L,
                     min.dist = 0.3,
                     seed.use = seed_use)
## plot
out_path <- "03_output/03_clustering/raw/"
tiff(file = glue("{out_path}UMAP_RNA.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(scRNA_obj, 
        reduction = "umap", 
        raster = F,
        label = T) + NoLegend()
dev.off()
tiff(file = glue("{out_path}UMAP_RNA_section.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(scRNA_obj, 
        reduction = "umap", 
        group.by = "Section", 
        raster = F,
        label = T) + NoLegend()
dev.off()
tiff(file = glue("{out_path}UMAP_RNA_condition.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(scRNA_obj, 
        reduction = "umap", 
        group.by = "Condition", 
        raster = F,
        label = T) + NoLegend()
dev.off()
##### 5. rough annotation for immune and non-immune cells #####
DefaultAssay(scRNA_obj) <- "RNA"
scRNA_obj <- NormalizeData(scRNA_obj)
feature_list <- list(
  "T_features_list" = c("CD2", "CD3D", "CD3E", "CD3G"),
  "NK_features_list" = c("KLRF1","KLRD1","NCAM1","CD160"),
  "MP_features_list" = c("CD14", "CD163", "CD68", "CSF1R") ,
  "DC_features_list" = c("CD1C", "LILRA4", "CD1E", "FCER1A"),
  "B_features_list" = c("CD79A", "MS4A1", "JSRP1", "MZB1", "CD38"),
  "PLT_features_list" = c("PF4", "PPBP", "GNG11"),
  "NonIM_features_list" = c("PECAM1", "EPCAM", "KRT18", "ACTA2", 
                            "TPSAB1", "TPSAB2", "CPA3",
                            "DCN", "LUM", "COL1A2", "COL1A1"))
dot_ref <- DotPlot(scRNA_obj, 
                   features = unlist(feature_list) %>% unique) + 
  theme(axis.text = element_text(angle = 90))
tiff(glue("{out_path}dot_test.tiff"), 
     height = 6, width = 10, units = "in", res = 300, compression = "lzw")
dot_ref
dev.off()
##
scRNA_obj$ann_raw <- "Immune"
scRNA_obj$ann_raw[scRNA_obj$SCT_snn_res.0.6 == 16] <- "non_Immune"
tiff(file = glue("{out_path}UMAP_RNA_ann.tiff"),
     width = 6, height = 6, units = "in", res = 600, compression = "lzw")
DimPlot(scRNA_obj, 
        group.by = "ann_raw",
        raster = F,
        label = T) + 
  NoLegend()
dev.off()

##### 6. save meta data and scRNA obj for immune cells #####
sc_meta <- scRNA_obj@meta.data %>% rownames_to_column(., var = "cell")
sc_meta <- cbind(sc_meta,
                 Embeddings(scRNA_obj, reduction = "umap") %>% as.data.frame())
saveRDS(sc_meta, file = glue("{out_path}sc_meta_raw.rds"))

##
scRNA_obj_raw$ann_raw <- scRNA_obj@meta.data$ann_raw
saveRDS(subset(scRNA_obj_raw, ann_raw == "Immune"), 
        file = "03_output/02_clean/scRNA_obj_immune.rds")

# ##### 7. clustering for check (PASSED) #####
# scRNA_obj_sub <- subset(scRNA_obj,
#                         cells = rownames(sc_meta)[sc_meta$ann_raw == "Immune"])
# scRNA_obj_sub <- FindNeighbors(scRNA_obj_sub,
#                                reduction = "harmony_SCT",
#                                n.trees = 300,
#                                k.param = 20,
#                                dims = 1:pc_clust_gex)
# scRNA_obj_sub <- FindClusters(scRNA_obj_sub, 
#                               n.iter = 300,
#                               graph.name = "SCT_snn",
#                               resolution = 0.6,
#                               random.seed = seed_use)
# scRNA_obj_sub <- RunUMAP(scRNA_obj_sub,
#                          dims = 1:pc_clust_gex,
#                          reduction = "harmony_SCT",
#                          n.neighbors = 20L,
#                          umap.method = "uwot",
#                          n.epochs = 300,
#                          negative.sample.rate = 20L,
#                          min.dist = 0.3,
#                          seed.use = seed_use)
# ## plot
# tiff(file = glue("{out_path}UMAP_RNA_immune.tiff"),
#      width = 6, height = 6, units = "in", res = 600, compression = "lzw")
# DimPlot(scRNA_obj_sub,
#         reduction = "umap",
#         raster = F,
#         label = T) + NoLegend()
# dev.off()
# dot_ref <- DotPlot(scRNA_obj_sub,
#                    features = unlist(feature_list) %>% unique) +
#   theme(axis.text = element_text(angle = 90))
# tiff(glue("{out_path}dot_test_immune.tiff"),
#      height = 6, width = 10, units = "in", res = 300, compression = "lzw")
# dot_ref
# dev.off()
# ##
# saveRDS(scRNA_obj, file = glue("{out_path}scRNA_obj_test.rds"))
# saveRDS(scRNA_obj@meta.data, file = glue("{out_path}scRNA_meta_test.rds"))
# saveRDS(scRNA_obj[["SCT"]], file = glue("{out_path}scSCT_obj_test.rds"))
# saveRDS(scRNA_obj@reductions, file = glue("{out_path}reducRNA_test.rds"))
# saveRDS(scRNA_obj@commands, file = glue("{out_path}cmdRNA_test.rds"))

