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
data_path <- "03_output/03_clustering/"
out_path <- "03_output/03_clustering/CD4T/WNN_ADT_RNA/"
sc_meta_CD4 <- readRDS(glue("{data_path}CD4T/WNN_ADT_RNA/sc_meta.rds"))
#
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
scRNA_obj <- subset(scRNA_obj, cells = rownames(sc_meta_CD4)[sc_meta_CD4$ann_level3_final == "CD4_Treg"])
DefaultAssay(scRNA_obj) <- "RNA"
scRNA_obj <- NormalizeData(scRNA_obj)

msig_path <- "/ix1/wchen/xiangyu/Ref_data/pathway/MSigDB/Hs.symbols.v2024.1/"
go_gmt <- read.gmt(glue("{msig_path}c5.go.bp.v2024.1.Hs.symbols.gmt"))
reactome_gmt <- read.gmt(glue("{msig_path}c2.cp.reactome.v2024.1.Hs.symbols.gmt"))
kegg_gmt <- read.gmt(glue("{msig_path}c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt"))
imsigdb_gmt <- read.gmt(glue("{msig_path}c7.immunesigdb.v2024.1.Hs.symbols.gmt"))

path_gmt <- rbind(go_gmt, reactome_gmt, kegg_gmt, imsigdb_gmt)
path_gmt <- subset(path_gmt, gene %in% rownames(scRNA_obj[["RNA"]]))
path_list <- split(path_gmt$gene, f = path_gmt$term)
gsva_param <- gsvaParam(exprData = as.matrix(scRNA_obj[["RNA"]]$data),
                        geneSets = path_list, 
                        minSize = 20, 
                        maxSize = 500,
                        kcdf = "Gaussian",
                        maxDiff = T)
gsva_df <- gsva(gsva_param,
                BPPARAM = MulticoreParam(workers = 40))
saveRDS(gsva_df, file = glue("{out_path}gsva_treg_df.rds"))

#####
gsva_df <- readRDS(glue("{out_path}gsva_treg_df.rds"))
scRNA_obj[["GSVA"]] <- CreateAssayObject(gsva_df)
# scRNA_obj@meta.data <- readRDS(glue("{out_path}sc_meta.rds"))
DefaultAssay(scRNA_obj) <- "GSVA"
Idents(scRNA_obj) <- scRNA_obj$ann_level4_refine
diff_gsva <- FindAllMarkers(scRNA_obj, 
                            only.pos = T)
saveRDS(diff_gsva, "03_output/03_clustering_test/CD4T/test153000/test302015/diff_gsva_treg.rds")

##
top_df <- diff_gsva %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  slice_head(n = 10) %>%
  ungroup()
top_df <- top_df[order(top_df$p_val_adj, -abs(top_df$avg_log2FC)),]
top_df <- top_df[!duplicated(top_df$gene),]
top_df <- top_df[order(top_df$cluster, top_df$p_val_adj, -abs(top_df$avg_log2FC)),]

scRNA_obj_bulk <- Seurat::AverageExpression(scRNA_obj,
                                         features = unique(top_df$gene),
                                         group.by = "ann_level4_refine",
                                         assays = "GSVA",
                                         return.seurat = T,
                                         layer = "data")
##
rownames(scRNA_obj_bulk) <- gsub("\\-", " ", rownames(scRNA_obj_bulk)) %>%
  gsub("KEGG", "KEGG:", .) %>%
  gsub("REACTOME", "REACTOME:", .) %>%
  gsub("GOBP", "GOBP:", .)
ann_level <- gsub("\\-", " ", colnames(scRNA_obj_bulk))
ann_col <- data.frame(Celltype = factor(ann_level,
                                        levels = ann_level),
                      row.names = colnames(scRNA_obj_bulk))
anno_colors <- c("#E61737", "#EC6E73", "#F8AAC0")
names(anno_colors) <- ann_level
##
png(glue("{out_path}gsva_heat_treg.png"),
     height = 5.5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(scRNA_obj_bulk[["GSVA"]]$data,
                   scale = "row",
                   annotation_col = ann_col,
                   annotation_colors = list(Celltype = anno_colors),
                   border_color = NA,
                   fontsize = 8,
                   legend = T,
                   annotation_legend = F,
                   show_colnames = F,
                   cluster_cols = F, 
                   cluster_rows = F)
dev.off()


library(Seurat)
library(clusterProfiler)
library(GSVA)
library(glue)
library(limma)
library(SummarizedExperiment)
library(BiocParallel)
#
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/"
setwd(project_path)
##
data_path <- "03_output/03_clustering_test/CD4T/test153000/"
out_path <- "03_output/03_clustering_test/CD4T/test153000/test302015/"

scRNA_obj <- readRDS(glue("{data_path}scRNA_obj_test.rds"))
gsva_df <- readRDS(glue("{out_path}gsva_treg_df.rds"))
scRNA_obj[["GSVA"]] <- CreateAssayObject(gsva_df)
DefaultAssay(scRNA_obj) <- "GSVA"
##### limma #####
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/sc_DEG.R")
de_method <- "limma"
cov_use = c("Age", "Gender", 
            "use_Aminosalic", "use_Immunomodu", "use_MTX",
            "use_AntiTNF", "use_AntiIL23", "use_AntiIntegrin")
block_col = "sample_exp"
sample_col = "sample_id_exp"
ann_col <- "ann_level4_refine"
group_col <- "condition"
subgroup_col <- "section_comb"

sc_meta_all <- readRDS(glue("03_output/03_clustering_test/CD4T/test153000/test302015/sc_meta_refine0726.rds"))
clinic_df <- readRDS(file = "03_output/02_clean/clinic_df.rds")
sc_meta <- cbind(sc_meta_all[, c(block_col, sample_col, 
                                 ann_col, group_col, subgroup_col)],
                 clinic_df[match(sc_meta_all$sample_id_exp, clinic_df$sample_id_exp), cov_use])
scRNA_obj@meta.data <- sc_meta[colnames(scRNA_obj),]
Idents(scRNA_obj) <- scRNA_obj@meta.data$condition

all_group <- c("I_I", "N_U", "N_N")
comb_group <- combn(all_group, 2)
out_path_it <- "03_output/04_Diff/GSVA/small_intestine/CD4T/"
out_path_cr <- "03_output/04_Diff/GSVA/colorect/CD4T/"

######
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/pseudoBulk.R")
for (ngx in c(2, 3)) {
  
  g1 <- comb_group[1, ngx]
  g2 <- comb_group[2, ngx]
  
  scRNA_obj_sub <- subset(scRNA_obj, condition %in% c(g1, g2))
  scRNA_obj_sub$condition <- factor(scRNA_obj_sub$condition,
                                 levels = c(g2, g1))
  out_prefix <- glue("lmm_{ann_col}_{g1}_vs_{g2}.rds")
  deg_list_it <- pseudo.bulk.deg(seurat_obj = subset(scRNA_obj_sub,
                                                     section_comb == "TI"),
                                 use_assay = "GSVA",
                                 ident_col = ann_col,
                                 group_col = "condition",
                                 cov_col = cov_use,
                                 vif_thre = 3,
                                 sample_col = "sample_id_exp",
                                 block_col = block_col,
                                 paired_only = F,
                                 do_norm = F,
                                 agg_strategy = "average",
                                 min_cell_per_sample = 0,
                                 min_percent = 0,
                                 filter_ByExpr = F,
                                 test_use = "limma",
                                 n_core = 4)
  
  saveRDS(deg_list_it, file = glue("{out_path_it}/{de_method}/{out_prefix}"))
  #
  deg_list_cr <- pseudo.bulk.deg(seurat_obj = subset(scRNA_obj_sub,
                                                     section_comb == "Colon"),
                                 use_assay = "GSVA",
                                 ident_col = ann_col,
                                 group_col = "condition",
                                 cov_col = cov_use,
                                 vif_thre = 3,
                                 sample_col = "sample_id_exp",
                                 block_col = block_col,
                                 paired_only = F,
                                 do_norm = F,
                                 agg_strategy = "average",
                                 min_cell_per_sample = 0,
                                 min_percent = 0,
                                 filter_ByExpr = F,
                                 test_use = "limma",
                                 n_core = 4)
  
  saveRDS(deg_list_cr, file = glue("{out_path_cr}/{de_method}/{out_prefix}"))
}

