##
library(Seurat)
library(tidyverse)
library(magrittr)
library(glue)
library(monocle)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(corrplot)
library(cols4all)
#
Project_path="/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(Project_path)
## set parameters
ct <- "CD4T"
cluster_col <- "ann_level4_final"
group_col <- "Condition"
sed_use <- 20250528
ct_cols <- c("#F8AAC0",  "#D81B60",  "#F06292")
gp_cols <- c("#0073C2FF","#EFC000FF", "#E41A1C")
n_sample <- 20000

##### load data #####
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
sc_meta_CD4 <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta.rds")
scRNA_obj@meta.data <- sc_meta_CD4

out_path <- "03_output/06_Traj/Monocle2/CD4_Treg/"
system(glue("mkdir -p {out_path}"))

source("code/FUNCTION/monocle2/monocle_fun_modify.R")

##### monocle 2 in Tregs #####
scRNA_treg <- subset(scRNA_obj, ann_level3_final == "CD4_Treg")
## select DEGs as order_genes
DefaultAssay(scRNA_treg) <- "RNA"
scRNA_treg <- NormalizeData(scRNA_treg) %>% ScaleData()
Idents(scRNA_treg) <- scRNA_treg$ann_level4_final

ct_markers <- FindAllMarkers(scRNA_treg, 
                             logfc.threshold = 1,
                             only.pos = T)
top_markers <- ct_markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_head(n = 100) %>%
  ungroup()
order_genes <- unique(top_markers$gene)

## run monocle 2
monocle_obj <- monocle2.fun(seurat_obj = scRNA_treg,
                            order_genes = order_genes,
                            cluster_col = cluster_col,
                            group_col = group_col,
                            n_sample = n_sample,
                            sed_use = sed_use,
                            out_path = out_path)
mono_plt <- monocle.plot(monocle_obj = monocle_obj,
                         cluster_col = cluster_col)
ggsave(glue("{out_path}/traj_monocle2.png"), 
       mono_plt$traj_plt1 | mono_plt$traj_plt2,
       height = 4, width = 11, units = "in", dpi = 300)
saveRDS(monocle_obj, file = glue("{out_path}/cds_monocle2.rds"))
##### plot #####
cds <- readRDS(glue("{out_path}/cds_monocle2.rds"))
cds$Condition <- factor(cds$Condition, 
                        levels = c("NN", "NU", "II"))
#
traj_plt1 <- monocle::plot_cell_trajectory(cds, 
                                           color_by = "State",
                                           cell_size = 1,
                                           show_branch_points = F) +
  theme_bw()
ggsave(glue("{out_path}/traj_monocle_state.png"), 
       traj_plt1,
       height = 4, width = 5, units = "in", dpi = 300)

##
for (xx in 4:5) {
  cds_sub <- cds[, pData(cds)$State %in% c(1, 2, 3, xx)]
  cds_sub <- detectGenes(cds_sub, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(cds_sub), num_cells_expressed >= 10))
  message("Start running differential test!")
  res_sub <- differentialGeneTest(
    cds_sub[expressed_genes, ],
    fullModelFormulaStr = "~ sm.ns(Pseudotime)",
    cores = 10,
    verbose = T
  )
  
  saveRDS(res_sub, glue("{out_path}/deg_pseudo_monocle2_sub{xx}.rds"))
}
####
xx <- 4
cds_sub <- cds[, pData(cds)$State %in% c(1, 2, 3, xx)]
deg_pseudo <- readRDS(glue("{out_path}/deg_pseudo_monocle2_sub{xx}.rds"))
top_gene_pseudo <- dplyr::top_n(deg_pseudo, n = 60, wt = -qval) %>%
  dplyr::top_n(., n = 60, wt = vf_vst_counts_variance)
ann_col <- make_add_annotation_col(cds_sub, 
                                   cols = c(cluster_col, group_col))
colnames(ann_col) <- c("Cell type", "Condition")
ann_col[["Cell type"]] <- factor(ann_col[["Cell type"]],
                                 levels = c("CD4_Treg_naive", "CD4_IKZF2low_Treg", "CD4_Treg"),
                                 labels = c("CD4 Treg naive", "CD4 IKZF2low Treg", "CD4 Treg"))
##
my_ann_colors <- list(
  Condition = c(NN = gp_cols[1], 
                NU = gp_cols[2], 
                II = gp_cols[3]),
  `Cell type` = c(`CD4 Treg naive` = ct_cols[1], 
                  `CD4 IKZF2low Treg` = ct_cols[2], 
                  `CD4 Treg` = ct_cols[3])
)
#

png(glue("{out_path}/ptime_heat_monocle2_sub{xx}.png"),
     height = 6, width = 6, units = "in", res = 300)
plot_pseudotime_heatmap_modify(cds_sub[top_gene_pseudo$gene_short_name,],
                               num_clusters = 2,
                               add_annotation_col = ann_col,
                               annotation_colors = my_ann_colors,
                               show_rownames = T,
                               return_heatmap = F)
dev.off()


##
gene_sel <- c("IKZF2", "CCR6", "RORA", 
              "IL23R", "IL26", "MAF",
              "PPARG", "KLRB1", "ADAM19")
xx=6
cds_sub <- cds[gene_sel, pData(cds)$State %in% c(1, 2, 3, xx, 7)]
plt_gene <- plot_genes_in_pseudotime(cds_sub[gene_sel,], 
                                     cell_size = 0.75,
                                     panel_order = gene_sel,
                                     ncol = 3,
                                     color_by = "ann_level4_final") + 
  scale_color_manual(values = ct_cols) +
  theme(legend.title = element_blank())
ggsave(glue("{out_path}/plt_gene.png"), 
       plt_gene,
       height = 6, width = 8, units = "in", dpi = 300)

plt_jitter <- plot_genes_jitter(cds_sub[gene_sel, 
                                        pData(cds_sub)$ann_level4_final == "CD4 IKZF2low Treg"],
                                cell_size = 0.3,
                                grouping = "condition", 
                                color_by = "ann_level4_final", 
                                panel_order = gene_sel,
                                plot_trend = F) +
  stat_summary(aes_string(color = "ann_level2_final"),
               fun.data = "mean_cl_boot", size = 0.35) +
  stat_summary(aes_string(x = "condition", 
                          y = "expression", 
                          color = "ann_level2_final", 
                          group = "ann_level2_final"), 
               fun.data = "mean_cl_boot", 
               size = 0.35, geom = "line") + 
  scale_color_manual(values = c("CD4 IKZF2low Treg" = ct_cols[2],
                                "CD4T" = "grey40")) +
  facet_wrap( ~ feature_label, scales= "free_y", ncol = 3) + 
  theme(legend.position = "none")
ggsave(glue("{out_path}/plt_jitter.png"), 
       plt_jitter,
       height = 6, width = 6, units = "in", dpi = 300)

