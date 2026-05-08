##
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(reshape2)
library(cols4all)
library(glue)
library(tibble)
library(ggpubr)
library(rstatix)

###data input and select parameters
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/"
setwd(project_path)
# set parameters
source("code/FUNCTION/process/PROCESS_FUN.R")
feature_list <- list(
  "T_features_list" = c("CD2", "CD3G"),
  "CD4_features_list" = c("CD4"),
  "CD8_features_list" = c("CD8A","CD8B"),
  "NK_features_list" = c("KLRF1","KLRD1"),
  "B_features_list" = c("CD79A", "MS4A1", "MZB1"),
  "MP_features_list" = c("CD14", "CD163") ,
  "DC_features_list" = c("FCER1A", "CD1C")
)
pair_list <- list("II_vs_NN_TI" = c("II:TI", "NN:TI"),
                  "NU_vs_NN_TI" = c("NU:TI", "NN:TI"),
                  "II_vs_NN_Colon" = c("II:Colon", "NN:Colon"),
                  "NU_vs_NN_Colon" = c("NU:Colon", "NN:Colon"),
                  "TI_vs_Colon_NN" = c("NN:TI", "NN:Colon"))
seed_use <- 20250528

# set path
data_path = "03_output/03_clustering/"
plt_path = "04_plot/Fig1/"

# load data
sc_obj <- readRDS("03_output/03_clustering/all/WNN_ADT_RNA/scWNN_obj.rds")
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)
sc_meta_all <- readRDS("03_output/03_clustering/sc_meta_ann_comb_all1205.rds")
sc_meta_all$Condition <- factor(sc_meta_all$Condition, levels = c("NN", "NU", "II"))
sc_meta_all$Section <- factor(sc_meta_all$Section, levels = c("Colon", "TI"), 
                              labels = c("Colon", "Ileum"))

sc_obj@meta.data <- sc_meta_all
Idents(sc_obj) <- sc_obj$ann_level2_refine

# set plot parameters
ct_cols <- c("#098476", "#03ADF0", "#016FAD", "#E61737", "#F8AAC0",
             "#FED037", "#CFC554", "#F57C33", "#C02020", "#904A30")
names(ct_cols) <- levels(sc_meta_all$ann_level2_refine)

mini_theme <- theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(size = 15, face = "bold", 
                                  color = "black", hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        panel.spacing.x = unit(0, "cm"), 
        panel.spacing.y = unit(0, "cm"), 
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

##### Fig. 1B: UMAP #####
umap_plt <- DimPlot(sc_obj, 
                    reduction = "wnn.umap",
                    cols = ct_cols,
                    label = T, 
                    label.size = 8,
                    raster = F)+
  xlab("wnnUMAP1") + ylab("wnnUMAP2")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#
ggsave(file = glue("{plt_path}UmapPlot_all.png"),
       umap_plt,
       width = 10, height = 10,units = "in", dpi = 300)

##### Fig. 2C plot markers #####
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/process_ATAC.R")
adt_ft <- c("CD3", "CD4", "CD8", "CD19", "HLA-DR","CD64", "CD123")
rna_ft <- c("CD3E", "CD8A", "TRDC", "CD79A", "MS4A1", "MZB1", "KLRF1", "KIT", "CD163", "CD1C", "CLEC4C")

dot_adt_all <- vln.plt(seurat_obj = sc_obj,
                       clus_col = "ann_level2_refine",
                       assay_use = "ADT",
                       feature_list = adt_ft,
                       color_use = ct_cols) + 
  ggtitle("ADT markers") + 
  theme(plot.title = element_text(size = 11, face = "bold"))
dot_rna_all <- dot.minmax(seurat_obj = sc_obj,
                          assay = "RNA",
                          features = rna_ft,
                          col_min = 0,
                          col_max = 1,
                          group = "ann_level2_refine") + 
  scale_fill_distiller(type = 'div', palette = 'RdYlBu')+
  scale_y_discrete(limits = levels(sc_obj$ann_level2_refine) %>% rev)  + 
  ggtitle("RNA markers") + 
  theme(legend.position = "none",
        axis.text.y =  element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 11, face = "bold"))
dot_plt <- (dot_adt_all | dot_rna_all) + 
  plot_layout(widths = c(length(adt_ft) + 2, 
                         length(rna_ft)))

ggsave(file = glue("{plt_path}dot_plt_all.png"),
       dot_plt,
       width = 8, height = 5,units = "in", dpi = 300)

##### Fig. 1D: pie plot #####
ann_col <-  "ann_level2_refine"
group_col <- "Condition"
section_col <- "Section"
##
sc_meta_all$group_comb <- paste0(sc_meta_all[[group_col]], ":", sc_meta_all[[section_col]])
tab_prop <- table(sc_meta_all[[ann_col]], sc_meta_all$group_comb) %>% 
  prop.table(2) %>% as.data.frame()
colnames(tab_prop) <- c("Cell type", "group_comb", "Proportion")
tab_prop$Proportion <- tab_prop$Proportion * 100
tab_prop$Section <- str_split_i(tab_prop$group_comb, ":", 2)
tab_prop$Group <- str_split_i(tab_prop$group_comb, ":", 1)

tab_prop$Group <- factor(tab_prop$Group, 
                         levels = c("NN", "NU", "II"))
##
pie_plt <- ggplot(tab_prop, aes(x = 3, y = Proportion, fill = `Cell type`)) + 
  geom_col(width = 1.5, color = NA) + 
  facet_grid(Group ~ Section, switch = "y") + 
  coord_polar(theta = "y") + 
  xlim(c(0.2, 3.8)) + 
  scale_fill_manual(values = ct_cols) + 
  theme_void()+
  theme(strip.text = element_text(size = 15, face = "bold"), 
        legend.title = element_text(size = 15, face = "bold"), 
        legend.text = element_text(size = 12))
ggsave(file = glue("{plt_path}pie_plt_all.png"),
       pie_plt,
       width = 6, height = 6,units = "in", dpi = 300)

##### Fig. 1E, G and Extended Fig. 1: miloR #####
source("code/FUNCTION/miloR.R")
pair_sample_set <- readRDS("03_output/03_clustering/pair_sample_set_1205.rds")
sc_obj$group_pair <- pair_sample_set$group_pair[match(sc_obj$Sample_ID_exp,
                                                      pair_sample_set$Sample_ID_exp)]
## compare between Conditions
all_pair_set <- c("II vs NN in Ileum", "NU vs NN in Ileum",
                  "II vs NN in Colon", "NU vs NN in Colon")
names(all_pair_set) <- c("TI_II_NN", "TI_NU_NN",
                         "Colon_II_NN", "Colon_NU_NN")
for (pair_setx in names(all_pair_set)) {
  
  sel_cellx <- colnames(sc_obj)[grep(pair_setx, sc_obj@meta.data$group_pair)]
  sc_obj_sub <- subset(sc_obj, 
                       cells = sel_cellx) %>%
    DietSeurat(., 
               graphs = "wsnn", 
               dimreducs = c("wnn.umap", "harmony_lsi", "harmony_SCT"))
  bulk_meta_sub <- sc_obj_sub@meta.data
  bulk_meta_sub <- bulk_meta_sub[!duplicated(bulk_meta_sub$Sample_ID_exp),]
  pair_sample_sub <- bulk_meta_sub$Sample_exp[duplicated(bulk_meta_sub$Sample_exp)]
  sc_obj_sub <- subset(sc_obj_sub, 
                       Sample_exp %in% pair_sample_sub)
  sc_obj_sub$Condition <- droplevels(sc_obj_sub$Condition)
  sc_obj_sub$ann_level2_refine <- factor(sc_obj_sub$ann_level2_refine,
                                         levels = levels(sc_obj_sub$ann_level2_refine),
                                         labels = gsub("_", " ",
                                                       levels(sc_obj_sub$ann_level2_refine)))
  #
  res_milo <- milor.seurat(seurat_obj = sc_obj_sub,
                           n_dims = 30,
                           k_nn = 50,
                           group_col = "Condition",
                           cov_col = "Sample_exp",
                           block_col = NULL,
                           sample_col = "Sample_ID_exp",
                           ann_col = "ann_level2_refine",
                           use_reduc = "harmony_SCT",
                           plot_reduc = "wnn.umap",
                           seed_use = seed_use,
                           n_core = 10)
  res_milo$plt_milo <- res_milo$plt_milo + 
    ggtitle(all_pair_set[pair_setx]) +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  ##
  plt_milo_comb <- (res_milo$plt_milo | res_milo$plt_milo2) + 
    plot_layout(widths = c(6, 4))
  ggsave(glue("{plt_path}traj_milo_{pair_setx}.png"),
         plt_milo_comb,
         height = 7, width = 14, units = "in", dpi = 300)
}

## compare between Sections
pair_setx <- "TI_Colon_NN"
sel_cellx <- colnames(sc_obj)[grep(pair_setx, sc_obj@meta.data$group_pair)]
sc_obj_sub <- subset(sc_obj, 
                     cells = sel_cellx) %>%
  DietSeurat(., 
             graphs = "wsnn", 
             dimreducs = c("wnn.umap", "harmony_lsi", "harmony_SCT"))
bulk_meta_sub <- sc_obj_sub@meta.data
bulk_meta_sub <- bulk_meta_sub[!duplicated(bulk_meta_sub$Sample_ID_exp),]
pair_sample_sub <- bulk_meta_sub$Sample_exp[duplicated(bulk_meta_sub$Sample_exp)]
sc_obj_sub <- subset(sc_obj_sub, 
                     Sample_exp %in% pair_sample_sub)
sc_obj_sub$Section <- factor(sc_obj_sub$Section, levels = c("TI", "Colon"))
sc_obj_sub$ann_level2_refine <- factor(sc_obj_sub$ann_level2_refine,
                                       levels = levels(sc_obj_sub$ann_level2_refine),
                                       labels = gsub("_", " ",
                                                     levels(sc_obj_sub$ann_level2_refine)))
#
res_milo <- milor.seurat(seurat_obj = sc_obj_sub,
                         n_dims = 30,
                         k_nn = 50,
                         group_col = "Section",
                         cov_col = "Sample_exp",
                         block_col = NULL,
                         sample_col = "Sample_ID_exp",
                         ann_col = "ann_level2_refine",
                         use_reduc = "harmony_SCT",
                         plot_reduc = "wnn.umap",
                         seed_use = seed_use,
                         n_core = 10)
#
res_milo$plt_milo <- res_milo$plt_milo + 
  ggtitle("NN in Ileum vs Colon") + 
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
##
plt_milo_comb <- (res_milo$plt_milo | res_milo$plt_milo2) + 
  plot_layout(widths = c(6, 4))
ggsave(glue("{plt_path}traj_milo_{pair_setx}.png"),
       plt_milo_comb,
       height = 7, width = 14, units = "in", dpi = 300)

##### Fig. 1F and Extended Fig. 1: Diff prop ######
sample_col <- "Sample_ID_exp"
block_col <- "Sample_exp"
ann_col <-  "ann_level2_refine"
group_col <- "Condition"
section_col <- "Section"

## 0. format proportion table
bulk_meta <- sc_meta_all[!duplicated(sc_meta_all[[sample_col]]),]
bulk_meta$group_comb  <- paste0(bulk_meta[[group_col]], ":", bulk_meta[[section_col]])
sample_list <- split(bulk_meta$Sample_exp, f = bulk_meta$group_comb)
#
tab_prop <- table(sc_meta_all[[ann_col]], sc_meta_all[[sample_col]]) %>% 
  prop.table(2) %>% as.data.frame()
colnames(tab_prop) <- c("cluster", "sample", "Proportion")
tab_prop$Proportion <- tab_prop$Proportion * 100
tab_prop$block <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                            block_col]
tab_prop$group <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                            group_col]
tab_prop$Section <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                              section_col] %>%
  factor(., levels = c("TI", "Colon"),
         labels = c("Ileum", "Colon"))
tab_prop$group <- factor(tab_prop$group, 
                         levels = c("NN", "NU", "II"))

## 1. paired test between Conditions
all_ct <- unique(tab_prop$cluster)
ppair_prop_df1 <- lapply(1:4, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  tab_prop_use <- subset(tab_prop, block %in% sample_pair_use)
  tab_prop_use$group <- droplevels(tab_prop_use$group)
  group_level <- levels(tab_prop_use$group)
  ppairx <- lapply(all_ct, function(ctx){
    #
    tab_prop_usex <- subset(tab_prop_use, cluster == ctx)
    pair_dfx <- split(tab_prop_usex, f = tab_prop_usex$block) %>%
      lapply(., function(x){
        propx <- x$Proportion[match(group_level, x$group)]
      }) %>% Reduce("rbind", .)
    wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
    
  }) %>% unlist()
  names(ppairx) <- c(all_ct)
  return(ppairx)
  
}) %>% Reduce("cbind", .) %>% as.data.frame()

# adjust p value for each set of comparison
ppair_prop_adj_df1 <- apply(ppair_prop_df1, 2, function(x){
  p.adjust(x, method = "BH")
}) %>% as.data.frame()
colnames(ppair_prop_df1) <- colnames(ppair_prop_adj_df1) <- names(pair_list)[1:4]

## 2.1 plot for Ileum
lab_facet <- c(1, 6)
box_pair_prop_list1_TI <- lapply(all_ct, function(ctx){
  #
  plt_dfx <- subset(tab_prop, cluster == ctx & Section == "Ileum")
  plt_dfx$Condition <- plt_dfx$group
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df1[ctx,2:1], scientific = T, digits = 3) %>%
                   gsub("e", "E", .) %>% as.vector()
  #
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1, max_propx * 1.2),
                         x = c(1.5, 2.5),
                         xmin = c(1, 1),
                         xmax = c(2, 3),
                         group2 = c("NN", "NN"),
                         group1 = c("NU", "II"),
                         padj_sig = padj)
  ylabx <- ifelse(ctx %in% all_ct[lab_facet], "Proportion (%) in Ileum", "")
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Condition, y = Proportion, color = Condition)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(ctx) +
    xlab("") + ylab(ylabx) + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))

  
})

ggsave(glue("{plt_path}box_pair_level2_Ileum.png"),
       patchwork::wrap_plots(box_pair_prop_list1_TI, nrow = 2),
       height = 8, width = 12, units = "in", dpi = 300)

## 2.2 plot for Colon
lab_facet <- c(1, 6)
box_pair_list1_Colon <- lapply(all_ct, function(ctx){
  #
  plt_dfx <- subset(tab_prop, cluster == ctx & Section == "Colon")
  plt_dfx$Condition <- plt_dfx$group
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df1[ctx, 4:3], scientific = T, digits = 3) %>%
                   gsub("e", "E", .) %>% as.vector()
  #
  ylabx <- ifelse(ctx %in% all_ct[lab_facet], "Proportion (%) in Colon", "")
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1, max_propx * 1.2),
                         x = c(1.5, 2.5),
                         xmin = c(1, 1),
                         xmax = c(2, 3),
                         group2 = c("NN", "NN"),
                         group1 = c("NU", "II"),
                         padj_sig = padj)
  
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Condition, y = Proportion, color = Condition)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(ctx) +
    xlab("") + ylab(ylabx) + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})

ggsave(glue("{plt_path}box_pair_level2_Colon.png"),
       patchwork::wrap_plots(box_pair_list1_Colon, nrow = 2),
       height = 8, width = 12, units = "in", dpi = 300)

## 3. paired test between Sections
ppair_prop_df2 <- lapply(5, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  tab_prop_use <- subset(tab_prop, block %in% sample_pair_use)
  tab_prop_use$Section <- droplevels(tab_prop_use$Section)
  group_level <- levels(tab_prop_use$Section)
  ppairx <- lapply(all_ct, function(ctx){
    #
    tab_prop_usex <- subset(tab_prop_use, cluster == ctx)
    pair_dfx <- split(tab_prop_usex, f = tab_prop_usex$block) %>%
      lapply(., function(x){
        propx <- x$Proportion[match(group_level, x$Section)]
      }) %>% Reduce("rbind", .)
    wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
    
  }) %>% unlist()
  names(ppairx) <- c(all_ct)
  return(ppairx)
  
}) %>% Reduce("cbind", .) %>% as.data.frame()

# adjust p value for each set of comparison
ppair_prop_adj_df2 <- apply(ppair_prop_df2, 2, function(x){
  p.adjust(x, method = "BH")
}) %>% as.data.frame()
colnames(ppair_prop_df2) <- colnames(ppair_prop_adj_df2) <- names(pair_list)[5]

## 3.1 plot
lab_facet <- c(1, 6)
box_pair_list2 <- lapply(all_ct, function(ctx){
  #
  plt_dfx <- subset(tab_prop, cluster == ctx & group == "NN")
  plt_dfx$Condition <- plt_dfx$Section
  max_propx <-max(plt_dfx$Proportion)
  #
  padj <- format(ppair_prop_adj_df2[ctx, 1], scientific = T, digits = 3) %>%
                   gsub("e", "E", .) %>% as.vector()
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1),
                         x = c(1.5),
                         xmin = c(1),
                         xmax = c(2),
                         group2 = c("Ileum"),
                         group1 = c("Colon"),
                         padj_sig = padj)
  
  ylabx <- ifelse(ctx %in% all_ct[lab_facet], "Proportion (%) in NN", "")
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Section, y = Proportion, color = Section)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("#0073C2FF", "#E41A1C")) + 
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(ctx) +
    xlab("") + ylab(ylabx) + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})

ggsave(glue("{plt_path}box_pair_level2_NN.png"),
       patchwork::wrap_plots(box_pair_list2, nrow = 2),
       height = 7, width = 10, units = "in", dpi = 300)

##### Fig. 1H inflammatory scores #####
## 1. load NI score
sample_col <- "Sample_ID_exp"
block_col <- "Sample_exp"
group_col <- "Condition"
section_col <- "Section"
#
bulk_meta <- sc_meta_all[!duplicated(sc_meta_all[[sample_col]]),]
bulk_meta$group_comb  <- paste0(bulk_meta[[group_col]], ":", bulk_meta[[section_col]])
rownames(bulk_meta) <- gsub("_", "-", bulk_meta[[sample_col]])
bulk_meta$subgroup <- factor(bulk_meta[[section_col]], 
                             levels = c("TI", "Colon"),
                             labels = c("Ileum", "Colon"))
bulk_meta$group <- factor(bulk_meta[[group_col]], 
                          levels = c("NN", "NU", "II"))
bulk_meta$block <- bulk_meta[[block_col]]
sample_list <- split(bulk_meta$block, f = bulk_meta$group_comb)
# load NI inflammation score
gsva_df_NI <- readRDS("ssgsva_NI_score.rds")
bulk_meta$NI_score <- gsva_df_NI[rownames(bulk_meta), 1] %>% scales::rescale(., c(0, 10))

## 2. test across Conditions
ppair_df <- lapply(1:4, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  bulk_metax <- subset(bulk_meta, block %in% sample_pair_use)
  bulk_metax$group <- droplevels(bulk_metax$group)
  group_level <- levels(bulk_metax$group)
  
  pair_dfx <- split(bulk_metax, f = bulk_metax$block) %>%
    lapply(., function(x){
      propx <- x$NI_score[match(group_level, x$group)]
    }) %>% Reduce("rbind", .)
  ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
  
  return(ppairx)
  
}) %>% unlist()
padjpair_df <- p.adjust(ppair_df, method = "BH")
names(ppair_df) <- names(padjpair_df) <- names(pair_list)[1:4]

## 3. plot
padj <- format(padjpair_df, scientific = T, digits = 3) %>%
                     gsub("e", "E", .) %>% as.vector()
diff_df <- data.frame(y.position = c(10, 9, 10, 9),
                      x = c(1.5, 2.5, 1.5, 2.5),
                      xmin = c(2, 2, 2, 2),
                      xmax = c(1, 3, 1, 3),
                      group2 = c("NN", "NN", "NN", "NN"),
                      group1 = c("II", "NU", "II", "NU"),
                      subgroup = c("Ileum", "Ileum", "Colon", "Colon"),
                      padj_sig = padj)
if_box <- ggplot(bulk_meta, aes(x = group, y = NI_score, color = group)) +
  geom_boxplot(fill = NA, width = 0.5) +
  geom_point(aes(color = group), position = position_jitter(width = 0.1)) +
  geom_line(aes(group = Sample_exp), color = "gray", alpha = 0.5) + 
  scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
  xlim(c("II", "NN", "NU")) + 
  stat_pvalue_manual(diff_df,
                     label = "padj_sig",
                     tip.length = 0.01,
                     label.size = 4,
                     bracket.size = 1) +
  facet_wrap(~ subgroup, scales = 'fixed', nrow = 1) +
  xlab("") + ylab("Inflammation score") + 
  theme_bw() + 
  theme(legend.position = "none",
        title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold"))
ggsave(glue("{plt_path}if_score_box.png"),
       if_box,
       height = 6, width = 6, units = "in", dpi = 300)

##### Fig. 1I: scDRS #####
scDRS_out_path <- glue("03_output/01_scDRS/04_downstream/cov/")
##
scDRS_test_comb <- fread2(glue("{scDRS_out_path}CD.scdrs_group.ann_level2_refine"))
scDRS_test_comb$section <- "Combined"
scDRS_test_comb$celltype <- factor(scDRS_test_comb$group, 
                                   levels = levels(sc_meta_all$ann_level2_refine) %>% rev,
                                   labels = gsub("_", " ", levels(sc_meta_all$ann_level2_refine) %>% rev))
##
scDRS_test <- fread2(glue("{scDRS_out_path}CD.scdrs_group.sample_ct_level2"))
scDRS_test$section <- str_split_i(scDRS_test$group, ":", 1)
scDRS_test$celltype <- str_split_i(scDRS_test$group, ":", 2) %>% 
  factor(., levels = levels(sc_meta_all$ann_level2_refine) %>% rev,
         labels = gsub("_", " ", levels(sc_meta_all$ann_level2_refine) %>% rev))
## 
scDRS_test_plt <- rbind(scDRS_test_comb, scDRS_test)
scDRS_test_plt$section <- factor(scDRS_test_plt$section, 
                                 levels = c("Combined", "TI", "RC"),
                                 labels = c("Combined", "Ileum", "Colon"))
scDRS_test_plt <- split(scDRS_test_plt, f = scDRS_test_plt$section) %>%
  lapply(., function(scDRS_test_pltx){
    scDRS_test_pltx$FDR <- p.adjust(scDRS_test_pltx$assoc_mcp, method = "BH")
    scDRS_test_pltx$hetero_padj <- p.adjust(scDRS_test_pltx$hetero_mcp, method = "BH")
    return(scDRS_test_pltx)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
scDRS_test_plt$assoc_sig <- scDRS_test_plt$FDR < 0.05
scDRS_test_plt$hetero_sig <- scDRS_test_plt$hetero_padj < 0.05

##
scdrs_plt <- ggplot() + 
  geom_tile(data = scDRS_test_plt, 
            aes(x = section, y = celltype, fill = -log10(FDR))) + 
  scale_fill_viridis_c() + 
  geom_tile(data = subset(scDRS_test_plt, assoc_sig), 
            aes(x = section, y = celltype),
            color = "black", fill = NA, size = 1, show.legend = T) +
  geom_text(data = subset(scDRS_test_plt, hetero_sig),
            aes(x = section, y = celltype, label = "\u2716"),
            color = "black", size = 6, show.legend = T) + 
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10))

#
ggsave(file = glue("{plt_path}CD_sample_ct_level2.png"),
       scdrs_plt,
       width = 5, height = 7,units = "in", dpi = 300, limitsize = F)
##### Supp Fig. 1 #####
### T cells ###
sc_obj_t <- readRDS("03_output/03_clustering/T/WNN_ADT_RNA/scWNN_obj.rds")
adt_ft_t <- c("CD3", "CD4", "CD8")
rna_ft_t <- c("CD3E", "TRDC")
ann_col_t <- "ann_level2_final"
ann_level_t <- c("CD4T", "CD8T", "gdT")
sc_obj_t@meta.data[[ann_col_t]] <- factor(sc_obj_t@meta.data[[ann_col_t]],
                                          levels = ann_level_t)
## umap
umap_plt_t <- DimPlot(sc_obj_t, 
                      cols = ct_cols[ann_level_t],
                      reduction = "wnn.umap",
                      label = T, 
                      label.size = 8,
                      raster = F)+
  xlab("wnnUMAP1") + ylab("wnnUMAP2")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank())
## feature plot
DefaultAssay(sc_obj_t) <- "ADT"
ft_plt_adt_t <- FeaturePlot(sc_obj_t, 
                            reduction = "wnn.umap",
                            features = adt_ft_t, 
                            cols = c("lightgrey", "darkgreen"),
                            min.cutoff = "q5", 
                            max.cutoff = "q95",
                            raster = F) %>%
  lapply(., function(x) x + mini_theme)
DefaultAssay(sc_obj_t) <- "RNA"
ft_plt_rna_t <- FeaturePlot(sc_obj_t, 
                            reduction = "wnn.umap",
                            features = rna_ft_t, 
                            cols = c("lightgrey", "#082C8C"),
                            min.cutoff = "q1", 
                            max.cutoff = "q99",
                            raster = F) %>%
  lapply(., function(x) x + mini_theme)
ft_plt_t <- wrap_plots(c(list(umap_plt_t), ft_plt_adt_t, ft_plt_rna_t), ncol = 3)
#
ggsave(file = glue("{plt_path}ft_plt_t.png"),
       ft_plt_t,
       width = 12, height = 8,units = "in", dpi = 300, limitsize = T)

## pie plot
sc_obj_t$Section <- factor(sc_obj_t$Section, 
                           levels = c("Colon", "TI"), 
                           labels = c("Colon", "Ileum"))
pie_plt_t <- pt.plt.fun(sc_meta = sc_obj_t@meta.data,
                        ann_col = ann_col_t,
                        wd_group_col = "Section",
                        ht_group_col = "Condition",
                        color_use = ct_cols[ann_level_t])
ggsave(file = glue("{plt_path}pie_plt_t.png"),
       pie_plt_t,
       width = 6, height = 6,units = "in", dpi = 300, limitsize = T)

### B cells ###
sc_obj_b <- readRDS("03_output/03_clustering/B/WNN_ADT_RNA/scWNN_obj.rds")
adt_ft_b <- c("CD19", "CD319", "CD20")
rna_ft_b <- c("MS4A1", "MZB1")
ann_col_b <- "ann_level2_refine"
ann_level_b <- c("B", "Plasma")
#
Idents(sc_obj_b) <- sc_obj_b@meta.data[[ann_col_b]]
umap_plt_b <- DimPlot(sc_obj_b, 
                      cols = ct_cols[ann_level_b],
                      reduction = "wnn.umap",
                      label = T, 
                      label.size = 8,
                      raster = F)+
  xlab("wnnUMAP1") + ylab("wnnUMAP2")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank())
## plt2 feature plot
DefaultAssay(sc_obj_b) <- "ADT"
ft_plt_adt_b <- FeaturePlot(sc_obj_b, 
                            reduction = "wnn.umap",
                            features = adt_ft_b, 
                            cols = c("lightgrey", "darkgreen"),
                            min.cutoff = "q5", 
                            max.cutoff = "q95",
                            raster = F) %>%
  lapply(., function(x) x + mini_theme)
DefaultAssay(sc_obj_b) <- "RNA"
ft_plt_rna_b <- FeaturePlot(sc_obj_b, 
                            reduction = "wnn.umap",
                            features = rna_ft_b, 
                            cols = c("lightgrey", "#082C8C"),
                            min.cutoff = "q1", 
                            max.cutoff = "q99",
                            raster = F) %>%
  lapply(., function(x) x + mini_theme)
ft_plt_b <- wrap_plots(c(list(umap_plt_b), ft_plt_rna_b), ncol = 3)
#
ggsave(file = glue("{plt_path}ft_plt_b.png"),
       ft_plt_b,
       width = 12, height = 4,units = "in", dpi = 300, limitsize = T)

## pie plot
sc_obj_b$Section <- factor(sc_obj_b$Section, 
                           levels = c("Colon", "TI"), 
                           labels = c("Colon", "Ileum"))
pie_plt_b <- pt.plt.fun(sc_meta = sc_obj_b@meta.data,
                        ann_col = ann_col_b,
                        wd_group_col = "Section",
                        ht_group_col = "Condition",
                        color_use = ct_cols[ann_level_b])
ggsave(file = glue("{plt_path}pie_plt_b.png"),
       pie_plt_b,
       width = 6, height = 6,units = "in", dpi = 300, limitsize = T)

### Naive immune cells ###
sc_obj_n <- readRDS("03_output/03_clustering/all/WNN_ADT_RNA/sc_obj_sub_ILC.rds")
adt_ft_n <- c("CD64", "CD123")
rna_ft_n <- c("CD1C", "KIT", "KLRF1")
ann_col_n <- "ann_level2_final"
ann_level_n <- c("NK", "ILC", "Macrophage", "cDC", "pDC")
sc_obj_n@meta.data[[ann_col_n]] <- factor(sc_obj_n@meta.data[[ann_col_n]],
                                          levels = ann_level_n)
## umap
Idents(sc_obj_n) <- sc_obj_n@meta.data[[ann_col_n]] 
umap_plt_n <- DimPlot(sc_obj_n, 
                      reduction = "wnn.umap",
                      cols = ct_cols[ann_level_n],
                      label = T, 
                      label.size = 8,
                      raster = F)+
  xlab("wnnUMAP1") + ylab("wnnUMAP2")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank())
## feature plot
DefaultAssay(sc_obj_n) <- "ADT"
ft_plt_adt_n <- FeaturePlot(sc_obj_n, 
                            reduction = "wnn.umap",
                            features = adt_ft_n, 
                            cols = c("lightgrey", "darkgreen"),
                            min.cutoff = "q5", 
                            max.cutoff = "q95",
                            raster = F) %>%
  lapply(., function(x) x + mini_theme)
DefaultAssay(sc_obj_n) <- "RNA"
ft_plt_rna_n <- FeaturePlot(sc_obj_n, 
                            reduction = "wnn.umap",
                            features = rna_ft_n, 
                            cols = c("lightgrey", "#082C8C"),
                            min.cutoff = "q1", 
                            max.cutoff = "q99",
                            raster = F) %>%
  lapply(., function(x) x + mini_theme)
ft_plt_n <- wrap_plots(c(list(umap_plt_n), ft_plt_adt_n, ft_plt_rna_n), ncol = 3)
#
ggsave(file = glue("{plt_path}ft_plt_n.png"),
       ft_plt_n,
       width = 12, height = 8,units = "in", dpi = 300, limitsize = T)
## pie plot
sc_obj_n$Section <- factor(sc_obj_n$Section, 
                           levels = c("Colon", "TI"), 
                           labels = c("Colon", "Ileum"))
pie_plt_n <- pt.plt.fun(sc_meta = sc_obj_n@meta.data,
                        ann_col = ann_col_n,
                        wd_group_col = "Section",
                        ht_group_col = "Condition",
                        color_use = ct_cols[ann_level_n])
ggsave(file = glue("{plt_path}pie_plt_n.png"),
       pie_plt_n,
       width = 6, height = 6,units = "in", dpi = 300, limitsize = T)

