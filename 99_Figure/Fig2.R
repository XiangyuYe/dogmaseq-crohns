## module load r/4.5.0
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(reshape2)
library(cols4all)
library(tibble)
library(ggpubr)
library(rstatix)
###data input and select parameters
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/"
setwd(project_path)

# set parameters
source("code/FUNCTION/process/PROCESS_FUN.R")
pair_list <- list("II_vs_NN_TI" = c("II:Ileum", "NN:Ileum"),
                  "NU_vs_NN_TI" = c("NU:Ileum", "NN:Ileum"),
                  "II_vs_NN_Colon" = c("II:Colon", "NN:Colon"),
                  "NU_vs_NN_Colon" = c("NU:Colon", "NN:Colon"),
                  "TI_vs_Colon_NN" = c("NN:Ileum", "NN:Colon"))
ct_cols <- c("#00897B", "#06C3F2", "#016FAD", "#004F78", "#00A0E9",
             "#FFEB3B", "#FFB300", "#F57C33", "#FF5722", "#D32F2F",  "#A14141",  "#5D4037",
             "#F8AAC0", "#D81B60", "#F06292")
seed_use <- 20250528

# set path
plt_path = "04_plot/Fig2/"

##### load data #####
sc_obj <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/scWNN_obj_drcl.rds")
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
sc_obj[["RNA"]] <- subset(scRNA_obj, cells = colnames(sc_obj))[["RNA"]]
scATAC_obj <- readRDS("03_output/03_clustering/recall_comb_ann_level2_refine_addmotif.rds")
sc_obj[["peaks"]] <- subset(scATAC_obj, cells = colnames(sc_obj))[["peaks"]]

DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj)
sc_meta_CD4 <- readRDS("03_output/03_clustering/CD4T/WNN_ADT_RNA/sc_meta.rds")
sc_meta_CD4$Condition <- factor(sc_meta_CD4$Condition, levels = c("NN", "NU", "II"))
sc_meta_CD4$Section <- factor(sc_meta_CD4$Section, levels = c("Colon", "TI"), 
                              labels = c("Colon", "Ileum"))
sc_obj@meta.data <- sc_meta_CD4
Idents(sc_obj) <- sc_obj$ann_level4_final <- sc_meta_CD4$ann_level4_final <- 
  factor(sc_meta_CD4$ann_level4_final,
         levels = levels(sc_meta_CD4$ann_level4_final),
         labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final)))

##### Fig. 2A: Circle Manha plot #####
library(GenomicRanges)
library(CMplot)

project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

## 0.1 diff path
diff_path <- "03_output/04_Diff/"
deg_path_it <- glue("{diff_path}DEG/TI/all/DESeq2/")
dar_path_it <- glue("{diff_path}DAR/TI/all/DESeq2/")
p2g_path_it <- glue("{diff_path}peak2gene/TI/")

## 0.2 SNP position GRCh 38
summ_df <- fread2("/ix1/wchen/xiangyu/Ref_data/GWAS/Raw/28067908-GCST004132-EFO_0000384.h.tsv.gz")
summ_dff <- subset(summ_df, 
                   hm_chrom %in% c(1:22) &
                     hm_other_allele %in% c("T", "G", "C", "A") &
                     hm_effect_allele %in% c("T", "G", "C", "A") &
                     standard_error > 0,
                   select = c("hm_chrom", "hm_pos", "hm_rsid", "p_value"))
colnames(summ_dff) <- c("chr", "pos", "rsid", "pval_gwas")
gr_gwas <- GRanges(seqnames = paste0("chr", summ_dff$chr), 
                   ranges = IRanges(start = summ_dff$pos, 
                                    width = 1))
## 0.3 Gene coord
gene_coord <- fread2("/ix1/wchen/xiangyu/Ref_data/genome/cellranger/arc_GRCh38/genes/gene_coordinates.tsv")
gene_coord$chr <- gene_coord$seqnames %>% 
  gsub("chr", "", .) %>% as.integer()

## 1. DEG
deg_CD4 <- readRDS(glue("{deg_path_it}pseudo_bulk_ann_level2_refine_II_vs_NN.rds"))[["CD4T"]]
deg_CD4$chr <- gene_coord$seqnames[match(deg_CD4$Term, gene_coord$gene_name)] %>%
  gsub("chr", "", .) %>% as.integer()
# map to nearest SNP pos
gr_gene <- GRanges(seqnames = gene_coord$seqnames,
                   ranges = IRanges(start = gene_coord$start,
                                    end = gene_coord$end))
idx_gene <- IRanges::nearest(gr_gene, gr_gwas)
gene_coord$pos <- summ_dff$pos[idx_gene]
gene_coord$rsid_nearest <- summ_dff$rsid[idx_gene]
deg_CD4$rsid_nearest <- gene_coord$rsid_nearest[match(deg_CD4$Term, gene_coord$gene_name)]

## 2. MAGMA
magma_df <- fread2("03_output/01_scDRS/01_magma/CD.genes.out")
magma_zdf <- fread2("03_output/01_scDRS/01_magma/zscore_CD.tsv")
magma_df$Term <- magma_zdf$GENE
deg_CD4$p_magma <- magma_df$P[match(deg_CD4$Term, magma_df$Term)]

## 3. DAR
dar_CD4 <- readRDS(glue("{dar_path_it}pseudo_bulk_ann_level2_refine_II_vs_NN.rds"))[["CD4T"]]
# map to nearest SNP pos
gr_peak <- StringToGRanges(dar_CD4$Term)
idx_peak <- IRanges::nearest(gr_peak, gr_gwas)
dar_CD4$pos <- summ_dff$pos[idx_peak]
dar_CD4$rsid_nearest <- summ_dff$rsid[idx_peak]
#
dar_CD4 <- subset(dar_CD4, !is.na(padj) & !is.na(pos))
dar_CD4 <- dar_CD4[order(dar_CD4$padj),]
dar_CD4 <- dar_CD4[!duplicated(dar_CD4$rsid_nearest),]

## 4. combine DEG, DAR, and MAGMA
summ_dff$Gene <- deg_CD4$Term[match(summ_dff$rsid, deg_CD4$rsid_nearest)]
summ_dff$log2FC_gene <- deg_CD4$log2FC[match(summ_dff$rsid, deg_CD4$rsid_nearest)]
summ_dff$padj_gene <- deg_CD4$padj[match(summ_dff$rsid, deg_CD4$rsid_nearest)]
summ_dff$p_magma <- deg_CD4$p_magma[match(summ_dff$rsid, deg_CD4$rsid_nearest)]

summ_dff$Peak <- dar_CD4$Term[match(summ_dff$rsid, dar_CD4$rsid_nearest)]
summ_dff$log2FC_peak <- dar_CD4$log2FC[match(summ_dff$rsid, dar_CD4$rsid_nearest)]
summ_dff$padj_peak <- dar_CD4$padj[match(summ_dff$rsid, dar_CD4$rsid_nearest)]

plt_df <- data.frame(SNP = summ_dff$rsid,
                     Chromosome = summ_dff$chr,
                     Position = summ_dff$pos,
                     MAGMA = summ_dff$p_magma,
                     DAP = summ_dff$padj_peak,
                     DEG = summ_dff$padj_gene)
plt_df <- subset(plt_df, 
                 !is.na(MAGMA) | !is.na(DAP) | !is.na(DEG))

## 5. parameters for plot
# dot color for three modalities
col_mat <- matrix(c("#0073C2FF", "#EFC000FF", # MAGMA
                    "royalblue4", "darksalmon", # DAP
                    "grey", "lightblue"), # DEG
                  nrow = 3, byrow = T)

# labels genes: DEGs & MAGMA Sig
top_gene_magma <- deg_CD4$Term[which(deg_CD4$padj < 0.05 & 
                                       abs(deg_CD4$log2FC) > 0.25 &
                                       deg_CD4$p_magma < 1E-5)]
top_gene_magma_pos <- deg_CD4$rsid_nearest[match(top_gene_magma, deg_CD4$Term)]
gene_high_col_magma <- ifelse(deg_CD4$log2FC[match(top_gene_magma, deg_CD4$Term)] > 0,
                              "#F02720", "#2C69B0")

# other highlight genes: DEGs & link to DAPs
links_gr <- readRDS(glue("{p2g_path_it}LinkPeaks_links_CD4T_DEG_II_vs_NN.rds"))
df_links <- as.data.frame(links_gr)
df_links$padj <- p.adjust(df_links$pvalue, method = "BH")
dff_links <- subset(df_links, df_links$padj < 0.05 & df_links$score > 0.1)
dff_links_sig <- subset(dff_links, peak %in% 
                          dar_CD4$Term[which(dar_CD4$padj < 0.05 & abs(dar_CD4$log2FC) > 0.25)] &
                          gene %in% 
                          deg_CD4$Term[which(deg_CD4$padj < 0.05 & abs(deg_CD4$log2FC) > 0.25)])
top_gene_link <- setdiff(unique(dff_links_sig$gene), top_gene_magma)
top_gene_link_pos <- deg_CD4$rsid_nearest[match(top_gene_link, deg_CD4$Term)]
gene_high_col_link <- rep("orange", length(top_gene_link))

# highlight peaks: DEGs & link to DAPs
top_peak_link <- unique(dff_links_sig$peak)
top_peak_link_pos <- dar_CD4$rsid_nearest[match(top_peak_link, dar_CD4$Term)]

# label & bold genes: DEGs & link to DAPs & MAGMA Sig
top_gene_link_magma <- intersect(unique(dff_links_sig$gene), top_gene_magma)

## 6. Circle plot
CMplot(plt_df,
       type = "p",
       plot.type = "c",
       r = 1,
       col = col_mat,
       threshold = c(0.05, 1E-5),
       cir.chr.h = 1.5,
       amplify = F,
       signal.line = 0.2,
       threshold.col = c("blue", "red"),
       threshold.lty = 2,   
       highlight = list("MAGMA" = top_gene_magma_pos,
                        "DAP" = top_gene_link_pos,
                        "DEG" = c(top_gene_magma_pos, top_gene_link_pos)),
       highlight.col = list("MAGMA" = "brown",
                            "DAP" = "orange",
                            "DEG" = c(gene_high_col_magma, gene_high_col_link)),
       highlight.cex = 1,
       highlight.text = list("MAGMA" = NA,
                             "DAP" = NA,
                             "DEG" = c(gene_high_col_magma, gene_high_col_link)),
       highlight.text.cex = 2,
       highlight.text.col = "black",
       outward = T,
       file = "png",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 15,
       height = 15)
system(glue("mv -f Cir_Manhtn.MAGMA_DAP_DEG.png {plt_path}"))

## 99. manha plot for gene modality only for manual labeling ref
CMplot(plt_df[,c(1:3, 6)],
       type = "p",
       plot.type = "m",
       r = 1,
       col = c("grey20", "grey80"),
       threshold = c(0.05, 1E-5),
       cir.chr.h = 1.5,
       amplify = F,
       signal.line = 0.2,
       threshold.col = c("blue", "red"),
       threshold.lty = 2,   
       highlight = list("DEG" = c(top_gene_magma_pos)),
       highlight.col = list("DEG" = c(gene_high_col_magma)),
       highlight.cex = 1,
       highlight.text = list("DEG" = c(top_gene_magma)),
       highlight.text.cex = 2,
       highlight.text.col = "black",
       outward = T,
       file = "png",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 15,
       height = 10)
system(glue("mv -f Rect_Manhtn.DEG.png {plt_path}"))
print(top_gene_link_magma)

##### Fig. 2B and Extended Fig. 2: Vene DEGs #####
use_de <- "DESeq2"
ann_col <- "ann_level2_refine"
deg_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{deg_path}TI/all/{use_de}/")
deg_path_cr <- glue("{deg_path}Colon/all/{use_de}/")
ctx <- "CD4T"

## 
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
deg_df_it2 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[[ctx]]
deg_df_it3 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NU.rds"))[[ctx]]
deg_list_it <- list("II vs NN in Ileum" = deg_df_it1$Term[which(deg_df_it1$padj < 0.05 & abs(deg_df_it1$log2FC) > 0.25)], 
                    "NU vs NN in Ileum" = deg_df_it2$Term[which(deg_df_it2$padj < 0.05 & abs(deg_df_it2$log2FC) > 0.25)],
                    "II vs NU in Ileum" = deg_df_it3$Term[which(deg_df_it3$padj < 0.05 & abs(deg_df_it3$log2FC) > 0.25)])
val_idx <- lapply(deg_list_it, function(x){length(x) > 0}) %>% unlist

vd <- VennDiagram::venn.diagram(
  deg_list_it[val_idx],
  filename = NULL, 
  fill=c("#0073C2FF","#EFC000FF", "#E41A1C")[val_idx], 
  cex  = 1.5, 
  cat.cex = 1.8,
  height = 6, 
  width = 6, 
  units = "in",
  resolution = 300
)

pdf(glue("{plt_path}/venn_{ann_col}_{ctx}_TI2.pdf"),
    height = 8, width = 8)
grid::grid.draw(vd)
dev.off()

##
deg_df_cr1 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
deg_df_cr2 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[[ctx]]
deg_df_cr3 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NU.rds"))[[ctx]]
deg_list_cr <- list("II vs NN in Colon" = deg_df_cr1$Term[which(deg_df_cr1$padj < 0.05 & abs(deg_df_cr1$log2FC) > 0.25)], 
                    "NU vs NN in Colon" = deg_df_cr2$Term[which(deg_df_cr2$padj < 0.05 & abs(deg_df_cr2$log2FC) > 0.25)],
                    "II vs NU in Colon" = deg_df_cr3$Term[which(deg_df_cr3$padj < 0.05 & abs(deg_df_cr3$log2FC) > 0.25)])
val_idx <- lapply(deg_list_it, function(x){length(x) > 0}) %>% unlist

vd <- VennDiagram::venn.diagram(
  deg_list_cr[val_idx],
  filename = NULL, 
  fill=c("#0073C2FF","#EFC000FF", "#E41A1C")[val_idx], 
  cex  = 1.5, 
  cat.cex = 1.8,
  height = 6, 
  width = 6, 
  units = "in",
  resolution = 300
)

pdf(glue("{plt_path}/venn_{ann_col}_{ctx}_Colon.pdf"),
    height = 6, width = 6)
grid::grid.draw(vd)
dev.off()

##
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
deg_df_cr1 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
vd <- VennDiagram::venn.diagram(
  list("II vs NN in ileum" = deg_df_it1$Term[which(deg_df_it1$padj < 0.05 & abs(deg_df_it1$log2FC) > 0.25)], 
       "II vs NN in Colon" = deg_df_cr1$Term[which(deg_df_cr1$padj < 0.05 & abs(deg_df_cr1$log2FC) > 0.25)]),
  filename = NULL, 
  fill=c("#0073C2FF", "#E41A1C"), 
  cex  = 1.5, 
  cat.cex = 1.8,
  height = 6, 
  width = 6, 
  units = "in",
  resolution = 300
)

pdf(glue("{plt_path}/venn_{ann_col}_{ctx}_TI_Colon.pdf"),
    height = 6, width = 6)
grid::grid.draw(vd)
dev.off()


##### Fig. 2C and D: Inflammation_score #####
library(GSVA)
library(GSEABase)

## 
use_de <- "DESeq2"
ann_col <- "ann_level2_refine"
ctx <- "CD4T"

deg_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{deg_path}TI/all/{use_de}/")
deg_path_cr <- glue("{deg_path}Colon/all/{use_de}/")
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]
deg_df_cr1 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[[ctx]]

##
inter_deg_up <- intersect(deg_df_it1$Term[which(deg_df_it1$padj < 0.05 &
                                                  deg_df_it1$log2FC > 0.25)],
                          deg_df_cr1$Term[which(deg_df_cr1$padj < 0.05 &
                                                  deg_df_cr1$log2FC > 0.25)])
inter_deg_down <- intersect(deg_df_it1$Term[which(deg_df_it1$padj < 0.05 & 
                                                    deg_df_it1$log2FC < -0.25)],
                            deg_df_cr1$Term[which(deg_df_cr1$padj < 0.05 & 
                                                    deg_df_cr1$log2FC < -0.25)])
inter_deg <- c(inter_deg_up, inter_deg_down)

##
sample_col <- "Sample_ID_exp"
bulk_meta <- sc_meta_CD4[!duplicated(sc_meta_CD4[[sample_col]]),
                         !duplicated(colnames(sc_meta_CD4))]
bulk_meta$group_comb  <- paste0(bulk_meta$Condition, ":", bulk_meta$Section)
sample_list <- split(bulk_meta$Sample_exp, f = bulk_meta$group_comb)
rownames(bulk_meta) <- gsub("_", "-", bulk_meta$Sample_ID_exp)
#
clinic_df <- readRDS("03_output/02_clean/clinic_df.rds")
exp_bulk <- AggregateExpression(sc_obj, 
                                assays = "RNA", 
                                return.seurat = F, 
                                group.by = "Sample_ID_exp")[["RNA"]]
exp_bulk_lcpm <- edgeR::cpm(exp_bulk, 
                            normalized.lib.sizes = T,
                            log = T)
## GSVA 
inter_set <- GeneSetCollection(GeneSet(inter_deg, setName = "inter_deg"))
ssgsva_param <- ssgseaParam(exprData = exp_bulk_lcpm,
                            geneSets = inter_set,
                            normalize = F)
gsva_df <- gsva(ssgsva_param) %>%
  t %>% as.data.frame()
#
bulk_meta$score_CD4 <- gsva_df[rownames(bulk_meta), "inter_deg"] %>% scales::rescale(., c(0, 10))
# add NI_score
ssgsea_df  <- readRDS("ssgsva_NI_score.rds")
bulk_meta$if_score_scale <- ssgsea_df[rownames(bulk_meta), "NIScore"] %>% scales::rescale(., c(0, 10))

## plot
bulk_meta$Section <- factor(bulk_meta$Section, 
                            levels = c("Colon", "Ileum"))
bulk_meta$group <- factor(bulk_meta$Condition, 
                          levels = c("NN", "NU", "II"))
##
ppair_df <- lapply(1:4, function(x){
  
  sample_listx <- pair_list[[x]]
  sample_pair_use <- intersect(sample_list[[sample_listx[1]]], sample_list[[sample_listx[2]]])
  bulk_metax <- subset(bulk_meta, Sample_exp %in% sample_pair_use)
  bulk_metax$group <- droplevels(bulk_metax$group)
  group_level <- levels(bulk_metax$group)
  
  pair_dfx <- split(bulk_metax, f = bulk_metax$Sample_exp) %>%
    lapply(., function(x){
      propx <- x$score_CD4[match(group_level, x$group)]
    }) %>% Reduce("rbind", .)
  ppairx <- wilcox.test(pair_dfx[,1], pair_dfx[,2], paired = T)$p.value
  
  return(ppairx)
  
}) %>% unlist()

padjpair_df <- p.adjust(ppair_df, method = "BH")
names(ppair_df) <- names(padjpair_df) <- names(pair_list)[1:4]

##
padj <- format(padjpair_df, scientific = T, digits = 3) %>%
  gsub("e", "E", .) %>% as.vector()
diff_df <- data.frame(y.position = c(10, 9, 10, 9),
                      x = c(1.5, 2.5, 1.5, 2.5),
                      xmin = c(2, 2, 2, 2),
                      xmax = c(1, 3, 1, 3),
                      group2 = c("NN", "NN", "NN", "NN"),
                      group1 = c("II", "NU", "II", "NU"),
                      Section = c("Ileum", "Ileum", "Colon", "Colon"),
                      padj_sig = padj)
# dd <- bulk_meta[,c("group", "score_CD4", "Sample_exp", "Section")]
if_box <- ggplot(bulk_meta, aes(x = group, y = score_CD4, color = group)) +
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
  facet_wrap(~ Section, scales = 'fixed', nrow = 1) +
  xlab("") + ylab("Inflammation score (CD4+T cells)") + 
  theme_bw() + 
  theme(legend.position = "none",
        title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold"))
ggsave(glue("{plt_path}/if_box_CD4.png"),
       if_box,
       height = 6, width = 6, units = "in", dpi = 300)
##
if_cor <- ggplot(bulk_meta, aes(x = if_score_scale, y = score_CD4)) + 
  geom_point(aes(color = group)) + 
  geom_smooth(method = "lm", se = T) + 
  facet_wrap(~Section, ncol = 1, scales = "free") +
  scale_color_manual(values = c("II" = "#E41A1C", "NN" = "#0073C2FF", "NU" = "#EFC000FF")) +
  xlab("Inflammation score (Thomas)") + 
  ylab("Inflammation score (CD4+T cells)") + 
  theme_bw() + 
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 13, face = "bold"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15, face = "bold"))
ggsave(glue("{plt_path}/if_cor.png"),
       if_cor,
       height = 8, width = 5, units = "in", dpi = 300)
# add test
cor_test <- split(bulk_meta, f = bulk_meta$Section) %>% 
  lapply(., function(x){
    corx <- cor.test(x$score_CD4, x$if_score_scale, method = "spearman")
    return(c("rho" = corx$estimate, "p" = corx$p.value))
    }) %>%
  Reduce("rbind", .) %>% as.data.frame()
write.table(cor_test, file = glue("{plt_path}/if_cor.txt"))

##### Fig. 2E: motif entichment #####
motif_enrich_it1 <- readRDS("03_output/04_Diff/DAR/TI/all/DESeq2/motif_enrich_df_CD4T_II_vs_NN.rds")
motif_enrich_it1 <- subset(motif_enrich_it1, fold.enrichment > 1)
motif_enrich_it1$log2FC <- ifelse(motif_enrich_it1$dir == "Up", 
                                 log2(motif_enrich_it1$fold.enrichment), 
                                 -log2(motif_enrich_it1$fold.enrichment))
motif_enrich_it1 <- motif_enrich_it1[order(motif_enrich_it1$p.adjust),]
motif_enrich_it1 <- motif_enrich_it1[!duplicated(motif_enrich_it1$motif),]
#
motif_enrich_it2 <- readRDS("03_output/04_Diff/DAR/TI/all/DESeq2/motif_enrich_df_CD4T_NU_vs_NN.rds")
motif_enrich_it2 <- subset(motif_enrich_it2, fold.enrichment > 1)
motif_enrich_it2$log2FC <- ifelse(motif_enrich_it2$dir == "Up", 
                                    log2(motif_enrich_it2$fold.enrichment), 
                                    -log2(motif_enrich_it2$fold.enrichment))
motif_enrich_it2 <- motif_enrich_it2[order(motif_enrich_it2$p.adjust),]
motif_enrich_it2 <- motif_enrich_it2[!duplicated(motif_enrich_it2$motif),]
## plot
source("code/FUNCTION/diff/volca_plot.R")
plt_inter <- bi.voca.plot(deg_df1 = motif_enrich_it1,
                          deg_df2 = motif_enrich_it2,
                          name_df1 = "Ileal II vs NN",
                          name_df2 = "Ileal NU vs NN",
                          p_col = "p.adjust",
                          FC_col = "log2FC",
                          term_col = "motif.name",
                          target_gene = c("BATF", "BACH1"),
                          n_top = 30,
                          p_thre = 0.05,
                          logFC_thre = 1,
                          col_use = c("#F02720", "#E9C39B", "#B5C8E2", "grey"))
ggsave(glue("{plt_path}voca_inter_motif_enrich_CD4T.png"),
       plt_inter,
       height = 6, width = 8, units = "in", dpi = 300)

##### Fig. 2F: Heatmap of sel SCENIC+ regulon #####
## load SCENIC+ object
sp_outpath <- glue("03_output/07_SCENIC/CD4T/scplus_pipeline/Snakemake/")
scAUC_obj <- readRDS(glue("{sp_outpath}scAUC_obj.rds"))
## load differential test results
diff_auc_ti <- rbind(readRDS("03_output/04_Diff/AUC/TI/CD4T/scAUC_pseudo_pair_ann_level2_final_II_vs_NN.rds"),
                            readRDS("03_output/04_Diff/AUC/TI/CD4T/scAUC_pseudo_pair_ann_level2_final_NU_vs_NN.rds"))
diff_auc_ti <- diff_auc_ti[!is.na(diff_auc_ti$padj),]
# selected Regulons
diff_auc_ti_sel <- diff_auc_ti[grepl("ROR|BACH1|BATF|MAF", diff_auc_ti$Term) &
                                 diff_auc_ti$padj < 0.05,]
sel_er <- unique(diff_auc_ti_sel$Term)
scAUC_obj_sel <- subset(scAUC_obj, 
                     Section == "TI", 
                     features = sel_er)
scAUC_obj_sel <- ScaleData(scAUC_obj_sel)
# pseudobulk AUC matrix
AUC_mat_agg <- AverageExpression(scAUC_obj_sel, 
                                  group.by = "Condition", 
                                  return.seurat = F, 
                                  assays = "AUC", 
                                  layer = "scale.data")[["AUC"]] %>% 
  as.data.frame()
AUC_mat_agg <- AUC_mat_agg[,c("NN", "NU", "II"), drop = F]
rownames(AUC_mat_agg) <- gsub("\\-\\(", "_\\(", rownames(AUC_mat_agg))
## create symbol mat
symbol_matx <- matrix("", 
                      ncol = ncol(AUC_mat_agg), 
                      nrow = nrow(AUC_mat_agg))
dimnames(symbol_matx) <- dimnames(AUC_mat_agg)
for (nx in 1:nrow(diff_auc_ti_sel)) {
  
  termx <- diff_auc_ti_sel[nx, "Term"] %>%
    gsub("\\-\\(", "_\\(", .)
  posx <- diff_auc_ti_sel[nx, "avg_diff"] > 0
  sig_g <- diff_auc_ti_sel[nx, "contrast"] %>% 
    stringr::str_split(., "\\-", simplify = T)
  if (posx) {
    sig_g <- sig_g[,1:2, drop = T]
  } else {
    sig_g <- sig_g[,2:1, drop = T]
  }
  #
  if (symbol_matx[termx, sig_g[1]] == "-") {
    symbol_matx[termx, sig_g[1]] <- "+/-"
  } else {
    symbol_matx[termx, sig_g[1]] <- "+"
  }
  #
  if (symbol_matx[termx, sig_g[2]] == "+") {
    symbol_matx[termx, sig_g[2]] <- "+/-"
  } else {
    symbol_matx[termx, sig_g[2]] <- "-"
  }
}
symbol_matx[,"NN"] <- ""

## plot heatmap
clust_row <- nrow(AUC_mat_agg) > 1
png(glue("{plt_path}/heatAUC_sel_CD4T.png"),
    height = 3, width = 6, 
    units = "in", res = 300)
pheatmap::pheatmap(t(AUC_mat_agg),
                   scale = "none", 
                   legend_labels = c("min","max"),
                   fontsize = 12,
                   display_numbers = t(symbol_matx),
                   fontsize_number = 15, 
                   treeheight_col = 0,
                   cluster_rows = F,
                   treeheight_row = 0,
                   cluster_cols = T, 
                   angle_col = 90) %>% print()
dev.off()

##### Fig. 2G: UMAP #####
umap_plt <- DimPlot(sc_obj, 
                    cols = ct_cols,
                    label = F, 
                    raster = F, 
                    pt.size = 0.1, 
                    alpha = 0.8)+
  xlab("wnnUMAP1") + ylab("wnnUMAP2")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#
ggsave(file = glue("{plt_path}UmapPlot_CD4T.png"),
       umap_plt,
       width = 9.5, height = 6,units = "in", dpi = 300)

##### Fig. 2H: Makers #####
# markers
adt_cd4 <- c("CD45RA","CD45RO", 
             "CD62L", "CD49a", "CD69", "CD103", "CD278", "CD279", 
             "CD25", "CD39")
rna_cd4 <- c("CCR7", "SELL", 
             "CXCR5", "CXCL13", "GZMK", "NKG7", "GZMB",
             "CCL5", "CCR6", "IL23R", "IL17A",
             "IKZF2", "FOXP3")
atac_cd4 <- c("CCL5", "IFNG", "IL17A", "RORC", "GZMB", "IKZF2", "FOXP3")
sel_motif <- c("MA0690.1", "MA1151.1", "MA0071.1", "MA0072.1")
names_motif <- c("TBX21", "RORC", "RORA.1", "RORA.2")

Idents(sc_obj) <- sc_obj$ann_level4_final
# ADT
dot_adt_cd4 <- vln.plt(seurat_obj = sc_obj,
                       clus_col = "ann_level4_final",
                       assay_use = "ADT",
                       min_exp = 5,
                       feature_list = adt_cd4,
                       color_use = ct_cols)
# RNA
dot_rna_cd4 <- dot.minmax(seurat_obj = sc_obj,
                          assay = "RNA",
                          features = rna_cd4,
                          col_min = 0,
                          col_max = 1,
                          group = "ann_level4_final") + 
  scale_fill_distiller(type = 'div', palette = 'RdYlBu')+
  scale_y_discrete(limits = levels(sc_obj$ann_level4_final) %>% rev)  + 
  theme(legend.position = "none",
        axis.text.y =  element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())
# ATAC
cv_plt_bk <- cv_plt <- peak.set.plt(sc_obj = sc_obj,
                                    assay_use = "peaks",
                                    mk_list = atac_cd4,
                                    extend_kb_up = 3000,
                                    extend_kb_down = 0,
                                    color_use = ct_cols)
cv_plt[[1]] <- cv_plt[[1]]  + 
  theme(strip.text.y.left = element_blank())
# motif
dot_mt <- DotPlot(sc_obj, 
                  assay = "chromvar", 
                  col.min = -2,
                  col.max = 2,
                  features = sel_motif,
                  group.by = "ann_level4_final",
                  cols = c("lightgrey", "brown")) + 
  scale_x_discrete(limits = sel_motif, labels = names_motif) + 
  scale_y_discrete(limits = levels(sc_obj$ann_level4_final) %>% rev) + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black",
                                   angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# plot
dot_plt <- (dot_adt_cd4 | dot_rna_cd4 | cv_plt | dot_mt) + 
  plot_layout(widths = c(length(adt_cd4) + 2, 
                         length(rna_cd4),
                         length(atac_cd4) * 2.5,
                         length(sel_motif)))
ggsave(file = glue("{plt_path}dot_plt_CD4.png"),
       dot_plt,
       height = 6, width = 18,units = "in", dpi = 300)

##### Fig. 2I and Extended Fig. 2: pie plot #####
pie_plt <- pt.plt.fun(sc_meta = sc_meta_CD4,
                      ann_col = "ann_level4_final",
                      wd_group_col = "Condition",
                      ht_group_col = "Section",
                      color_use = ct_cols) + 
  theme(legend.position = "none")
ggsave(file = glue("{plt_path}pie_plt_CD4T.png"),
       pie_plt,
       width = 6, height = 4,units = "in", dpi = 300)

# restrict to non-treat patients
clinic_df <- readRDS("03_output/02_clean/clinic_df.rds")
sample_notreat <- clinic_df$Sample_ID_exp[rowSums(clinic_df[,grep("^use_", colnames(clinic_df))]) == 0]
sc_meta_notreat <- subset(sc_meta_CD4, Sample_ID_exp %in% sample_notreat)

pie_plt_notreat <- pt.plt.fun(sc_meta = sc_meta_notreat,
                              ann_col = "ann_level4_final",
                              wd_group_col = "Condition",
                              ht_group_col = "Section",
                              color_use = ct_cols) + 
  guides(fill = guide_legend(ncol = 2)) + 
  theme(legend.position = "right")
ggsave(file = glue("{plt_path}pie_plt_CD4T_notreat.png"),
       pie_plt_notreat,
       width = 10, height = 4,units = "in", dpi = 300)

##### Fig. 2J and Extended Fig. 2:  miloR #####
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
  sc_obj_sub$ann_level4_final <- factor(sc_obj_sub$ann_level4_final,
                                         levels = levels(sc_obj_sub$ann_level4_final),
                                         labels = gsub("_", " ",
                                                       levels(sc_obj_sub$ann_level4_final)))
  #
  res_milo <- milor.seurat(seurat_obj = sc_obj_sub,
                           n_dims = 30,
                           k_nn = 50,
                           group_col = "Condition",
                           cov_col = "Sample_exp",
                           block_col = NULL,
                           sample_col = "Sample_ID_exp",
                           ann_col = "ann_level4_final",
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
         height = 7, width = 15, units = "in", dpi = 300)
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
sc_obj_sub$Section <- factor(sc_obj_sub$Section, levels = c("Ileum", "Colon"))
sc_obj_sub$ann_level4_final <- factor(sc_obj_sub$ann_level4_final,
                                       levels = levels(sc_obj_sub$ann_level4_final),
                                       labels = gsub("_", " ",
                                                     levels(sc_obj_sub$ann_level4_final)))
#
res_milo <- milor.seurat(seurat_obj = sc_obj_sub,
                         n_dims = 30,
                         k_nn = 50,
                         group_col = "Section",
                         cov_col = "Sample_exp",
                         block_col = NULL,
                         sample_col = "Sample_ID_exp",
                         ann_col = "ann_level4_final",
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
       height = 7, width = 15, units = "in", dpi = 300)

##### Fig. 2K and Extended Fig. 2: Diff prop ######
sample_col <- "Sample_ID_exp"
block_col <- "Sample_exp"
ann_col <-  "ann_level4_final"
group_col <- "Condition"
section_col <- "Section"

## 0. format proportion table
bulk_meta <- sc_meta_CD4[!duplicated(sc_meta_CD4[[sample_col]]),]
bulk_meta$group_comb  <- paste0(bulk_meta[[group_col]], ":", bulk_meta[[section_col]])
sample_list <- split(bulk_meta$Sample_exp, f = bulk_meta$group_comb)
#
tab_prop <- table(sc_meta_CD4[[ann_col]], sc_meta_CD4[[sample_col]]) %>% 
  prop.table(2) %>% as.data.frame()
colnames(tab_prop) <- c("cluster", "sample", "Proportion")
tab_prop$Proportion <- tab_prop$Proportion * 100
tab_prop$block <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                            block_col]
tab_prop$group <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                            group_col]
tab_prop$Section <- bulk_meta[match(tab_prop$sample, bulk_meta[[sample_col]]),
                              section_col] %>%
  factor(., levels = c("Ileum", "Colon"))
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
sel_ct <- c("CD4_naiveT", "CD4_CD69low_TCM", "CD4_TCM_Tfh_like", 
            "CD4_CD103_TRM_Th17", "CD4_Treg_naive", "CD4_IKZF2low_Treg", "CD4_Treg")
box_pair_prop_list1_TI <- lapply(sel_ct, function(ctx){
  
  ctx <- gsub("_", " ", ctx)
  plt_dfx <- subset(tab_prop, cluster == ctx & Section == "Ileum")
  plt_dfx$Condition <- plt_dfx$group
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df1[ctx,2:1], scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
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
    xlab("") + ylab("") + 
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

box_pair_prop_list1_TI[[1]] <- box_pair_prop_list1_TI[[1]] + 
  ylab("Proportion (%) in Ileum")
ggsave(glue("{plt_path}box_pair_level4_Ileum_sel.png"),
       patchwork::wrap_plots(box_pair_prop_list1_TI, nrow = 1),
       height = 4, width = 14, units = "in", dpi = 300)

## 2.2 plot for Colon
box_pair_list1_Colon <- lapply(all_ct, function(ctx){
  
  ctx <- gsub("_", " ", ctx)
  plt_dfx <- subset(tab_prop, cluster == ctx & Section == "Colon")
  plt_dfx$Condition <- plt_dfx$group
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df1[ctx, 4:3], scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
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
    xlab("") + ylab("") + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 9, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})
box_pair_list1_Colon[[6]] <- box_pair_list1_Colon[[6]] + ylab("Proportion (%) in Colon")
ggsave(glue("{plt_path}box_pair_level4_Colon.png"),
       patchwork::wrap_plots(box_pair_list1_Colon, nrow = 3),
       height = 11, width = 11, units = "in", dpi = 300)

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
box_pair_list2 <- lapply(all_ct, function(ctx){
  
  plt_dfx <- subset(tab_prop, cluster == ctx & group == "NN")
  plt_dfx$Condition <- plt_dfx$Section
  max_propx <-max(plt_dfx$Proportion)
  padj <- format(ppair_prop_adj_df2[ctx, 1], scientific = T, digits = 3) %>%
    gsub("e", "E", .) %>% as.vector()
  diff_dfx <- data.frame(y.position = c(max_propx * 1.1),
                         x = c(1.5),
                         xmin = c(1),
                         xmax = c(2),
                         group2 = c("Ileum"),
                         group1 = c("Colon"),
                         padj_sig = padj)
  
  box_pairx <- ggplot(data = plt_dfx, 
                      aes(x = Section, y = Proportion, color = Section)) + 
    geom_boxplot(width = 0.5) +
    geom_jitter(width = 0.2) + 
    scale_color_manual(values = c("#0073C2FF", "#E41A1C")) + 
    scale_y_continuous(limits = c(0, max_propx * 1.3)) + 
    ggtitle(ctx) +
    xlab("") + ylab("") + 
    stat_pvalue_manual(diff_dfx,
                       label = "padj_sig",
                       tip.length = 0.01,
                       label.size = 4,
                       bracket.size = 1) +
    theme_bw() + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          title = element_text(size = 9, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face = "bold"))
  
})
box_pair_list2[[6]] <- box_pair_list2[[6]] + ylab("Proportion (%) in NN")
ggsave(glue("{plt_path}box_pair_level4_NN.png"),
       patchwork::wrap_plots(box_pair_list2, nrow = 3),
       height = 11, width = 11, units = "in", dpi = 300)

##### Fig. 2L: Distribution of DEGs #####
use_de <- "DESeq2"
ann_col <- "ann_level2_refine"
deg_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{deg_path}TI/all/{use_de}/")
deg_df_it1 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4T"]]
deg_df_it2 <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_NU_vs_NN.rds"))[["CD4T"]]
inter_gene_up <- intersect(deg_df_it1$Term[which(deg_df_it1$padj < 0.05 & deg_df_it1$log2FC > 0.25)],
                         deg_df_it2$Term[which(deg_df_it2$padj < 0.05 & deg_df_it2$log2FC > 0.25)])
inter_gene_dowm <- intersect(deg_df_it1$Term[which(deg_df_it1$padj < 0.05 & deg_df_it1$log2FC < -0.25)],
                         deg_df_it2$Term[which(deg_df_it2$padj < 0.05 & deg_df_it2$log2FC < -0.25)])

DefaultAssay(sc_obj) <- "RNA"
sc_obj <- AddModuleScore(sc_obj, 
                         features = list(inter_gene_up, 
                                         inter_gene_dowm), 
                         name = c("Up in II/NU vs NN",
                                  "Down in II/NU vs NN"))
ft_theme <- theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ft_deg <- FeaturePlot(sc_obj,
                      reduction = "wnn.umap",
                      raster = F,
                      cols  = c("lightgrey","#FF9500", "#F66C2D"), 
                      features = c("Up in II/NU vs NN1", "Down in II/NU vs NN2"),
                      ncol = 2)
ft_deg[[1]] <- ft_deg[[1]] + 
  ggtitle("Up in Ileal II/NU vs NN") + 
  xlab("wnnUMAP1") + ylab("wnnUMAP2")+
  ft_theme
ft_deg[[2]] <- ft_deg[[2]] + 
  ggtitle("Down in Ileal II/NU vs NN") + 
  xlab("wnnUMAP1") + ylab("wnnUMAP2") +
  ft_theme

ggsave(file = glue("{plt_path}ft_deg.png"),
       ft_deg,
       width = 8, height = 4,units = "in", dpi = 300, limitsize = T)

##### Fig. 2M: scDRS #####
scDRS_out_path <- glue("03_output/01_scDRS/04_downstream/cov/")
##
scDRS_comb_level4 <- fread2(glue("{scDRS_out_path}CD.scdrs_group.ann_level4_final"))
scDRS_comb_level4_CD4 <- subset(scDRS_comb_level4, group %in% levels(sc_meta_CD4$ann_level4_final))
scDRS_comb_level4_CD4$section <- "Combined"
scDRS_comb_level4_CD4$celltype <- factor(scDRS_comb_level4_CD4$group, 
                                   levels = levels(sc_meta_CD4$ann_level4_final) %>% rev,
                                   labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final) %>% rev))
##
scDRS_level4 <- fread2(glue("{scDRS_out_path}CD.scdrs_group.sample_ct_level4"))
scDRS_level4$section <- str_split_i(scDRS_level4$group, ":", 1)
scDRS_level4$celltype <- str_split_i(scDRS_level4$group, ":", 2)
scDRS_level4_CD4 <- subset(scDRS_level4, celltype %in% levels(sc_meta_CD4$ann_level4_final))
scDRS_level4_CD4$celltype <- factor(scDRS_level4_CD4$celltype, 
                                    levels = levels(sc_meta_CD4$ann_level4_final) %>% rev,
                                    labels = gsub("_", " ", levels(sc_meta_CD4$ann_level4_final) %>% rev))
## 
scDRS_CD4_plt <- rbind(scDRS_comb_level4_CD4, scDRS_level4_CD4)
scDRS_CD4_plt$section <- factor(scDRS_CD4_plt$section, 
                                 levels = c("Combined", "TI", "RC"),
                                 labels = c("Combined", "Ileum", "Colon"))
scDRS_CD4_plt <- split(scDRS_CD4_plt, f = scDRS_CD4_plt$section) %>%
  lapply(., function(scDRS_CD4_pltx){
    scDRS_CD4_pltx$FDR <- p.adjust(scDRS_CD4_pltx$assoc_mcp, method = "BH")
    scDRS_CD4_pltx$hetero_padj <- p.adjust(scDRS_CD4_pltx$hetero_mcp, method = "BH")
    return(scDRS_CD4_pltx)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
scDRS_CD4_plt$assoc_sig <- scDRS_CD4_plt$FDR < 0.05
scDRS_CD4_plt$hetero_sig <- scDRS_CD4_plt$hetero_padj < 0.05


scdrs_plt <- ggplot() + 
  geom_tile(data = scDRS_CD4_plt, 
            aes(x = section, y = celltype, fill = -log10(FDR))) + 
  scale_fill_viridis_c() + 
  geom_tile(data = subset(scDRS_CD4_plt, assoc_sig), 
            aes(x = section, y = celltype),
            color = "black", fill = NA, size = 1, show.legend = T) +
  geom_text(data = subset(scDRS_CD4_plt, hetero_sig),
            aes(x = section, y = celltype, label = "\u2716"),
            color = "black", size = 6, show.legend = T) + 
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 10))
#
ggsave(file = glue("{plt_path}CD_sample_ct_CD4T.png"),
       scdrs_plt,
       width = 6, height = 8,units = "in", dpi = 300, limitsize = F)
##### Extended Fig. 2A: Volcano plot on DEGs between Colonic II vs NN #####
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_cr <- glue("{out_path}Colon/all/{use_de}/")
ann_col <- "ann_level2_refine"
source("/ix1/wchen/xiangyu/Projects/03_CD_DOGMA/code/volca_plot.R")
## 1. Ileum
deg_df_it_cd4 <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds"))[["CD4T"]]
deg_df_it_cd4 <- subset(deg_df_it_cd4, !is.na(deg_df_it_cd4$padj))
vc_plt_cd4 <- volca_function(deg_df = deg_df_it_cd4,
                             term_col = "Term",
                             fc_col = "log2FC",
                             p_col = "padj",
                             # target_gene = NULL,
                             n_top = 30,
                             log_thresh = 0.25, 
                             top_log_thresh = 0.5,
                             group_label = c("Down regulation", "Others", "Up regulation"),
                             valco_col = c4a("classic_blue_red12", 3),
                             highlight_col = "orange") + 
  xlab(bquote(~Log[2]~"(Fold Change) in Colonic II vs NN"))

ggsave(file = glue("{plt_path}/valcano_CD4_II_vs_NN_Colon.png"),
       vc_plt_cd4,
       width = 10, height = 6,units = "in", dpi = 300, limitsize = T)

##### Extended Fig. 2C: starCAT #####
tcat_mat <- fread2("03_output/03_clustering/CD4T/TCAT/usage.txt") %>%
  tibble::column_to_rownames("V1")
sc_obj[["TCAT"]] <- CreateAssayObject(t(tcat_mat))
DefaultAssay(sc_obj) <- "TCAT"
ft_tcat <- FeaturePlot(sc_obj, 
                       reduction = "wnn.umap", 
                       raster = F, 
                       features = c("CD4-Naive", "CD4-CM", 
                                    "Th1-Like", "Th17-Resting", "Th17-Activated", 
                                    "Tfh-1", "Tfh-2", "Treg"),
                       max.cutoff = "q99", min.cutoff = "q1",
                       ncol = 2)
ggsave(glue("{plt_path}feature_tcat.png"), 
       ft_tcat,
       height = 20, width = 12, units = "in", dpi = 300)
##
DefaultAssay(sc_obj) <- "RNA"
sc_obj[["TCAT"]] <- NULL

##### Supp. Fig. 2 #####
#### 1. DEGs
use_de <- "DESeq2"
out_path <- "03_output/04_Diff/DEG/"
deg_path_it <- glue("{out_path}TI/CD4T/{use_de}/")
deg_path_cr <- glue("{out_path}Colon/CD4T/{use_de}/")
ann_col <- "ann_level4_final"
source("code/FUNCTION/diff/volca_plot.R")
## 1.1 DEGs bewteen II and NN in Ileal CD4T
deg_df_ti <- readRDS(glue("{deg_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
deg_df_ti$cluster <- gsub("_", " ", deg_df_ti$cluster)
deg_df_ti$cluster <- factor(deg_df_ti$cluster, 
                            levels = levels(sc_meta_CD4[[ann_col]]))
plt_deg_ti <- diff.num.plot(diff_df = deg_df_ti,
                            cluster_col = "cluster",
                            p_adj_col = "padj",
                            fc_col = "log2FC",
                            p_adj_thresh = 0.05,
                            logfc_thresh = 0.25,
                            col_use = c("Down" = "#2C69B0", "Up" = "#F02720")) + 
  xlab("Number of DEGs in Ileal II vs NN")
#
ggsave(file = glue("{plt_path}/ndeg_{ann_col}_II_vs_NN_ti.pdf"),
       plt_deg_ti,
       height = 7, width = 8,units = "in")
## 1.2 DEGs bewteen II and NN in Colonic CD4T
deg_df_cr <- readRDS(glue("{deg_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
deg_df_cr$cluster <- gsub("_", " ", deg_df_cr$cluster)
deg_df_cr$cluster <- factor(deg_df_cr$cluster, 
                            levels = levels(sc_meta_CD4[[ann_col]]))
plt_deg_cr <- diff.num.plot(diff_df = deg_df_cr,
                            cluster_col = "cluster",
                            p_adj_col = "padj",
                            fc_col = "log2FC",
                            p_adj_thresh = 0.05,
                            logfc_thresh = 0.25,
                            col_use = c("Down" = "#2C69B0", "Up" = "#F02720")) + 
  xlab("Number of DEGs in Colonic II vs NN")
#
ggsave(file = glue("{plt_path}/ndeg_{ann_col}_II_vs_NN_cr.pdf"),
       plt_deg_cr,
       height = 7, width = 8,units = "in")
#### 2. DARs
## 2.1 DARs bewteen II and NN in Ileal CD4T
out_path <- "03_output/04_Diff/DAR/"
dar_path_it <- glue("{out_path}TI/CD4T/{use_de}/")
dar_path_cr <- glue("{out_path}Colon/CD4T/{use_de}/")
ann_col <- "ann_level4_final"
##
dar_df_ti <- readRDS(glue("{dar_path_it}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
dar_df_ti$cluster <- gsub("_", " ", dar_df_ti$cluster)
dar_df_ti$cluster <- factor(dar_df_ti$cluster, 
                            levels = levels(sc_meta_CD4[[ann_col]]))
plt_dar_ti <- diff.num.plot(diff_df = dar_df_ti,
                            cluster_col = "cluster",
                            p_adj_col = "padj",
                            fc_col = "log2FC",
                            p_adj_thresh = 0.05,
                            logfc_thresh = 0.25,
                            col_use = c("Down" = "#2C69B0", "Up" = "#F02720")) + 
  xlab("Number of DARs in Ileal II vs NN")
#
ggsave(file = glue("{plt_path}/ndar_{ann_col}_II_vs_NN_ti.pdf"),
       plt_dar_ti,
       height = 7, width = 8,units = "in")

## 2.2 DARs bewteen II and NN in Colonic CD4T
dar_df_cr <- readRDS(glue("{dar_path_cr}/pseudo_bulk_{ann_col}_II_vs_NN.rds")) %>%
  Reduce("rbind", .)
dar_df_cr$cluster <- gsub("_", " ", dar_df_cr$cluster)
dar_df_cr$cluster <- factor(dar_df_cr$cluster, 
                            levels = levels(sc_meta_CD4[[ann_col]]))
plt_dar_cr <- diff.num.plot(diff_df = dar_df_cr,
                            cluster_col = "cluster",
                            p_adj_col = "padj",
                            fc_col = "log2FC",
                            p_adj_thresh = 0.05,
                            logfc_thresh = 0.25,
                            col_use = c("Down" = "#2C69B0", "Up" = "#F02720")) + 
  xlab("Number of DARs in Colonic II vs NN")
#
ggsave(file = glue("{plt_path}/ndar_{ann_col}_II_vs_NN_cr.pdf"),
       plt_dar_cr,
       height = 7, width = 8,units = "in")


