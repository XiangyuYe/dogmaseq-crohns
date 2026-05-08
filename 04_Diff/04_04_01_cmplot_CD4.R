##### Fig. 2A: Circle Manha plot #####
library(GenomicRanges)
library(CMplot)

project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

## diff path
diff_path <- "03_output/04_Diff/"
deg_path_it <- glue("{diff_path}DEG/TI/all/DESeq2/")
dar_path_it <- glue("{diff_path}DAR/TI/all/DESeq2/")
p2g_path_it <- glue("{diff_path}peak2gene/TI/")

## SNP position GRCh 38
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
## Gene coord
gene_coord <- fread2("/ix1/wchen/xiangyu/Ref_data/genome/cellranger/arc_GRCh38/genes/gene_coordinates.tsv")
gene_coord$chr <- gene_coord$seqnames %>% 
  gsub("chr", "", .) %>% as.integer()

## DEG
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

## MAGMA
magma_df <- fread2("03_output/01_scDRS/01_magma/CD.genes.out")
magma_zdf <- fread2("03_output/01_scDRS/01_magma/zscore_CD.tsv")
magma_df$Term <- magma_zdf$GENE
deg_CD4$p_magma <- magma_df$P[match(deg_CD4$Term, magma_df$Term)]

## DAR
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

## combine DEG, DAR, and MAGMA
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
                 !is.na(trait1) | !is.na(trait2) | !is.na(trait3))

## parameters for plot
# dot color for three modalities
col_mat <- matrix(c("#0073C2FF", "#EFC000FF", # MAGMA
                    "royalblue4", "darksalmon", # DAP
                    "grey", "lightblue"), # DEG
                  nrow = 3, byrow = T)

# labels genes: DEGs & MAGMA Sig
top_gene_magma <- deg_CD4$Term[which(deg_CD4$padj < 0.05 & 
                                 abs(deg_CD4$log2FC) > 0.25 &
                                 deg_CD4$p_magma < 1E-5)]
top_gene_magma_pos <- deg_CD4$rsid_nearest[match(top_gene, deg_CD4$Term)]
gene_high_col_magma <- ifelse(deg_CD4$log2FC[match(top_gene, deg_CD4$Term)] > 0,
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
top_gene_link <- setdiff(unique(dff_links_sig$gene), top_gene)
top_gene_link_pos <- deg_CD4$rsid_nearest[match(top_gene22, deg_CD4$Term)]
gene_high_col_link <- rep("orange", length(top_gene_link))

# highlight peaks: DEGs & link to DAPs
top_peak_link <- unique(dff_links_sig$peak)

# label & bold genes: DEGs & link to DAPs & MAGMA Sig
top_gene_link_magma <- intersect(unique(dff_links_sig$gene), top_gene)

## Circle plot
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
                        "DAP" = top_peak_link,
                        "DEG" = c(top_gene_magma_pos, top_gene_link_pos)),
       highlight.col = list("MAGMA" = "brown",
                            "DAP" = "orange",
                            "DEG" = c(gene_high_col_magma, gene_high_col_link)),
       highlight.cex = 1,
       highlight.text = list("MAGMA" = NA,
                             "DAP" = NA,
                             "DEG" = NA),
       highlight.text.cex = 2,
       highlight.text.col = "black",
       outward = T,
       file = "png",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 15,
       height = 15)
system(glue("mv -f Cir_Manhtn.trait1_trait2_trait3.png {plt_path}"))

## manha plot for gene modality only for manual labeling ref
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
       highlight = list("trait3" = c(top_snp)),
       highlight.col = list("trait3" = c(gene_high_col)),
       highlight.cex = 1,
       highlight.text = list("trait3" = c(top_gene)),
       highlight.text.cex = 2,
       highlight.text.col = "black",
       outward = T,
       file = "png",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 15,
       height = 10)
system(glue("mv -f Rect_Manhtn.trait3.png {plt_path}"))
print(top_gene_link_magma)
