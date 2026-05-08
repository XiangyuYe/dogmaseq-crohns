##
library(Signac)
library(GenomicRanges)
library(CMplot)

##
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
## DEG
gene_coord <- fread2("/ix1/wchen/xiangyu/Ref_data/genome/cellranger/arc_GRCh38/genes/gene_coordinates.tsv")
gene_coord$chr <- gene_coord$seqnames %>% 
  gsub("chr", "", .) %>% as.integer()
deg_cd103trm17 <- readRDS("03_output/04_Diff/DEG/small_intestine/CD4T/DESeq2/pseudo_bulk_ann_level4_refine_I_I_vs_N_N.rds")[["CD4_CD103_TRM_Th17"]]
deg_cd103trm17$chr <- gene_coord$seqnames[match(deg_cd103trm17$Term, gene_coord$gene_name)] %>%
  gsub("chr", "", .) %>% as.integer()
# map to nearest SNP pos
gr_gene <- GRanges(seqnames = gene_coord$seqnames,
                   ranges = IRanges(start = gene_coord$start,
                                    end = gene_coord$end))
idx_gene <- IRanges::nearest(gr_gene, gr_gwas)
gene_coord$pos <- summ_dff$pos[idx_gene]
gene_coord$rsid_nearest <- summ_dff$rsid[idx_gene]
deg_cd103trm17$rsid_nearest <- gene_coord$rsid_nearest[match(deg_cd103trm17$Term, gene_coord$gene_name)]
#
# deg_cd103trm17 <- subset(deg_cd103trm17, !is.na(padj) & !is.na(chr) & !is.na(pos))
# deg_cd103trm17 <- deg_cd103trm17[order(deg_cd103trm17$padj),]
# deg_cd103trm17 <- deg_cd103trm17[!duplicated(deg_cd103trm17$rsid_nearest),]

## MAGMA
magma_df <- fread2("03_output/01_scDRS/01_magma/CD.genes.out")
magma_zdf <- fread2("03_output/01_scDRS/01_magma/zscore_CD.tsv")
magma_df$Term <- magma_zdf$GENE
deg_cd103trm17$p_magma <- magma_df$P[match(deg_cd103trm17$Term, magma_df$Term)]

## DAR
dar_cd103trm17 <- readRDS("03_output/04_Diff/DAR/small_intestine/CD4T/DESeq2/pseudo_bulk_ann_level4_refine_I_I_vs_N_N.rds")[["CD4_CD103_TRM_Th17"]]
# map to nearest SNP pos
gr_peak <- StringToGRanges(dar_cd103trm17$Term)
idx_peak <- IRanges::nearest(gr_peak, gr_gwas)
dar_cd103trm17$pos <- summ_dff$pos[idx_peak]
dar_cd103trm17$rsid_nearest <- summ_dff$rsid[idx_peak]
#
dar_cd103trm17 <- subset(dar_cd103trm17, !is.na(padj) & !is.na(pos))
dar_cd103trm17 <- dar_cd103trm17[order(dar_cd103trm17$padj),]
dar_cd103trm17 <- dar_cd103trm17[!duplicated(dar_cd103trm17$rsid_nearest),]

links_gr <- readRDS("03_output/03_clustering_test/CD4T/test153000/test302015/LinkPeaks_links_cd103trm17.rds")
df_links <- as.data.frame(links_gr)
df_links$padj <- p.adjust(df_links$pvalue, method = "BH")
dff_links <- subset(df_links, df_links$padj < 0.05 & df_links$score > 0.1)
dff_links_sig <- subset(dff_links, peak %in% 
                          dar_cd103trm17$Term[which(dar_cd103trm17$padj < 0.05 & abs(dar_cd103trm17$log2FC) > 0.25)])
##
summ_dff$Gene <- deg_cd103trm17$Term[match(summ_dff$rsid, deg_cd103trm17$rsid_nearest)]
summ_dff$log2FC_gene <- deg_cd103trm17$log2FC[match(summ_dff$rsid, deg_cd103trm17$rsid_nearest)]
summ_dff$padj_gene <- deg_cd103trm17$padj[match(summ_dff$rsid, deg_cd103trm17$rsid_nearest)]
summ_dff$p_magma <- deg_cd103trm17$p_magma[match(summ_dff$rsid, deg_cd103trm17$rsid_nearest)]

summ_dff$Peak <- dar_cd103trm17$Term[match(summ_dff$rsid, dar_cd103trm17$rsid_nearest)]
summ_dff$log2FC_peak <- dar_cd103trm17$log2FC[match(summ_dff$rsid, dar_cd103trm17$rsid_nearest)]
summ_dff$padj_peak <- dar_cd103trm17$padj[match(summ_dff$rsid, dar_cd103trm17$rsid_nearest)]

# summ_dff$p_magma[which(summ_dff$p_magma >= 1E-5)] <- NA
# summ_dff$padj_peak[which(summ_dff$padj_peak >= 0.05 | 
# !summ_dff$Peak %in% dff_links$peak)] <- NA
# summ_dff$padj_gene[which(summ_dff$padj_gene >= 0.05)] <- NA
##
# SNP Chromosome Position    trait1     trait2     trait3
plt_df <- data.frame(SNP = summ_dff$rsid,
                     Chromosome = summ_dff$chr,
                     Position = summ_dff$pos,
                     DARs = summ_dff$padj_peak,
                     DEGs = summ_dff$padj_gene)
plt_df <- subset(plt_df, 
                 !is.na(DARs) | !is.na(DEGs))

col_mat <- matrix(c("royalblue4", "darksalmon",
                    "grey", "lightblue"), 
                  nrow = 2, byrow = T)
##
top_gene1 <- deg_cd103trm17$Term[which(deg_cd103trm17$padj < 0.05 & 
                                 abs(deg_cd103trm17$log2FC) > 0.25 &
                                 deg_cd103trm17$p_magma < 1E-5)]
top_gene2 <- deg_cd103trm17$Term[which(deg_cd103trm17$padj < 0.05 & 
                                        abs(deg_cd103trm17$log2FC) > 0.25 &
                                        deg_cd103trm17$Term %in% dff_links_sig$gene)]

top_snp1 <- deg_cd103trm17$rsid_nearest[match(top_gene1, deg_cd103trm17$Term)]
gene_high_col1 <- rep("brown", length(top_gene1))

top_snp2 <- deg_cd103trm17$rsid_nearest[match(top_gene2, deg_cd103trm17$Term)]
top_peak2 <- summ_dff$Peak[which(summ_dff$padj_peak < 0.05 & 
                                   abs(summ_dff$log2FC_peak) > 0.25 &
                                   summ_dff$Peak %in% dff_links_sig$peak)]
top_snp22 <- summ_dff$rsid[which(summ_dff$padj_peak < 0.05 & 
                                   abs(summ_dff$log2FC_peak) > 0.25 &
                                   summ_dff$Peak %in% dff_links_sig$peak)]
gene_high_col2 <- rep("orange", length(top_gene2))

CMplot(plt_df,
       type = "p",
       plot.type = "m",
       r = 1,
       col = col_mat,
       cir.chr.h = 1.5,
       amplify = F,
       signal.line = 0.2,
       threshold = 0.05,
       threshold.col = "red",
       threshold.lty = 2,   
       highlight = list("DARs" = top_snp22,
                        "DEGs" = c(top_snp1, top_snp2)),
       highlight.col = list("DARs" = "orange",
                            "DEGs" = c(gene_high_col1, gene_high_col2)),
       highlight.text = list("DARs" = top_peak2,
                             "DEGs" = c(top_gene1, top_gene2)),
       highlight.cex = 1,
       highlight.text.cex = 2,
       highlight.text.col = "black",
       outward = T,
       file = "pdf",
       dpi = 300,
       file.output = T,
       verbose = T,
       width = 15,
       height = 7,
       multracks = T)
