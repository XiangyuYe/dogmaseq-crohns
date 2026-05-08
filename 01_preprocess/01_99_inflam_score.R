## module load r/4.5.0
library(clusterProfiler)
library(org.Hs.eg.db)

## set path
Project_path="/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(Project_path)

gene_list <- read.table("Inflammation_score_genes.txt")[,1, drop = T]
eg <- try(bitr(gene_list, 
               fromType="ENSEMBL", 
               toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
               OrgDb="org.Hs.eg.db"),silent = T)
saveRDS(eg, file = "Inflammation_score_genes.rds")

## inflammatory scores
library(Seurat)
library(dplyr)
library(GSVA)
library(GSEABase)
library(BiocParallel)

## 0. load gene list
gene_list_ni <- readRDS("Inflammation_score_genes.rds")
NI <- GeneSetCollection(GeneSet(gene_list_ni$SYMBOL, setName="NIScore"))

## 1. generate score via GSVA
scRNA_obj <- readRDS("03_output/02_clean/scRNA_obj_immune.rds")
exp_bulk <- AggregateExpression(scRNA_obj, 
                                assays = "RNA", 
                                return.seurat = F, 
                                group.by = "Sample_ID_exp")[["RNA"]]
exp_bulk_lcpm <- edgeR::cpm(exp_bulk, 
                            normalized.lib.sizes = T,
                            log = T)
ssgsva_param <- ssgseaParam(exprData = exp_bulk_lcpm,
                            geneSets = NI,
                            normalize = F)
gsva_df <- gsva(ssgsva_param,
                BPPARAM = MulticoreParam(workers = 20)) %>%
  t %>% as.data.frame()
saveRDS(gsva_df, file = "ssgsva_NI_score.rds")
