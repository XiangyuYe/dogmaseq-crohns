library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ReactomePA)
library(cols4all)
library(bigreadr)

# bar_col <- brewer.pal(3, "Set2")
bar_col <- c4a("10", 3)
names(bar_col) <- c("GO", "KEGG", "Reactome")

Enrichment.pipline <- function(gene_list = NULL,
                               list_name = NULL,
                               num_show = 10,
                               plot = T,
                               GO = T,
                               KEGG = T,
                               Reactome = T,
                               outpath = NULL
                               
){
  
  eg <- try(bitr(gene_list, 
                 fromType="SYMBOL", 
                 toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
                 OrgDb="org.Hs.eg.db"),silent = T)
  if (!is.data.frame(eg)) {
    print("Error: no gene maped!")
  } else {
    if (GO == T) {
      ## 1. enrich
      go_term <- enrichGO(eg$ENTREZID, 
                          OrgDb = org.Hs.eg.db, 
                          ont='ALL',
                          pAdjustMethod = 'BH',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.2,
                          minGSSize = 10,
                          maxGSSize = 500,
                          keyType = 'ENTREZID')
      ## 2. output
      if (is.null(go_term)) {
        print("No significant GO terms were found!")
        write.table(NA, file = paste0(outpath, "/", list_name, "_go.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
      } else if (nrow(go_term) == 0) {
        print("No significant GO terms were found!")
        write.table(NA, file = paste0(outpath, "/", list_name, "_go.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
      } else {
        write.table(go_term, file = paste0(outpath, "/", list_name, "_go.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
        ## 3. plot
        if (plot == T) {
          
          dot_go <- dotplot(go_term, 
                            showCategory = num_show, 
                            orderBy = "x", 
                            font.size=12)
          tiff(paste0(outpath, "/", list_name, "_dot_go.tiff"),
               height = 7,width = 7,
               units="in",res = 600,compression = "lzw")
          print(dot_go)
          dev.off()
          print("GO plot saved!")
        }
      }
    }
    
    if (KEGG == T) {
      ## 1. enrich
      kegg_term <- enrichKEGG(eg$ENTREZID, 
                              # organism = 'hsa',   
                              keyType = 'kegg', 
                              pvalueCutoff = 0.05,
                              pAdjustMethod = 'BH',
                              minGSSize = 10,
                              maxGSSize = 500,
                              qvalueCutoff = 0.2,
                              use_internal_data = FALSE)
      ## 2. output
      if (is.null(kegg_term)) {
        print("No significant KEGG terms were found!")
        write.table(NA, file = paste0(outpath, "/", list_name, "_kegg.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
      } else if(nrow(kegg_term) == 0) {
        print("No significant KEGG terms were found!")
        write.table(NA, file = paste0(outpath, "/", list_name, "_kegg.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
      } else {
        write.table(kegg_term, file = paste0(outpath, "/", list_name, "_kegg.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
        ## 3. plot
        if (plot == T) {
          
          dot_kegg <- dotplot(kegg_term, 
                              showCategory=num_show,
                              orderBy = "x",
                              font.size=12) 
          
          tiff(paste0(outpath, "/", list_name, "_dot_kegg.tiff"),
               height = 7,width = 7,
               units="in",res = 600,compression = "lzw")
          print(dot_kegg)
          dev.off()
          print("KEGG plot saved!")
        }
      }
    }
    
    if (Reactome == T) {
      ## 1. enrich
      reactome_term <- enrichPathway(eg$ENTREZID,
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2,
                                     minGSSize = 10,
                                     maxGSSize = 500,
                                     readable = FALSE
      )
      ## 2. output
      if (is.null(reactome_term)) {
        print("No significant Reactome terms were found!")
        write.table(NA, file = paste0(outpath, "/", list_name, "_reactome.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
      } else if (nrow(reactome_term) == 0) {
        print("No significant Reactome terms were found!")
        write.table(NA, file = paste0(outpath, "/", list_name, "_reactome.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
      } else {
        write.table(reactome_term, file = paste0(outpath, "/", list_name, "_reactome.csv"),
                  quote = F,row.names = F, col.names = T, sep = "\t")
        ## 3. plot
        if (plot == T) {
          
          dot_reactome <- dotplot(reactome_term, 
                                  showCategory = num_show, 
                                  orderBy = "x", 
                                  font.size=12)
          tiff(paste0(outpath, "/", list_name, "_dot_reactome.tiff"),
               height = 7,width = 7,
               units="in",res = 600,compression = "lzw")
          print(dot_reactome)
          dev.off()
          print("Reactome plot saved!")
        }
      }
    }
  }
}



#########################################
enrich_barplt <- function(prefix = NULL,
                          df_path = "./Enrichment/",
                          nn = 8,
                          min_gene = 10,
                          max_nchar = 50){
  go_file <- paste0(df_path, prefix, "_go.csv")
  kegg_file <- paste0(df_path, prefix, "_kegg.csv")
  reactome_file <- paste0(df_path, prefix, "_reactome.csv")
  if (file.exists(go_file)) {
    go_df <- fread2(go_file)
    if (ncol(go_df) > 1) {
      go_df <- go_df[go_df$p.adjust < 0.05,
                     c("ID", "Description", "p.adjust", "GeneRatio", "Count", "geneID")]
      go_df$Group <- "GO"
    } else {
      go_df <- NULL
    }
  } else {
    go_df <- NULL
    message("GO file not found!")
  }
  #
  if (file.exists(kegg_file)) {
    kegg_df <- fread2(kegg_file)
    if (ncol(kegg_df) > 1) {
      kegg_df <- kegg_df[kegg_df$p.adjust < 0.05,
                         c("ID", "Description", "p.adjust", "GeneRatio", "Count", "geneID")]
      kegg_df$Group <- "KEGG"
    } else {
      kegg_df <- NULL
    }
  } else {
    kegg_df <- NULL
    message("KEGG file not found!")
  }
  #
  if (file.exists(reactome_file)) {
    reactome_df <- fread2(reactome_file)
    if (ncol(reactome_df) > 1) {
      reactome_df <- reactome_df[reactome_df$p.adjust < 0.05,
                                 c("ID", "Description", "p.adjust", "GeneRatio", "Count", "geneID")]
      reactome_df$Group <- "Reactome"
    } else {
      reactome_df <- NULL
    }
  } else {
    reactome_df <- NULL
    message("Reactome file not found!")
  }
  #
 enrich_df <- Reduce("rbind", list(go_df, kegg_df, reactome_df))
  if (is.null(enrich_df)) {
    message("No significant terms!")
    return(NA)
  } else {
    enrich_df$Group <- factor(enrich_df$Group,
                               levels = c("GO", "KEGG", "Reactome"))
    
    barplot_df <- subset(enrich_df, nchar(Description) <= max_nchar & Count >= min_gene) %>%
      group_by(., Group) %>%
      top_n(., nn, -p.adjust) %>%
      top_n(., nn, GeneRatio) %>%
      top_n(., nn, Count)
    barplot_df$Description <- factor(barplot_df$Description,
                                     levels = barplot_df$Description[order(barplot_df$Group,
                                                                           -barplot_df$p.adjust)])
    
    bar_plt <- ggplot(barplot_df)+
      geom_bar(aes(x = Description, y = -log10(p.adjust), fill = Group),
               stat = "identity") +
      geom_hline(yintercept = -log10(0.05), linetype = 2, size = 1, color = "gold")+
      scale_fill_manual(values = bar_col[unique(barplot_df$Group)]) +
      ylab(bquote("-"~Log[10]~"(adjusted P-value)"))+
      coord_flip() +
      facet_grid(Group~., scales = "free_y", space="free_y", drop = T)+
      theme_bw()+
      theme(legend.position = "none",
            strip.placement = "outside",
            strip.background = element_blank(),
            strip.text = element_text(size = 12, face = "bold", color = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 15, face = "bold", color = "black"),
            axis.text = element_text(size = 12, color = "black"))
    return(list("bar_plt" = bar_plt,
                "enrich_df" = enrich_df))
  }
}

