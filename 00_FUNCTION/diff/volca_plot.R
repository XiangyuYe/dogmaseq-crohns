library(ggplot2)
library(ggrepel)
library(dplyr)
library(cols4all)

# valco_col <- c4a("classic_blue_red12", 3)
# names(valco_col) <- c("Down regulation", "Others", "Up regulation")
volca_function <- function(deg_df = NULL,
                           term_col = "Gene",
                           fc_col = "log2FC",
                           p_col = "adj.P.Val",
                           target_gene = NULL,
                           n_top = 20,
                           log_thresh = 0.5, 
                           top_log_thresh = 0.5,
                           group_label = c("Down regulation", "Others", "Up regulation"),
                           valco_col = c4a("classic_blue_red12", 3),
                           highlight_col = "orange"){
  ## set deg groups
  deg_plt <- data.frame("Gene" = deg_df[[term_col]],
                        "log2FC" = deg_df[[fc_col]],
                        "adj.P.Val" = deg_df[[p_col]])
  deg_plt$adj.P.Val[deg_plt$adj.P.Val == 0] <- min(deg_plt$adj.P.Val[deg_plt$adj.P.Val > 0])
  deg_plt$DT <- ifelse(deg_plt$adj.P.Val >= 0.05 | abs(deg_plt$log2FC) < log_thresh, group_label[2],
                      ifelse(deg_plt$log2FC < 0, group_label[1], group_label[3])) 
  deg_plt$sig <- -log10(deg_plt$adj.P.Val)
  ## set highlight genes
  if (n_top <= 0 & is.null(target_gene)){
    
    target_df <- NULL
    
  } else {
    
    if(is.null(target_gene)){
      
      top_df <- group_by(deg_plt, DT) %>%
        dplyr::filter(., DT != "Others" & abs(log2FC) > top_log_thresh) %>%
        dplyr::top_n(., n_top, sig) %>%
        dplyr::top_n(., n_top, abs(log2FC))
      target_df <- top_df
      
    } else {
      
      target_df <- subset(deg_plt, deg_plt$Gene %in% target_gene)
      if (n_top > 0) {
        
        target_df <- target_df %>%
          dplyr::top_n(., n_top, sig) %>%
          dplyr::top_n(., n_top, abs(log2FC))
        
      }
      
      
    }
    
  }
  
  max_log2FC <- max(abs(deg_plt$log2FC))
  names(valco_col) <- group_label
  
  volca_plt <- ggplot(deg_plt)+
    geom_point(aes(x = log2FC, y = sig, color = DT),
               shape = 19, size = 1.5) + 
    geom_hline(yintercept = -log10(0.05), 
               linetype = 2, size = 1, color = "grey50")+
    geom_vline(xintercept = c(-log_thresh, 0, log_thresh), 
               linetype = 2, size = 1, color = c("grey80", "grey50", "grey80"))
  
  if (!is.null(target_df)) {
    volca_plt <- volca_plt + 
      geom_point(data = target_df, 
                 aes(x = log2FC, y = sig), 
                 color = highlight_col, shape = 19, size = 1.5) +
      geom_label_repel(data = target_df, 
                       aes(label = Gene, x = log2FC, y = sig), 
                       segment.size = 0.25, 
                       size = 2.5, 
                       max.iter = 1E7, 
                       min.segment.length = 10,
                       max.overlaps = getOption("ggrepel.max.overlaps", default = 200))
  }
  
  volca_plt <- volca_plt+
    scale_x_continuous(limits = c(-max_log2FC, max_log2FC))+
    scale_y_continuous(limits = c(0, max(deg_plt$sig) * 1.1))+
    scale_color_manual(values = valco_col[unique(deg_plt$DT)])+
    theme_bw()+
    xlab(bquote(~Log[2]~"(Fold Change)"))+
    # xlab(expression(Delta * "(" * beta * ")"))+
    ylab(bquote("-"~Log[10]~"(adjusted P-value)"))+
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 12, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face= "bold", color = "black"))
  
  return(volca_plt)
}

#####
bi.voca.plot <- function(deg_df1 = NULL,
                         deg_df2 = NULL,
                         name_df1 = "Data 1",
                         name_df2 = "Data 2",
                         p_col = "padj",
                         FC_col = "log2FC",
                         term_col = "Term",
                         fc_name = "Log2(Fold Change)",
                         target_gene = NULL,
                         n_top = 20,
                         p_thre = 0.05,
                         logFC_thre = 0.5,
                         col_use = c("#F02720", "#E9C39B", "#B5C8E2", "grey")){
  
  if (length(col_use) != 4) {
    col_use <- c("#F02720", "#E9C39B", "#B5C8E2", "grey")
  }
  names(col_use) <- c("Both", 
                      paste0("Significant in ", name_df1), 
                      paste0("Significant in ", name_df2), 
                      "Neither")
  # format
  deg_df1$Term <- deg_df1[[term_col]]
  deg_df1$padj <- deg_df1[[p_col]]
  deg_df1$log2FC <- deg_df1[[FC_col]]
  deg_df1$Expression <- ifelse(deg_df1$padj >= p_thre | abs(deg_df1$log2FC) < logFC_thre,
                               "Not Significant",
                               ifelse(deg_df1$log2FC > 0, "Up-regulated", "Down-regulated"))
  deg_df2$Term <- deg_df2[[term_col]]
  deg_df2$padj <- deg_df2[[p_col]]
  deg_df2$log2FC <- deg_df2[[FC_col]]
  deg_df2$Expression <- ifelse(deg_df2$padj >= p_thre | abs(deg_df2$log2FC) < logFC_thre,
                               "Not Significant",
                               ifelse(deg_df2$log2FC > 0, "Up-regulated", "Down-regulated"))
  ##
  deg_df1 <- deg_df1[, c("Term", "log2FC", "padj", "Expression")]
  deg_df2 <- deg_df2[, c("Term", "log2FC", "padj", "Expression")]
  
  colnames(deg_df1) <- paste0(colnames(deg_df1), "_1")
  colnames(deg_df2) <- paste0(colnames(deg_df2), "_2")
  
  deg_comb <- merge(deg_df1, 
                    deg_df2, 
                    by.x = "Term_1", 
                    by.y = "Term_2")
  deg_comb$Group <- ifelse(deg_comb$Expression_1 == "Not Significant" & 
                             deg_comb$Expression_2 == "Not Significant",
                           names(col_use)[4],
                           ifelse(deg_comb$Expression_1 != "Not Significant" & 
                                    deg_comb$Expression_2 != "Not Significant",
                                  names(col_use)[1],
                                  ifelse(deg_comb$Expression_1 != "Not Significant",
                                         names(col_use)[2],
                                         names(col_use)[3]))) %>% 
    factor(., levels = names(col_use))
  
  # deg_comb <- subset(deg_comb, Group != "Neither")
  top_df <- top_n(deg_comb[deg_comb$Group == "Both",], n_top, 
                  abs(log2FC_1 + log2FC_2))
  if (!is.null(target_gene)) {
    top_df <- rbind(top_df,
                    deg_comb[deg_comb$Term_1 %in% target_gene,])
  }
  #
  max_X <- max(abs(deg_comb$log2FC_1)) %>% ceiling()
  max_Y <- max(abs(deg_comb$log2FC_2)) %>% ceiling()
  use_color <- c("#F02720", "#E9C39B", "#B5C8E2", "grey")
  plt_inter <- ggplot(deg_comb) + 
    geom_point(aes(x = log2FC_1, y = log2FC_2, color = Group),
               shape = 19, size = 1.5) + 
    geom_hline(yintercept = 0, linetype = 2, size = 1, color = "grey20")+ 
    geom_vline(xintercept = 0, linetype = 2, size = 1, color = "grey20")+ 
    scale_x_continuous(limits = c(-max_X, max_X)) + 
    scale_y_continuous(limits = c(-max_Y, max_Y)) + 
    xlab(glue("{fc_name} in {name_df1}"))+
    ylab(glue("{fc_name} in {name_df2}"))+
    geom_label_repel(data = top_df, 
                     aes(label = Term_1, x = log2FC_1, y = log2FC_2), 
                     segment.size = 0.25, size = 3.5, 
                     max.iter = 1E7, min.segment.length = 2,
                     max.overlaps = getOption("ggrepel.max.overlaps", default = 100))+
    scale_colour_manual(values = use_color) + 
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.text = element_text(size = 12, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 15, face= "bold", color = "black"))
  
  return(plt_inter)
} 
##### number of diff #####
diff.num.plot <- function(diff_df = NULL,
                          cluster_col = "cluster",
                          p_adj_col = "padj",
                          fc_col = "log2FC",
                          label_x = "Number of DEGs",
                          p_adj_thresh = 0.05,
                          logfc_thresh = 0.25,
                          col_use = c("Down" = "#2C69B0", "Up" = "#F02720")){
  ## format df for number of DEGs
  if (!is.factor(diff_df[[cluster_col]])) {
    diff_df[[cluster_col]] <- as.factor(diff_df[[cluster_col]])
  }
  ndiff_list <- lapply(levels(diff_df[[cluster_col]]), function(ctx){
    
    diff_dfx <- subset(diff_df, diff_df[[cluster_col]] == ctx)
    n_up <- sum(diff_dfx[[p_adj_col]] < p_adj_thresh & diff_dfx[[fc_col]] > logfc_thresh, na.rm = T)
    n_down <- sum(diff_dfx[[p_adj_col]] < p_adj_thresh & diff_dfx[[fc_col]] < -logfc_thresh, na.rm = T)
    ndiff_dfx <- data.frame(celltype = ctx,
                            Direction = c("Up", "Down"),
                            ndiff = c(n_up, n_down))
    return(ndiff_dfx)
  }) %>% Reduce("rbind", .) %>% as.data.frame()
  # rev levels for plot
  ndiff_list$celltype <- factor(ndiff_list$celltype, 
                                levels = rev(levels(diff_df[[cluster_col]])))
  ## set parameters for plot
  n_max <- max(ndiff_list$ndiff)
  n_log <- floor(log10(n_max))
  n_max_up <- ceiling(n_max/10^n_log)*10^n_log
  ndiff_list$ndiff_db <- ifelse(ndiff_list$Direction == "Up", ndiff_list$ndiff, -ndiff_list$ndiff)
  ndiff_list$label_pos <- ifelse(ndiff_list$Direction == "Up", 
                                 ndiff_list$ndiff_db + n_max * 0.1, ndiff_list$ndiff_db - n_max * 0.1)
  ## plot
  plt <- ggplot(ndiff_list) + 
    geom_bar(aes(x = ndiff_db, 
                 y = celltype, 
                 fill = Direction),
             stat = "identity", width = 0.8) +
    xlab(label_x) + ylab("") + 
    geom_text(aes(x = label_pos, 
                  y = celltype, 
                  label = ndiff)) +
    annotate("text",
             x = -n_max * 0.4, y = 1,
             size = 5, colour="grey10", 
             label = "bold(Down)", parse = TRUE) +
    annotate("text",
             x = n_max * 0.4, y = 1,
             size = 5, colour="grey10",
             label = "bold(Up)", parse = TRUE) +
    scale_fill_manual(values = col_use) +
    scale_x_continuous(limits = c(-n_max * 1.3, n_max * 1.3),
                       breaks = seq(-n_max_up, n_max_up, n_max_up/4) %>% round() %>% unique,
                       labels = abs(seq(-n_max_up, n_max_up, n_max_up/4) %>% round() %>% unique)) +
    geom_vline(xintercept = 0, linetype = 1, 
               color = "black", linewidth = 0.8) + 
    theme_bw() + 
    theme(axis.title = element_text(size = 15, face = "bold"),
          axis.text = element_text(size = 12, color = "black"),
          legend.position = "none")
  #
  return(plt)
}
