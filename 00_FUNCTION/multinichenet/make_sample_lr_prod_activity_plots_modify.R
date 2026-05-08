## module load r/4.5.0
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(nichenetr)
library(multinichenetr)
library(glue)

##### make_circos_group_comparison_modify #####
make_circos_group_comparison_modify <- 
  function(prioritized_tbl_oi, 
           colors_sender, 
           colors_receiver,
           title_goi = NULL){
    
    requireNamespace("dplyr")
    requireNamespace("ggplot2")
    suppressMessages(require("circlize"))
    
    prioritized_tbl_oi <- prioritized_tbl_oi %>% 
      dplyr::ungroup() # if grouped: things will be messed up downstream
    # Link each cell type to a color
    grid_col_tbl_ligand <- tibble::tibble(sender = colors_sender %>% names(), 
                                          color_ligand_type = colors_sender)
    grid_col_tbl_receptor <- tibble::tibble(receiver = colors_receiver %>% names(), 
                                            color_receptor_type = colors_receiver)
    # Make the plot for each condition
    groups_oi <- prioritized_tbl_oi$group %>% unique()
    all_plots <- groups_oi %>% 
      lapply(function(group_oi){
        
        # Make the plot for condition of interest - title of the plot
        if (is.null(title_goi)) {
          title <- group_oi
        } else {
          title <- title_goi[group_oi]
        }
        circos_links_oi <- prioritized_tbl_oi %>% 
          dplyr::filter(group == group_oi)
        # deal with duplicated sector names
        # dplyr::rename the ligands so we can have the same ligand in multiple senders (and receptors in multiple receivers)
        # only do it with duplicated ones!
        # this is in iterative process because changing one ligand to try to solve a mistake here, can actually induce another problem
        circos_links <- circos_links_oi %>% 
          dplyr::rename(weight = prioritization_score)
        df <- circos_links
        ##
        ligand.uni <- unique(df$ligand)
        for (i in 1:length(ligand.uni)) {
          df.i <- df[df$ligand == ligand.uni[i], ]
          sender.uni <- unique(df.i$sender)
          for (j in 1:length(sender.uni)) {
            df.i.j <- df.i[df.i$sender == sender.uni[j], ]
            df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
            df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
          }
        }
        receptor.uni <- unique(df$receptor)
        for (i in 1:length(receptor.uni)) {
          df.i <- df[df$receptor == receptor.uni[i], ]
          receiver.uni <- unique(df.i$receiver)
          for (j in 1:length(receiver.uni)) {
            df.i.j <- df.i[df.i$receiver == receiver.uni[j], ]
            df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
            df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
          }
        }
        
        intersecting_ligands_receptors <- generics::intersect(unique(df$ligand),unique(df$receptor))
        
        while(length(intersecting_ligands_receptors) > 0){
          df_unique <- df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
          df_duplicated <- df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
          df_duplicated <- df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
          df <- dplyr::bind_rows(df_unique, df_duplicated)
          intersecting_ligands_receptors <- generics::intersect(unique(df$ligand),unique(df$receptor))
        }
        
        circos_links <- df
        
        # Link ligands/Receptors to the colors of senders/receivers
        circos_links <- circos_links %>% 
          dplyr::inner_join(grid_col_tbl_ligand) %>% 
          dplyr::inner_join(grid_col_tbl_receptor)
        links_circle <- circos_links %>% 
          dplyr::distinct(ligand,receptor, weight)
        ligand_color <- circos_links %>% 
          dplyr::distinct(ligand,color_ligand_type)
        grid_ligand_color <- ligand_color$color_ligand_type %>% 
          magrittr::set_names(ligand_color$ligand)
        receptor_color <- circos_links %>% 
          dplyr::distinct(receptor,color_receptor_type)
        grid_receptor_color <- receptor_color$color_receptor_type %>% 
          magrittr::set_names(receptor_color$receptor)
        grid_col <- c(grid_ligand_color,grid_receptor_color)
        # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
        transparency <- circos_links %>% 
          dplyr::mutate(weight = weight/max(weight)) %>% 
          dplyr::mutate(transparency = 1 - weight) %>% 
          .$transparency
        # Define order of the ligands and receptors and the gaps
        ligand_order <- prioritized_tbl_oi$sender %>% 
          unique() %>% 
          sort() %>% 
          lapply(function(sender_oi){
            ligands <- circos_links %>% 
              dplyr::filter(sender == sender_oi) %>%  
              dplyr::arrange(ligand) %>% 
              dplyr::distinct(ligand)
          }) %>% unlist()
        receptor_order <- prioritized_tbl_oi$receiver %>% 
          unique() %>% 
          sort() %>% 
          lapply(function(receiver_oi){
            receptors <- circos_links %>% 
              dplyr::filter(receiver == receiver_oi) %>%  
              dplyr::arrange(receptor) %>% 
              dplyr::distinct(receptor)
          }) %>% unlist()
        
        order <- c(ligand_order,receptor_order)
        
        width_same_cell_same_ligand_type <- 0.275
        width_different_cell <- 3
        width_ligand_receptor <- 9
        width_same_cell_same_receptor_type <- 0.275
        
        sender_gaps <- intersect(prioritized_tbl_oi$sender %>% unique(),
                                 circos_links$sender %>% unique()) %>% 
          sort() %>% 
          lapply(function(sender_oi){
            sector <- rep(width_same_cell_same_ligand_type, 
                          times = (circos_links %>% 
                                     dplyr::filter(sender == sender_oi) %>% 
                                     dplyr::distinct(ligand) %>% 
                                     nrow() -1))
            gap <- width_different_cell
            return(c(sector,gap))
          }) %>% unlist()
        sender_gaps <- sender_gaps[-length(sender_gaps)]
        
        receiver_gaps <- intersect(prioritized_tbl_oi$receiver %>% unique(),
                                   circos_links$receiver %>% unique()) %>% 
          sort() %>% 
          lapply(function(receiver_oi){
            sector <- rep(width_same_cell_same_receptor_type, 
                          times = (circos_links %>% 
                                     dplyr::filter(receiver == receiver_oi) %>% 
                                     dplyr::distinct(receptor) %>% 
                                     nrow() -1))
            gap <- width_different_cell
            return(c(sector,gap))
          }) %>% unlist()
        receiver_gaps <- receiver_gaps[-length(receiver_gaps)]
        
        gaps <- c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
        
        if(length(gaps) != length(union(circos_links$ligand, circos_links$receptor) %>% unique())){
          warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
        }
        
        links_circle$weight[links_circle$weight == 0] <- 0.01
        circos.clear()
        circos.par(gap.degree = gaps)
        chordDiagram(links_circle,
                     directional = 1,
                     order = order,
                     link.sort = TRUE,
                     link.decreasing = TRUE,
                     grid.col = grid_col,
                     transparency = transparency,
                     diffHeight = 0.0075,
                     direction.type = c("diffHeight", "arrows"),
                     link.visible = links_circle$weight > 0.01,
                     annotationTrack = "grid",
                     preAllocateTracks = list(track.height = 0.175),
                     grid.border = "gray35", 
                     link.arr.length = 0.05, 
                     link.arr.type = "big.arrow", 
                     link.lwd = 1.25, 
                     link.lty = 1, 
                     link.border="gray35",
                     reduce = 0,
                     scale = TRUE)
        circos.track(track.index = 1, panel.fun = function(x, y) {
          circos.text(CELL_META$xcenter, 
                      CELL_META$ylim[1], 
                      CELL_META$sector.index,
                      facing = "clockwise", 
                      niceFacing = TRUE, 
                      adj = c(0, 0.5), 
                      cex = 1)
        }, bg.border = NA) #
        
        title(title)
        p_circos <- recordPlot()
        return(p_circos)
        
      })
    names(all_plots) = groups_oi
    
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    # grid_col_all = c(colors_receiver, colors_sender)
    legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$receiver %>% unique() %>% sort(),
                                    type = "grid",
                                    legend_gp = grid::gpar(fill = colors_receiver[prioritized_tbl_oi$receiver %>% unique() %>% sort()]),
                                    title_position = "topleft",
                                    title = "Receiver")
    ComplexHeatmap::draw(legend, just = c("left", "bottom"))
    
    legend = ComplexHeatmap::Legend(at = prioritized_tbl_oi$sender %>% unique() %>% sort(),
                                    type = "grid",
                                    legend_gp = grid::gpar(fill = colors_sender[prioritized_tbl_oi$sender %>% unique() %>% sort()]),
                                    title_position = "topleft",
                                    title = "Sender")
    ComplexHeatmap::draw(legend, just = c("left", "top"))
    
    p_legend = grDevices::recordPlot()
    
    all_plots$legend = p_legend
    
    return(all_plots)
  }
##### make_sample_lr_prod_activity_plots_modify #####
make_sample_lr_prod_activity_plots_modify <- 
  function(prioritization_tables, 
           prioritized_tbl_oi, 
           widths = NULL){
    requireNamespace("dplyr")
    requireNamespace("ggplot2")
    sample_data = prioritization_tables$sample_prioritization_tbl %>% 
      dplyr::filter(id %in% prioritized_tbl_oi$id) %>% 
      dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                    lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
      dplyr::arrange(receiver) %>% 
      dplyr::group_by(receiver) %>% 
      dplyr::arrange(sender, .by_group = TRUE)
    sample_data = sample_data %>% 
      dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sample_data$sender_receiver %>% unique()))
    group_data = prioritization_tables$group_prioritization_table_source %>% 
      dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                    lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
      dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, 
                      group, activity, activity_scaled, direction_regulation, prioritization_score) %>% 
      dplyr::filter(id %in% sample_data$id) %>% 
      dplyr::arrange(receiver) %>% 
      dplyr::group_by(receiver) %>% 
      dplyr::arrange(sender, .by_group = TRUE)
    group_data = group_data %>% 
      dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
    group_data_celltype_specificity = prioritization_tables$group_prioritization_tbl %>% 
      dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                    lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
      dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, 
                      group, scaled_pb_ligand, scaled_pb_receptor) %>% 
      dplyr::filter(id %in% sample_data$id) %>% 
      dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>% 
      dplyr::arrange(sender, .by_group = TRUE)
    group_data_celltype_specificity = group_data_celltype_specificity %>% 
      dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data_celltype_specificity$sender_receiver %>% unique()))
    group_data_frac_expression = prioritization_tables$group_prioritization_table_source %>% 
      dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                    lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
      dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, 
                      group, fraction_ligand_group, fraction_receptor_group) %>% 
      dplyr::filter(id %in% sample_data$id) %>% 
      dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>% 
      dplyr::arrange(sender, .by_group = TRUE)
    group_data_frac_expression = group_data_frac_expression %>% 
      dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
    group_data = group_data %>% inner_join(group_data_celltype_specificity) %>% 
      inner_join(group_data_frac_expression)
    group_data = group_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
    rm(group_data_celltype_specificity)
    rm(group_data_frac_expression)
    keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
    names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
    
    ## aggregate
    sample_data_agg <- lapply(unique(sample_data$group), function(gx){
      
      sample_datax <- subset(sample_data, group == gx)
      lapply(unique(sample_datax$id), function(idx){
        
        sample_dataxx <- subset(sample_datax, id == idx)
        data.frame(group = gx, 
                   lr_interaction = unique(sample_dataxx$lr_interaction),
                   scaled_LR_pb_prod = mean(sample_dataxx$scaled_LR_pb_prod),
                   prop_sender_receiver = sum(sample_dataxx$keep_sender_receiver == "Sender & Receiver present")/
                     nrow(sample_dataxx),
                   sender_receiver = unique(sample_dataxx$sender_receiver))
        
      }) %>% Reduce("rbind", .)
      
    }) %>% Reduce("rbind", .)
    
    p1 = sample_data_agg %>% 
      ggplot(aes(group, 
                 lr_interaction, 
                 color = scaled_LR_pb_prod, 
                 size = prop_sender_receiver)) + 
      geom_point() + 
      facet_grid(sender_receiver ~ ., 
                 scales = "free", 
                 space = "free", 
                 switch = "y") + 
      scale_x_discrete(position = "top") + 
      theme_light() + 
      theme(axis.ticks = element_blank(), 
            axis.title = element_blank(), 
            axis.text.y = element_text(face = "bold.italic", size = 9), 
            axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.key.size = unit(0.2, "inches"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.spacing.x = unit(0.4, "lines"), 
            panel.spacing.y = unit(0.25, "lines"),
            strip.text.x.top = element_text(size = 10, color = "black", face = "bold", angle = 0), 
            strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0), 
            strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
      labs(color = "Scaled L-R\npseudobulk exprs product", 
           size = "Proportion of sufficient\n sender & receiver presence")
    max_lfc = abs(sample_data_agg$scaled_LR_pb_prod) %>% max()
    custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(), 
                                              values = c(0, 0.35, 0.485, 0.5, 0.515, 0.65, 1), 
                                              limits = c(-1 * max_lfc, max_lfc))
    p1 = p1 + custom_scale_fill
    p2 = group_data %>% 
      ggplot(aes(direction_regulation, lr_interaction, fill = activity_scaled)) + 
      geom_tile(color = "whitesmoke") + 
      facet_grid(sender_receiver ~ group, 
                 scales = "free", 
                 space = "free") + 
      scale_x_discrete(position = "top") + 
      theme_light() + 
      theme(axis.ticks = element_blank(), 
            axis.title = element_blank(), 
            axis.text.y = element_text(face = "bold.italic", size = 9), 
            axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.key.size = unit(0.2, "inches"),
            strip.text.x.top = element_text(angle = 0), 
            # strip.text.y = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.spacing.x = unit(0.2, "lines"), 
            panel.spacing.y = unit(0.25, "lines"), 
            strip.text.x = element_text(size = 10, color = "black", face = "bold"), 
            strip.text.y = element_blank(), 
            strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
      labs(fill = "Scaled Ligand\nActivity in Receiver")
    max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
    custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]), 
                                             values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1), 
                                             limits = c(-1 * max_activity, max_activity))
    p2 = p2 + custom_scale_fill
    p3 = group_data %>% 
      ggplot(aes(direction_regulation, 
                 lr_interaction, 
                 fill = activity)) + 
      geom_tile(color = "whitesmoke") + 
      facet_grid(sender_receiver ~ group, 
                 scales = "free", 
                 space = "free") + 
      scale_x_discrete(position = "top") + 
      theme_light() + theme(axis.ticks = element_blank(), 
                            axis.title = element_blank(), 
                            axis.title.y = element_blank(), 
                            axis.text.y = element_blank(), 
                            axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
                            legend.text = element_text(size = 8),
                            legend.title = element_text(size = 8),
                            legend.key.size = unit(0.2, "inches"),
                            strip.text.x.top = element_text(angle = 0), 
                            # strip.text.y = element_blank(), 
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                            panel.spacing.x = unit(0.2, "lines"), panel.spacing.y = unit(0.25, "lines"), 
                            strip.text.x = element_text(size = 10, color = "black", face = "bold"), 
                            strip.text.y = element_blank(), 
                            strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
      labs(fill = "Ligand\nActivity in Receiver")
    max_activity = (group_data$activity) %>% max()
    min_activity = (group_data$activity) %>% min()
    custom_scale_fill = scale_fill_gradient2(low = "white", 
                                             mid = "white", 
                                             high = "darkorange", 
                                             midpoint = 0)
    p3 = p3 + custom_scale_fill
    cs_data = group_data %>% 
      distinct(sender_receiver, lr_interaction, 
               group, scaled_pb_ligand, scaled_pb_receptor) %>% 
      tidyr::gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
    cs_data$LR[cs_data$LR == "scaled_pb_ligand"] = "ligand"
    cs_data$LR[cs_data$LR == "scaled_pb_receptor"] = "receptor"
    frac_data = group_data %>% 
      distinct(sender_receiver, lr_interaction, 
               group, fraction_ligand_group, fraction_receptor_group) %>% 
      tidyr::gather(LR, fraction_expression, fraction_ligand_group:fraction_receptor_group)
    frac_data$LR[frac_data$LR == "fraction_ligand_group"] = "ligand"
    frac_data$LR[frac_data$LR == "fraction_receptor_group"] = "receptor"
    cs_data = cs_data %>% inner_join(frac_data)
    p_cs = cs_data %>% 
      ggplot(aes(LR, 
                 lr_interaction, 
                 color = celltype_specificity, 
                 size = fraction_expression)) + 
      geom_point() + 
      facet_grid(sender_receiver ~ group, scales = "free", space = "free") + 
      scale_x_discrete(position = "top") + 
      theme_light() + 
      viridis::scale_color_viridis() + 
      theme(axis.ticks = element_blank(), 
            axis.title = element_blank(), 
            axis.title.y = element_blank(), 
            axis.text.y = element_blank(), 
            axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.key.size = unit(0.2, "inches"),
            strip.text.x.top = element_text(angle = 0), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.spacing.x = unit(0.2, "lines"), 
            panel.spacing.y = unit(0.25, "lines"), 
            strip.text.x = element_text(size = 10, color = "black", face = "bold"), 
            strip.text.y = element_blank(), 
            strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
      labs(color = "Scaled celltype specificity") + labs(size = "Fraction of expression")
    if (!is.null(widths)) {
      p = patchwork::wrap_plots(p1, 
                                p2, 
                                # p3, 
                                p_cs, nrow = 1, 
                                guides = "collect", widths = widths)
    }
    else {
      p = patchwork::wrap_plots(p1, 
                                p2, 
                                # p3, 
                                p_cs, nrow = 1, 
                                guides = "collect", 
                                widths = c(sample_data$group %>% unique() %>% length(), 
                                           2 * (sample_data$group %>% unique() %>% length()), 
                                           # 2 * (sample_data$group %>% unique() %>% length()), 
                                           2 * (sample_data$group %>% unique() %>% length())))
    }
    return(p)
  }
##### make_DEgene_dotplot_pseudobulk_reversed_modify #####
make_DEgene_dotplot_pseudobulk_reversed_modify <- 
  function (genes_oi, 
            celltype_info, 
            prioritization_tables, 
            celltype_oi, 
            grouping_tbl, 
            groups_oi = NULL, 
            target_regulation_df = NULL){
    ##
    requireNamespace("dplyr")
    requireNamespace("ggplot2")
    if (is.null(target_regulation_df)) {
      keep_tbl = prioritization_tables$sample_prioritization_tbl %>% 
        dplyr::distinct(sample, group, receiver, keep_receiver) %>% 
        dplyr::rename(celltype = receiver) %>% 
        dplyr::mutate(keep_receiver = as.logical(keep_receiver))
      keep_sender_receiver_values = c(1, 4)
      names(keep_sender_receiver_values) = c(FALSE, TRUE)
      plot_data = celltype_info$pb_df %>% 
        dplyr::inner_join(grouping_tbl, by = c("sample")) %>% 
        dplyr::inner_join(keep_tbl, by = c("sample", "group", "celltype"))
      plot_data = plot_data %>% 
        dplyr::group_by(gene, celltype) %>% 
        dplyr::mutate(scaled_gene_exprs = nichenetr::scaling_zscore(pb_sample)) %>% 
        dplyr::ungroup()
      plot_data$gene = factor(plot_data$gene, levels = genes_oi)
      plot_data = plot_data %>% 
        dplyr::filter(gene %in% genes_oi & 
                        celltype %in% celltype_oi)
      if (!is.null(groups_oi)) {
        plot_data = plot_data %>% 
          dplyr::filter(group %in% groups_oi)
      }
      p1 = plot_data %>% 
        ggplot(aes(gene, sample, 
                   color = scaled_gene_exprs, 
                   size = keep_receiver)) + 
        geom_point() + 
        facet_grid(group ~ ., scales = "free", space = "free") + 
        scale_x_discrete(position = "top") + 
        theme_light() + 
        theme(axis.ticks = element_blank(), 
              axis.title.x = element_text(size = 10), 
              axis.title.y = element_text(size = 10), 
              axis.text.x = element_text(face = "italic", size = 9, angle = 90, hjust = 0), 
              axis.text.y = element_text(size = 9), 
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.key.size = unit(0.2, "inches"),
              strip.text.y.left = element_text(angle = 0), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.spacing.y = unit(1, "lines"), 
              panel.spacing.x = unit(0.5, "lines"), 
              strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0), 
              strip.text.x = element_text(size = 9, color = "black", face = "bold"), 
              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(color = "Scaled pseudobulk\nexpression", size = "Celltype present") + 
        xlab("Genes") + 
        ylab("Samples") + 
        scale_size_manual(values = keep_sender_receiver_values)
      max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
      custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(), 
                                                values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1), 
                                                limits = c(-1 * max_lfc, max_lfc))
      p1 = p1 + custom_scale_fill
      frq_df = celltype_info$frq_df %>% 
        dplyr::filter(gene %in% genes_oi & 
                        celltype %in% celltype_oi)
      plot_data = plot_data %>% 
        dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype", "group"))
      plot_data$gene = factor(plot_data$gene, levels = genes_oi)
      #
      p2 = plot_data %>% 
        ggplot(aes(gene, sample, 
                   color = scaled_gene_exprs, 
                   size = fraction_sample)) + 
        geom_point() + 
        facet_grid(group ~ ., scales = "free", space = "free") + 
        scale_x_discrete(position = "top") + 
        theme_light() + 
        theme(axis.ticks = element_blank(), 
              axis.title.x = element_text(size = 10), 
              axis.title.y = element_text(size = 10), 
              axis.text.x = element_text(face = "italic", size = 9, angle = 90, hjust = 0), 
              axis.text.y = element_text(size = 9), 
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.key.size = unit(0.2, "inches"),
              strip.text.y.left = element_text(angle = 0), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.spacing.y = unit(1, "lines"), 
              panel.spacing.x = unit(0.5, "lines"), 
              strip.text.y = element_text(size = 9, color = "black", face = "bold", angle = 0), 
              strip.text.x = element_text(size = 9, color = "black", face = "bold"), 
              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(color = "Scaled pseudobulk\nexpression", 
             size = "Fraction of\nexpressing cells") + 
        xlab("Genes") + 
        ylab("Samples")
      max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
      custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(), 
                                                values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1), 
                                                limits = c(-1 * max_lfc, max_lfc))
      p2 = p2 + custom_scale_fill
    } else {
      keep_tbl = prioritization_tables$sample_prioritization_tbl %>% 
        dplyr::distinct(sample, group, receiver, keep_receiver) %>% 
        dplyr::rename(celltype = receiver) %>% 
        dplyr::mutate(keep_receiver = as.logical(keep_receiver))
      keep_sender_receiver_values = c(1, 4)
      names(keep_sender_receiver_values) = c(FALSE, TRUE)
      plot_data = celltype_info$pb_df %>% 
        dplyr::inner_join(grouping_tbl, by = c("sample")) %>% 
        dplyr::inner_join(keep_tbl, by = c("sample", "group", "celltype"))
      plot_data = plot_data %>% 
        dplyr::group_by(gene, celltype) %>% 
        dplyr::mutate(scaled_gene_exprs = nichenetr::scaling_zscore(pb_sample)) %>% 
        dplyr::ungroup()
      plot_data$gene = factor(plot_data$gene, levels = genes_oi)
      plot_data = plot_data %>% 
        dplyr::filter(gene %in% genes_oi & 
                        celltype %in% celltype_oi) %>% 
        dplyr::inner_join(target_regulation_df) %>% 
        dplyr::mutate(gene = factor(gene, levels = genes_oi))
      if (!is.null(groups_oi)) {
        plot_data = plot_data %>% 
          dplyr::filter(group %in% groups_oi)
      }
      ##
      p1 = plot_data %>% 
        ggplot(aes(gene, sample, 
                   color = scaled_gene_exprs, 
                   size = keep_receiver)) + 
        geom_point() + 
        facet_grid(group ~ direction_regulation, scales = "free", space = "free") + 
        scale_x_discrete(position = "top") + 
        theme_light() + 
        theme(axis.ticks = element_blank(), 
              axis.title.x = element_text(size = 10), 
              axis.title.y = element_text(size = 10), 
              axis.text.x = element_text(face = "italic", size = 9, angle = 90, hjust = 0), 
              axis.text.y = element_text(size = 9), 
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.key.size = unit(0.2, "inches"),
              strip.text.y.left = element_text(angle = 0), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.spacing.y = unit(1, "lines"), 
              panel.spacing.x = unit(0.5, "lines"), 
              strip.text.y = element_text(size = 9, color = "black", angle = 0), 
              strip.text.x = element_text(size = 9, color = "black"), 
              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(color = "Scaled pseudobulk\nexpression", 
             size = "Celltype present") + 
        xlab("Genes") + 
        ylab("Samples") + 
        scale_size_manual(values = keep_sender_receiver_values)
      max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
      custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(), 
                                                values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1), 
                                                limits = c(-1 * max_lfc, max_lfc))
      p1 = p1 + custom_scale_fill
      frq_df = celltype_info$frq_df %>% 
        dplyr::filter(gene %in% genes_oi & 
                        celltype %in% celltype_oi)
      plot_data = plot_data %>% 
        dplyr::inner_join(frq_df, by = c("gene", "sample", "celltype", "group"))
      plot_data$gene = factor(plot_data$gene, levels = genes_oi)
      p2 = plot_data %>% 
        ggplot(aes(gene, sample, 
                   color = scaled_gene_exprs, 
                   size = fraction_sample)) + 
        geom_point() + 
        facet_grid(group ~ direction_regulation, scales = "free", space = "free") + 
        scale_x_discrete(position = "top") + 
        theme_light() + 
        theme(axis.ticks = element_blank(), 
              axis.title.x = element_text(size = 10), 
              axis.title.y = element_text(size = 10), 
              axis.text.x = element_text(face = "italic", size = 9, angle = 90, hjust = 0), 
              axis.text.y = element_text(size = 9), 
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.key.size = unit(0.2, "inches"),
              strip.text.y.left = element_text(angle = 0), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.spacing.y = unit(1, "lines"), 
              panel.spacing.x = unit(0.5, "lines"), 
              strip.text.x = element_text(size = 9, color = "black", face = "bold"), 
              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(color = "Scaled pseudobulk\nexpression", 
             size = "Fraction of\nexpressing cells") + 
        xlab("Genes") + 
        ylab("Samples")
      max_lfc = abs(plot_data$scaled_gene_exprs) %>% max()
      custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(), 
                                                values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1), 
                                                limits = c(-1 * max_lfc, max_lfc))
      p2 = p2 + custom_scale_fill
      ## aggregate
      plot_data_agg <- lapply(unique(plot_data$group), function(gx){
        
        plot_datax <- subset(plot_data, group == gx)
        lapply(unique(plot_datax$gene), function(genex){
          
          plot_dataxx <- subset(plot_datax, gene == genex)
          data.frame(gene = genex,
                     group = gx, 
                     pb_group = mean(plot_dataxx$pb_sample),
                     celltype = unique(plot_dataxx$celltype),
                     scaled_gene_exprs = mean(plot_dataxx$scaled_gene_exprs),
                     fraction_sample = mean(plot_dataxx$fraction_sample),
                     direction_regulation = unique(plot_dataxx$direction_regulation))
          
        }) %>% Reduce("rbind", .)
        
      }) %>% Reduce("rbind", .)
      p3 = plot_data_agg %>% 
        ggplot(aes(gene, group, 
                   color = scaled_gene_exprs, 
                   size = fraction_sample)) + 
        geom_point() + 
        facet_grid(. ~ direction_regulation, scales = "free", space = "free") + 
        scale_x_discrete(position = "top") + 
        theme_light() + 
        theme(axis.ticks = element_blank(), 
              axis.title.x = element_text(size = 10), 
              axis.title.y = element_text(size = 10), 
              axis.text.x = element_text(face = "italic", size = 9, angle = 90, hjust = 0), 
              axis.text.y = element_text(size = 9), 
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.key.size = unit(0.2, "inches"),
              strip.text.y.left = element_text(angle = 0), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.spacing.y = unit(1, "lines"), 
              panel.spacing.x = unit(0.5, "lines"), 
              strip.text.x = element_text(size = 9, color = "black", face = "bold"), 
              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(color = "Scaled pseudobulk\nexpression", 
             size = "Average fraction of\nexpressing cells") + 
        xlab("Genes") + 
        ylab("Groups")
      max_lfc = abs(plot_data_agg$scaled_gene_exprs) %>% max()
      custom_scale_fill = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(), 
                                                values = c(0, 0.35, 0.465, 0.5, 0.535, 0.65, 1), 
                                                limits = c(-1 * max_lfc, max_lfc))
      p3 = p3 + custom_scale_fill
    }
    return(list(pseudobulk_plot = p1, singlecell_plot = p2, group_plot = p3))
  }
##### make_ligand_activity_target_plot_modify #####
make_ligand_activity_target_plot_modify <- 
  function (group_oi, 
            receiver_oi, 
            prioritized_tbl_oi, 
            prioritization_tables, 
            ligand_activities_targets_DEgenes, 
            contrast_tbl, 
            grouping_tbl, 
            receiver_info, 
            ligand_target_matrix, 
            groups_oi = NULL, 
            plot_legend = TRUE, 
            heights = NULL, 
            widths = NULL) 
  {
    requireNamespace("dplyr")
    requireNamespace("ggplot2")
    if (is.null(groups_oi)) {
      groups_oi = contrast_tbl %>% 
        dplyr::pull(group) %>% 
        unique()
    }
    best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()
    active_ligand_target_links_df = ligand_activities_targets_DEgenes$ligand_activities %>% 
      dplyr::ungroup() %>% 
      dplyr::inner_join(contrast_tbl) %>% 
      dplyr::filter(ligand %in% best_upstream_ligands & 
                      receiver == receiver_oi & 
                      group == group_oi) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(ligand, target, ligand_target_weight, 
                    direction_regulation) %>% 
      dplyr::rename(weight = ligand_target_weight)
    active_ligand_target_links_df = active_ligand_target_links_df %>% 
      dplyr::filter(!is.na(weight))
    if (nrow(active_ligand_target_links_df) == 0) {
      message("No active ligand_target links retained!")
      return(NULL)
    } else {
      if (active_ligand_target_links_df$target %>% 
          unique() %>% 
          length() <= 2) {
        cutoff = 0
      } else {
        cutoff = 0.2
      }
    }
    active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                                ligand_target_matrix = ligand_target_matrix, 
                                                                                cutoff = cutoff)
    order_ligands_ = generics::intersect(best_upstream_ligands, 
                                         colnames(active_ligand_target_links)) %>% rev()
    order_targets_ = active_ligand_target_links_df$target %>% 
      unique() %>% 
      generics::intersect(rownames(active_ligand_target_links))
    order_ligands = order_ligands_ %>% make.names()
    order_targets = order_targets_ %>% make.names()
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% 
      make.names()
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% 
      make.names()
    if (!is.matrix(active_ligand_target_links[order_targets, order_ligands])) {
      vis_ligand_target = active_ligand_target_links[order_targets, order_ligands] %>% matrix(ncol = length(order_targets))
      rownames(vis_ligand_target) = order_ligands
      colnames(vis_ligand_target) = order_targets
    } else {
      vis_ligand_target = active_ligand_target_links[order_targets, order_ligands] %>% t()
    }
    vis_ligand_target_df = vis_ligand_target %>% data.frame() %>% 
      tibble::rownames_to_column("ligand") %>% 
      tidyr::gather("target", "score", -ligand) %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(ligand = factor(ligand, levels = order_ligands)) %>% 
      dplyr::inner_join(active_ligand_target_links_df %>% distinct(target, direction_regulation)) %>% 
      dplyr::mutate(target = factor(target, levels = order_targets))
    p_ligand_target_network = vis_ligand_target_df %>% 
      ggplot(aes(target, ligand, fill = score)) + 
      geom_tile(color = "whitesmoke", size = 0.5) + 
      facet_grid(. ~ direction_regulation, scales = "free", space = "free") + 
      scale_fill_gradient2(low = "white", mid = "purple", high = "darkred", midpoint = 0.14) + 
      theme_light() + 
      scale_x_discrete(position = "top") + 
      theme(axis.ticks = element_blank(), 
            axis.title.x = element_text(size = 10), 
            axis.title.y = element_text(size = 10), 
            axis.text.y = element_text(size = 9), 
            axis.text.x = element_text(size = 9, angle = 90, hjust = 0, face = "italic"), 
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.key.size = unit(0.2, "inches"),
            strip.text.x.top = element_text(angle = 0), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.spacing.x = unit(0.5, "lines"), 
            strip.text.x = element_text(size = 9, color = "black"), 
            strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
      labs(fill = "Regulatory Potential") + xlab("Predicted target genes") + 
      ylab("Prioritized ligands")
    custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 11, name = "PiYG") %>% .[1:5] %>% rev()), 
                                             values = c(0, 0.04, 0.12, 0.3, 0.4, 0.55, 1), 
                                             limits = c(0, max(ligand_target_matrix)))
    p_ligand_target_network = p_ligand_target_network + custom_scale_fill
    ligand_activity_df = ligand_activities_targets_DEgenes$ligand_activities %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(ligand %in% order_ligands_ & 
                      receiver == receiver_oi) %>% 
      dplyr::inner_join(contrast_tbl) %>% 
      dplyr::filter(group %in% groups_oi) %>% 
      dplyr::select(ligand, group, direction_regulation, activity_scaled) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(ligand = factor(ligand, levels = order_ligands))
    p_ligand_activity_scaled = ligand_activity_df %>% 
      ggplot(aes(direction_regulation, ligand, fill = activity_scaled)) + 
      geom_tile(color = "whitesmoke") + 
      facet_grid(. ~ group, scales = "free", space = "free") + 
      scale_x_discrete(position = "top") + theme_light() + 
      theme(axis.ticks = element_blank(), 
            axis.title.x = element_text(size = 10), 
            axis.title.y = element_text(size = 10), 
            axis.text.y = element_text(size = 9), 
            axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.key.size = unit(0.2, "inches"),
            strip.text.x.top = element_text(angle = 0), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.spacing.x = unit(0.5, "lines"), 
            panel.spacing.y = unit(0.5, "lines"), 
            strip.text.x = element_text(size = 10, color = "black"), 
            strip.text.y = element_blank(), 
            strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
      labs(fill = "Scaled Ligand\nActivity in Receiver") + 
      ylab("Prioritized ligands") + xlab("Scaled ligand activity")
    max_activity = abs(ligand_activity_df$activity_scaled) %>% 
      max(na.rm = TRUE)
    custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]), 
                                             values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1), 
                                             limits = c(-1 * max_activity, max_activity))
    p_ligand_activity_scaled = p_ligand_activity_scaled + custom_scale_fill
    ligand_activity_df = ligand_activities_targets_DEgenes$ligand_activities %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(ligand %in% order_ligands_ & 
                      receiver == receiver_oi) %>% 
      dplyr::inner_join(contrast_tbl) %>% 
      dplyr::filter(group %in% groups_oi) %>% 
      dplyr::select(ligand, group, direction_regulation, activity) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(ligand = factor(ligand, levels = order_ligands))
    p_ligand_activity = ligand_activity_df %>% 
      ggplot(aes(direction_regulation, ligand, fill = activity)) + 
      geom_tile(color = "whitesmoke") + 
      facet_grid(. ~ group, scales = "free", space = "free") + 
      scale_x_discrete(position = "top") + 
      theme_light() + 
      theme(axis.ticks = element_blank(), 
            axis.title.x = element_text(size = 10), 
            axis.title.y = element_blank(), 
            axis.text.y = element_text(size = 9), 
            axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8),
            legend.key.size = unit(0.2, "inches"),
            strip.text.x.top = element_text(angle = 0), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.spacing.x = unit(0.5, "lines"), 
            panel.spacing.y = unit(0.5, "lines"), 
            strip.text.x = element_text(size = 10, color = "black"), 
            strip.text.y = element_blank(), 
            strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
      labs(fill = "Ligand Activity\nin Receiver") + ylab("Prioritized ligands") + 
      xlab("Ligand activity")
    custom_scale_fill = scale_fill_gradient2(low = "white", 
                                             mid = "white", 
                                             high = "darkorange", 
                                             midpoint = 0)
    p_ligand_activity = p_ligand_activity + custom_scale_fill
    target_regulation_df = ligand_activities_targets_DEgenes$ligand_activities %>% 
      dplyr::ungroup() %>% 
      dplyr::inner_join(contrast_tbl) %>% 
      dplyr::filter(ligand %in% best_upstream_ligands & 
                      receiver == receiver_oi & group == group_oi) %>% 
      dplyr::ungroup() %>% 
      dplyr::distinct(target, direction_regulation) %>% dplyr::rename(gene = target)
    p_targets = make_DEgene_dotplot_pseudobulk_reversed_modify(genes_oi = order_targets_, 
                                                               celltype_info = receiver_info, 
                                                               prioritization_tables = prioritization_tables, 
                                                               celltype_oi = receiver_oi, 
                                                               grouping_tbl = grouping_tbl, 
                                                               groups_oi = groups_oi, 
                                                               target_regulation_df = target_regulation_df)
    n_groups = ligand_activity_df$group %>% unique() %>% length()
    n_targets = ncol(vis_ligand_target)
    n_ligands = nrow(vis_ligand_target)
    n_samples = grouping_tbl %>% 
      dplyr::filter(group %in% groups_oi) %>% 
      dplyr::pull(sample) %>% length()
    legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_activity_scaled)), 
                                    # ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_activity)), 
                                    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), 
                                    ncol = 1) %>% 
      patchwork::wrap_plots(., ggpubr::as_ggplot(ggpubr::get_legend(p_targets$group_plot)),
                            ncol = 1)
    if (is.null(heights)) {
      heights = c(n_ligands + 3, n_groups + 2)
    }
    if (is.null(widths)) {
      widths = c(n_groups * 2 + 0.75, 
                 # n_groups * 2, 
                 n_targets)
    }
    if (plot_legend == FALSE) {
      # design <- "AaB\n               ##C"
      design <- "AB\n               #C"
      combined_plot = patchwork::wrap_plots(A = p_ligand_activity_scaled + 
                                              theme(legend.position = "none", 
                                                    axis.ticks = element_blank()) + 
                                              theme(axis.title.x = element_text()), 
                                            # a = p_ligand_activity + 
                                            #   theme(legend.position = "none", 
                                            #         axis.ticks = element_blank()) + 
                                            #   ylab(""), 
                                            B = p_ligand_target_network + 
                                              theme(legend.position = "none", 
                                                    axis.ticks = element_blank()) + 
                                              ylab(""), 
                                            C = p_targets$group_plot + 
                                              theme(legend.position = "none") + 
                                              xlab(""), 
                                            nrow = 2, 
                                            design = design, 
                                            widths = widths, 
                                            heights = heights)
      return(list(combined_plot = combined_plot, legends = legends, 
                  "dim_lt" = c("n_ligand" = length(order_ligands), 
                               "n_target" = length(order_targets))))
    } else {
      # design <- "AaB\n               L#C"
      design <- "AB\n               LC"
      combined_plot = patchwork::wrap_plots(A = p_ligand_activity_scaled + 
                                              theme(legend.position = "none", 
                                                    axis.ticks = element_blank()) + 
                                              theme(axis.title.x = element_text()), 
                                            # a = p_ligand_activity + 
                                            #   theme(legend.position = "none", 
                                            #         axis.ticks = element_blank()) + 
                                            #   ylab(""), 
                                            B = p_ligand_target_network + 
                                              theme(legend.position = "none", 
                                                    axis.ticks = element_blank()) + 
                                              ylab(""), 
                                            C = p_targets$group_plot + 
                                              theme(legend.position = "none") + 
                                              xlab(""), 
                                            L = legends, 
                                            nrow = 2, 
                                            design = design, 
                                            widths = widths, 
                                            heights = heights)
      return(list(combined_plot = combined_plot, legends = legends, 
                  "dim_lt" = c("n_ligand" = length(order_ligands), 
                             "n_target" = length(order_targets))))
    }
  }
##### make_lr_target_correlation_plot_modify #####
make_lr_target_correlation_plot_modify <- 
  function(prioritization_tables, 
           prioritized_tbl_oi, 
           lr_target_prior_cor_filtered, 
           grouping_tbl, 
           receiver_info, 
           receiver_oi, 
           plot_legend = TRUE, 
           heights = NULL, 
           widths = NULL) 
  {
    requireNamespace("dplyr")
    requireNamespace("ggplot2")
    lr_target_prior_cor_filtered = lr_target_prior_cor_filtered %>% 
      dplyr::mutate(pearson = abs(pearson), 
                    direction_regulation = factor(direction_regulation, levels = c("up", "down")))
    sample_data = prioritization_tables$sample_prioritization_tbl %>% 
      dplyr::filter(id %in% prioritized_tbl_oi$id) %>% 
      dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), 
                    lr_interaction = paste(ligand, receptor, sep = " - ")) %>% 
      dplyr::arrange(sender)
    group_data = prioritization_tables$group_prioritization_table_source %>% 
      dplyr::distinct(id, group) %>% 
      dplyr::filter(id %in% sample_data$id)
    keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
    names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)
    sample_data = sample_data %>% 
      dplyr::mutate(sample = factor(sample)) %>% 
      dplyr::mutate(sample = factor(sample, levels = rev(levels(sample))))
    lr_target_data = lr_target_prior_cor_filtered %>% 
      dplyr::inner_join(sample_data, by = c("sender", "receiver", 
                                            "ligand", "receptor", 
                                            "id", "group")) %>% 
      dplyr::select(id, lr_interaction, sender, receiver, sender_receiver, 
                    group, target, prior_score, pearson, direction_regulation) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "))
    lr_target_data = lr_target_data %>% 
      dplyr::mutate(target = factor(target))
    order_targets = lr_target_data$target %>% 
      levels()
    lr_target_data = lr_target_data %>% 
      dplyr::mutate(lr_interaction = factor(lr_interaction))
    order_lr_interaction = lr_target_data$lr_interaction %>% 
      levels()
    if (nrow(lr_target_data) == 0) {
      return(NULL)
    } else {
      p2 = lr_target_data %>% 
        ggplot(aes(target, lr_interaction, 
                   color = prior_score, 
                   size = pearson)) + 
        geom_point() + 
        facet_grid(sender_receiver ~ direction_regulation, 
                   scales = "free", 
                   space = "free") + 
        scale_x_discrete(position = "top") + 
        theme_light() + 
        theme(axis.ticks = element_blank(), 
              axis.text.y = element_text(face = "bold.italic", size = 9), 
              axis.text.x = element_text(size = 8, angle = 90, hjust = 0, face = "italic"), 
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.key.size = unit(0.2, "inches"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.spacing.x = unit(0.5, "lines"), 
              panel.spacing.y = unit(0.25, "lines"), 
              strip.text.x.top = element_text(size = 9, color = "black", angle = 0), 
              # strip.text.y.right = element_text(size = 8, color = "black", angle = 0),
              strip.text.y.right = element_blank(),
              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(fill = "Prior knowledge score\nof correlated targets", 
             size = "Correlation of LR \nand target expression")
      custom_scale_color = scale_color_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 11, name = "PiYG") %>% .[1:5] %>% rev()), 
                                                 values = c(0, 0.0125, 0.004, 0.175, 0.35, 0.5, 1), 
                                                 limits = c(0, max(lr_target_prior_cor_filtered$prior_score)))
      p_lr_target = p2 + custom_scale_color + 
        xlab("Correlated target genes") + ylab("")
      groups_oi = group_data %>% 
        dplyr::pull(group) %>% 
        unique()
      target_regulation_df = lr_target_prior_cor_filtered %>% 
        dplyr::distinct(target, direction_regulation) %>% 
        dplyr::rename(gene = target)
      p_target_exprs = make_DEgene_dotplot_pseudobulk_reversed_modify(genes_oi = order_targets, 
                                                                      celltype_info = receiver_info, 
                                                                      prioritization_tables = prioritization_tables, 
                                                                      celltype_oi = receiver_oi, 
                                                                      grouping_tbl = grouping_tbl, 
                                                                      groups_oi = groups_oi, 
                                                                      target_regulation_df = target_regulation_df) %>% 
        .$group_plot + xlab("") + ylab("")
      ## aggregate
      sample_data = sample_data %>% 
        dplyr::filter(lr_interaction %in% order_lr_interaction) %>% 
        dplyr::mutate(lr_interaction = factor(lr_interaction))
      sample_data_agg <- lapply(unique(sample_data$group), function(gx){
        
        sample_datax <- subset(sample_data, group == gx)
        lapply(unique(sample_datax$id), function(idx){
          
          sample_dataxx <- subset(sample_datax, id == idx)
          data.frame(group = gx, 
                     lr_interaction = unique(sample_dataxx$lr_interaction),
                     scaled_LR_pb_prod = mean(sample_dataxx$scaled_LR_pb_prod),
                     prop_sender_receiver = sum(sample_dataxx$keep_sender_receiver == "Sender & Receiver present")/
                       nrow(sample_dataxx),
                     sender_receiver = unique(sample_dataxx$sender_receiver))
          
        }) %>% Reduce("rbind", .)
        
      }) %>% Reduce("rbind", .)
      
      sample_data_agg <- sample_data_agg[paste0(sample_data_agg$sender_receiver, 
                                                ":",
                                                sample_data_agg$lr_interaction) %in%
                                           paste0(lr_target_data$sender_receiver, 
                                                  ":",
                                                  lr_target_data$lr_interaction) ,]
      p1 = sample_data_agg %>% 
        ggplot(aes(group, 
                   lr_interaction, 
                   color = scaled_LR_pb_prod, 
                   size = prop_sender_receiver)) + 
        geom_point() + 
        facet_grid(sender_receiver ~ ., 
                   scales = "free", 
                   space = "free", 
                   switch = "y") + 
        scale_x_discrete(position = "top") + 
        theme_light() + 
        theme(axis.ticks = element_blank(), 
              axis.text.y = element_text(face = "bold.italic", size = 9), 
              axis.text.x = element_text(size = 9, angle = 90, hjust = 0), 
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 8),
              legend.key.size = unit(0.2, "inches"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.spacing.x = unit(0.4, "lines"), panel.spacing.y = unit(0.25, "lines"), 
              strip.text.x.top = element_text(size = 10, color = "black", angle = 0), 
              strip.text.y.left = element_text(size = 9, color = "black", angle = 0), 
              strip.background = element_rect(color = "darkgrey", fill = "whitesmoke", size = 1.5, linetype = "solid")) + 
        labs(color = "Scaled L-R\n pseudobulk product", 
             size = "Proportion of sufficient\n sender & receiver presence")
      max_lfc = abs(sample_data_agg$scaled_LR_pb_prod) %>% max()
      custom_scale_color = scale_color_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(), 
                                                 values = c(0, 0.35, 0.485, 0.5, 0.515, 0.65, 1), 
                                                 limits = c(-1 * max_lfc, max_lfc))
      p_lr_exprs = p1 + custom_scale_color + 
        xlab("Group") + ylab("Prioritized LR pairs")
      n_targets = length(order_targets)
      n_ligands = length(lr_target_data$id %>% unique())
      n_samples = grouping_tbl %>% 
        dplyr::filter(group %in% groups_oi) %>% 
        dplyr::pull(sample) %>% length()
      n_groups = grouping_tbl %>% 
        dplyr::pull(group) %>% unique() %>% length()
      n_groups_oi = grouping_tbl %>% 
        dplyr::filter(group %in% groups_oi) %>% 
        dplyr::pull(group) %>% unique() %>% length()
      legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_lr_exprs)), 
                                      ggpubr::as_ggplot(ggpubr::get_legend(p_lr_target)), 
                                      nrow = 1) %>% 
        patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_target_exprs)),
                              nrow = 2)
      if (is.null(heights)) {
        heights = c(n_ligands, n_groups_oi)
      }
      if (is.null(widths)) {
        widths = c(n_groups * 0.5 + 3, n_targets * 0.5 + 4)
      }
      if (plot_legend == FALSE) {
        design <- "AB\n               #C"
        combined_plot = patchwork::wrap_plots(A = p_lr_exprs + 
                                                theme(legend.position = "none", 
                                                      axis.ticks = element_blank()) + 
                                                theme(axis.title.x = element_text()), 
                                              B = p_lr_target + 
                                                theme(legend.position = "none", 
                                                      axis.ticks = element_blank()) + 
                                                ylab(""), 
                                              C = p_target_exprs + 
                                                theme(legend.position = "none"), 
                                              nrow = 2, 
                                              design = design, 
                                              widths = widths, 
                                              heights = heights)
        return(list(combined_plot = combined_plot, legends = legends))
      }
      else {
        design <- "AB\n               LC"
        combined_plot = patchwork::wrap_plots(A = p_lr_exprs + 
                                                theme(legend.position = "none", 
                                                      axis.ticks = element_blank()) + 
                                                theme(axis.title.x = element_text()), 
                                              B = p_lr_target + 
                                                theme(legend.position = "none", 
                                                      axis.ticks = element_blank()) + 
                                                ylab(""), 
                                              C = p_target_exprs + 
                                                theme(legend.position = "none"), 
                                              L = legends, 
                                              nrow = 2, 
                                              design = design, 
                                              widths = widths, 
                                              heights = heights)
        return(list(combined_plot = combined_plot, legends = legends))
      }
    }
    
  }

##### get top predictions across senders and receivers of interest #####
get.prior.tbl <- function(multinichenet_output = NULL,
                          n_top_lr = 50,
                          senders_oi = NULL,
                          receivers_oi = NULL){
  #
  prioritized_tbl_oi_all <- get_top_n_lr_pairs(
    prioritization_tables = multinichenet_output$prioritization_tables, 
    top_n = top_lr_circle,
    senders_oi = senders_oi, 
    receivers_oi = receivers_oi, 
    rank_per_group = F)
  prioritized_tbl_oi <- 
    multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
    filter(id %in% prioritized_tbl_oi_all$id) %>%
    distinct(id, 
             sender, 
             receiver, 
             ligand, 
             receptor, 
             group) %>% 
    left_join(prioritized_tbl_oi_all)
  prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] <- 0
  #
  return(prioritized_tbl_oi)
  
}


lr.target.filter <- function(multinichenet_output = NULL,
                             contrast_tbl = NULL,
                             cor_thresh = 0.5,
                             p_thresh = 0.05){
  
  if (nrow(multinichenet_output$lr_target_prior_cor) > 0) {
    
    lr_target_prior_cor_all <- multinichenet_output$lr_target_prior_cor %>%
      inner_join(
        multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% 
          distinct(ligand, 
                   target, 
                   direction_regulation, 
                   contrast)
      ) %>% 
      inner_join(contrast_tbl)
    
    lr_target_prior_cor_filtered_up <- lr_target_prior_cor_all %>% 
      filter(direction_regulation == "up") %>% 
      filter((pearson > cor_thresh | spearman > cor_thresh) &
               (pearson_pval < p_thresh | spearman_pval < p_thresh))
    lr_target_prior_cor_filtered_down <- lr_target_prior_cor_all %>% 
      filter(direction_regulation == "down") %>% 
      filter((pearson < -cor_thresh | spearman < -cor_thresh) &
               (pearson_pval < p_thresh | spearman_pval < p_thresh))
    lr_target_prior_cor_filtered <- bind_rows(
      lr_target_prior_cor_filtered_up, 
      lr_target_prior_cor_filtered_down)
  } else {
    lr_target_prior_cor_filtered <- data.frame(group = NA,
                                               sender = NA,
                                               receiver = NA)
  }
  
  return(lr_target_prior_cor_filtered)
}

##### mnn.lrt.plot #####
mnn.lrt.plot <- function(multinichenet_output = NULL,
                         contrast_tbl = NULL,
                         lr_target_cor_f = NULL,
                         lr_prod_activity = T,
                         ligand_target = T,
                         lr_target_cor = T,
                         n_top_lpa = 30,
                         receivers_oi = NULL,
                         senders_oi = NULL,
                         groups_oi = NULL){
  
  ## 0. get top LR pairs 
  prioritized_tbl_oi = get_top_n_lr_pairs(
    multinichenet_output$prioritization_tables, 
    top_n = n_top_lpa, 
    senders_oi = senders_oi, 
    receivers_oi = receivers_oi,
    groups_oi = groups_oi)
  
  #
  lr_prod_plot_list = NULL
  ligand_target_plot_list = NULL
  lr_target_cor_plot_list = NULL
  #
  if (nrow(prioritized_tbl_oi) > 0) {
    
    ## 1. lr_prod_plot
    if (lr_prod_activity) {
      
      ht_lr_prod_plot <- 10
      wd_lr_prod_plot <- 13
      lr_prod_plot <- make_sample_lr_prod_activity_plots_modify(
        multinichenet_output$prioritization_tables, 
        prioritized_tbl_oi)
      lr_prod_plot_list <- list("plt" = lr_prod_plot,
                                "ht" = ht_lr_prod_plot,
                                "wd" = wd_lr_prod_plot)
    }
    # 2. ligand_target
    if (ligand_target) {
      # get plot object
      ligand_target_plot_obj <- make_ligand_activity_target_plot_modify(
        group_oi = groups_oi, 
        receiver_oi = receivers_oi,
        prioritized_tbl_oi = prioritized_tbl_oi,
        prioritization_tables = multinichenet_output$prioritization_tables, 
        ligand_activities_targets_DEgenes = multinichenet_output$ligand_activities_targets_DEgenes, 
        contrast_tbl = contrast_tbl,
        grouping_tbl = multinichenet_output$grouping_tbl, 
        receiver_info = multinichenet_output$celltype_info, 
        ligand_target_matrix = ligand_target_matrix, 
        plot_legend = F)
      
      # format plot
      if (!is.null(ligand_target_plot_obj)) {
        
        ht_ligand_target <- ligand_target_plot_obj$dim_lt[["n_ligand"]] * 0.2 + 7
        wd_ligand_target <- ligand_target_plot_obj$dim_lt[["n_target"]] * 0.35 + 11
        ligand_target_plot <- (ligand_target_plot_obj[[1]] | 
                                 ligand_target_plot_obj[[2]]) + 
          patchwork::plot_layout(widths = c(wd_ligand_target - 3, 3))
        ligand_target_plot_list <- list("plt" = ligand_target_plot,
                                        "ht" = ht_ligand_target,
                                        "wd" = wd_ligand_target)
        
      }
    }
    # 3. lr_target_cor
    if (lr_target_cor) {
      #
      lr_target_cor_ff <- lr_target_cor_f %>%
        filter(group == groups_oi, 
               sender %in% senders_oi,
               receiver %in% receivers_oi)
      dir_sum <- (table(lr_target_cor_ff$target, 
                        lr_target_cor_ff$direction_regulation) > 0) %>% rowSums()
      bio_gene <- names(dir_sum)[dir_sum > 1]
      lr_target_cor_ff <- lr_target_cor_ff %>%
        filter(!target %in% bio_gene)
    }
    #
    if (nrow(lr_target_cor_ff) > 0) {
      # get plot object
      lr_target_cor_plot_obj <- 
        make_lr_target_correlation_plot_modify(
          prioritization_tables = multinichenet_output$prioritization_tables,
          prioritized_tbl_oi = prioritized_tbl_oi,
          lr_target_prior_cor_filtered = lr_target_cor_ff, 
          grouping_tbl = multinichenet_output$grouping_tbl,
          receiver_info = multinichenet_output$celltype_info,
          receiver_oi = receivers_oi,
          plot_legend = F)
      
      # format plot
      if (!is.null(lr_target_cor_plot_obj)) {
        
        ht_lr_target_cor <- sum(prioritized_tbl_oi$ligand %in% lr_target_cor_ff$ligand) * 0.27 + 6
        wt_lr_target_cor <- length(unique(lr_target_cor_ff$target)) * 0.2 + 14
        lr_target_cor_plot <- (lr_target_cor_plot_obj[[1]] | 
                                 lr_target_cor_plot_obj[[2]]) + 
          patchwork::plot_layout(widths = c(wt_lr_target_cor - 7, 7))
        lr_target_cor_plot_list <- list("plt" = lr_target_cor_plot,
                                        "ht" = ht_lr_target_cor,
                                        "wd" = wt_lr_target_cor)
      }
      
    }
  }
  
  return(list("lr_prod_plot_list" = lr_prod_plot_list,
              "ligand_target_plot_list" = ligand_target_plot_list,
              "lr_target_cor_plot_list" = lr_target_cor_plot_list))
}

##### mnn.lrt.plot.out #####
mnn.lrt.plot.out <- function(lr_prod_plot_list = NULL,
                             ligand_target_plot_list = NULL,
                             lr_target_cor_plot_list = NULL,
                             out_path = NULL,
                             plt_suffix = NULL){
  # 1. plot lr_prod
  if (!is.null(lr_prod_plot_list)) {
    
    pdf(glue("{out_path}lr_prod_activity_{plt_suffix}.pdf"),
        height = lr_prod_plot_list$ht, 
        width = lr_prod_plot_list$wd)
    print(lr_prod_plot_list$plt)
    dev.off()
    
  }
  # 2. plot ligand_target
  if (!is.null(ligand_target_plot_list)) {
    
    pdf(glue("{out_path}ligand_target_{plt_suffix}.pdf"),
        height = ligand_target_plot_list$ht, 
        width = ligand_target_plot_list$wd)
    print(ligand_target_plot_list$plt)
    dev.off()
    
  }
  # 4.1.3 plot lr_target_cor
  if (!is.null(lr_target_cor_plot_list)) {
    
    pdf(glue("{out_path}lr_target_cor_{plt_suffix}.pdf"),
        height = lr_target_cor_plot_list$ht, 
        width = lr_target_cor_plot_list$wd)
    print(lr_target_cor_plot_list$lr_prod_plot)
    dev.off()
  }
  
}
