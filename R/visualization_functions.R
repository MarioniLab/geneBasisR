

#' plot_coexpression
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector containing gene names to evaluate for co-expression.
#' @param title String to be passed as a title. Default title=NULL.
#' @param ... Additional arguments to pass.
#' @return Heatmap for co-expression.
#' @export
#' @import ggcorrplot
#'
plot_coexpression = function(sce , genes , title = NULL, ...){
  current.sce = sce[genes , ]
  current.counts = as.matrix(logcounts(current.sce))
  corr.stat = as.data.frame( cor(t(current.counts), method = "pearson") )
  p <- ggcorrplot(corr.stat, hc.order = TRUE, outline.col = "white", ggtheme = ggplot2::theme_gray, colors = c("#6D9EC1", "white", "#E46726")) +
    theme(axis.text.x = element_text(size=6), axis.text.y = element_text(size=6)) +
    ggtitle(title)
  p
  return(p)
}



#' plot_mapping_heatmap
#'
#' @param mapping data.frame containing celltype mapping information.
#' @param levels Character vector specifying the order for plotting.
#' @param title String to be passed as a title. Default title=NULL.
#' @param ... Additional arguments to pass.
#' @return Heatmap for confusion matrix.
#' @export
#' @import ggpubr ggplot2 viridis
#'
plot_mapping_heatmap = function(mapping , levels = unique(mapping$celltype) , title = NULL, ...){
  mapping$celltype = as.character(mapping$celltype)
  mapping$celltype_mapped = as.character(mapping$mapped_celltype)
  tab = table(mapping$celltype , mapping$celltype_mapped)
  tab = sweep(tab, 1, rowSums(tab), "/")
  tab = as.data.frame( tab )
  colnames(tab) = c("celltype", "celltype_mapped", "n")
  tab$celltype = factor(tab$celltype , levels = levels)
  tab$celltype_mapped = factor(tab$celltype_mapped , levels = levels)
  tab = tab[!is.na(tab$celltype) , ]
  p <- ggplot(tab, aes(x = celltype , y = celltype_mapped, fill = n)) +
    geom_tile() + viridis::scale_fill_viridis(discrete = F) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle(title)
  return(p)
}


#' plot_umaps_w_counts
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector containing gene names for which to plot UMAPs.
#' @param size Size of dots for geom_point
#' @param ncol Positive integer specifying number of columns for ggarrange
#' @return For each gene: scater plot with UMAP-coordinates colored by gene expression
#' @export
#' @import tibble ggpubr ggplot2 viridis SingleCellExperiment
#'
plot_umaps_w_counts = function(sce , genes, size = .25, ncol = NULL){
  # SCE should contain reducedDim = UMAP, which contains 2 columns: x and y
  umaps = as.data.frame(reducedDim(sce , "UMAP"))
  umaps = rownames_to_column(umaps, var = "cell")
  plots = lapply(genes, function(gene){
    counts = data.frame(cell = colnames(sce) , counts = as.numeric(logcounts(sce)[gene , ]))
    current.umaps = merge(umaps , counts)
    current.umaps = current.umaps[order(current.umaps$counts) , ]
    p <- ggplot(current.umaps , aes(x = x , y = y , col = counts)) +
      geom_point(size=size) +
      scale_color_gradient(low = "azure3" , high = "darkgreen") +
      theme_classic() +
      theme(legend.position="none") +
      ggtitle(gene) +
      labs(x = "UMAP-1" , y = "UMAP-2")
    return(p)
  })
  if (is.null(ncol)){
    p = ggarrange(plotlist = plots)
  }
  else {
    p = ggarrange(plotlist = plots, ncol = ncol)
  }
  return(p)
}


#' plot_expression_heatmap
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector containing gene names.
#' @param value.type String specifying whether to plot average expression (= "mean") or fraction of cells with non-zero counts(= "frac").
#'
#' @return Heatmap with average across (per gene/per celltype)
#' @export
#' @import ggpubr ggplot2 viridis SingleCellExperiment
#'
plot_expression_heatmap = function(sce , genes , value.type){
  sce = sce[genes, ]
  stat = lapply(unique(sce$celltype) , function(celltype){
    current.sce = sce[, sce$celltype == celltype]
    current.counts = as.matrix( assay(current.sce, "logcounts" ))
    if (value.type == "mean"){
      current.stat = data.frame(gene = rownames(sce) , value = apply(current.counts , 1 , mean))
    }
    else if (value.type == "frac"){
      current.stat = data.frame(gene = rownames(sce) , value = apply(current.counts , 1 , function(x) mean(x > 0)))
    }
    current.stat$celltype = celltype
    return(current.stat)
  })
  stat = do.call(rbind , stat)
  stat$gene = factor(stat$gene , levels = genes)
  p <- ggplot(data=stat , aes(x = celltype , y = gene , fill = value)) +
    geom_tile() +
    scale_fill_viridis(discrete = F) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.position = "none") +
    coord_flip()
  return(p)
}




#' plot_redundancy_stat
#'
#' @param redundancy_stat data.frame - an output of get_redundancy_stat
#' @param celltypes Character vector specifying for which celltypes the output should be plotted. Default celltypes = unique(redundancy_stat$celltype)
#' @param genes Character vector specifying for which genes the output should be plotted. Default genes = unique(redundancy_stat$gene)
#'
#' @return Heatmap of redundancy in celltype mapping
#' @export
#' @import wesanderson ggplot2 ggpubr

plot_redundancy_stat = function(redundancy_stat, celltypes = unique(redundancy_stat$celltype) , genes = unique(redundancy_stat$gene)){
  redundancy_stat = redundancy_stat[redundancy_stat$frac_correctly_mapped_all > 0 &
                                      redundancy_stat$celltype %in% celltypes &
                                      redundancy_stat$gene %in% genes,]
  pals = wes_palette("Darjeeling1")
  p = ggplot(redundancy_stat, aes(x = gene , y = celltype , fill = frac_correctly_mapped_ratio)) +
    geom_tile() +
    scale_fill_gradient2(low = pals[1] , high = pals[2] , midpoint = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = "Gene" , y = "Celltype")
  return(p)

}
