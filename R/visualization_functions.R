## This script contains various plotting functions to assist in the analysis of the results


#' Distribution of gene expression per celltype
#'
#' @param sce \code{SingleCellExperiment} object representing scRNA-seq counts matrix containing celltype field.
#' @param gene Gene which expression we want to plot.
#' @param assay A string specifying which assay from sce to use. Default=logcounts.
#' @param title A string specifying title. Default=gene.
#'
#' @return A \code{ggplot} object containing boxplots per celltype.
#' @export
#' @importFrom ggplot2 ggplot theme_classic theme ggtitle
#' @examples
plot_expr_distribution = function(sce , gene , assay = "logcounts" , title = gene){
  if (!gene %in% rownames(sce)){
    stop("Can not find gene in the counts matrix. Ensure that given entry exists.")
  }
  if (!assay %in% c("counts" , "logcounts")){
    stop("Option 'assay' have to be either 'counts' or 'logcounts'.")
  }
  if (!assay %in% names(assays(sce))){
    stop("Chosen assay option does not exist in counts matrix.")
  }
  counts = data.frame(cell = sce$cell ,
                      celltype = sce$celltype ,
                      counts = as.numeric( assay(sce[gene, ], assay)) )
  p <- ggplot(data=counts , aes(x = celltype , y = counts , fill = celltype)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.position = "none") +
    ggtitle(title)
  return(p)
}


#' Heatmap representing how well celltypes are mapped
#'
#' @param mapping A \code{data.frame} object representing mapping.
#' @param A string specifying title. Default=gene.
#'
#' @return A \code{ggplot} object containing boxplots per celltype. X-axis corresponds to actual celltype, Y-axis corresponds to mapped celltypes, and color gradient corresponds to fraction of cells from actual celltype that are mapped to mapped celltype.
#' @export
#' @importFrom ggplot2 ggplot theme_classic theme ggtitle
#' @import viridis

plot_mapping_heatmap = function(mapping , title = NULL){
  if (!.valid_mapping_df(mapping)) {
    stop()
  } else {
    if (!is.null(title) & !is(title, "string")){
      stop("Option 'title' should be either NULL or a string.")
    } else {
      mapping$celltype = as.character(mapping$celltype)
      mapping$celltype_mapped = as.character(mapping$celltype_mapped)
      tab = table(mapping$celltype , mapping$celltype_mapped)
      tab = sweep(tab, 1, rowSums(tab), "/")
      tab = as.data.frame( tab )
      colnames(tab) = c("celltype", "celltype_mapped", "n")
      tab$celltype = factor(tab$celltype , levels = unique(sce$celltype))
      tab$celltype_mapped = factor(tab$celltype_mapped , levels = c(unique(sce$celltype),"Unmapped"))
      tab = tab[!is.na(tab$celltype) , ]
      p <- ggplot(tab, aes(x = celltype , y = celltype_mapped, fill = n)) +
        geom_tile() + viridis::scale_fill_viridis(discrete = F) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggtitle(title)
      return(p)
    }
  }
}
