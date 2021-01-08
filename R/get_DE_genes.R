
#' Returns all genes that are deferentially expressed in at least one celltype
#'
#' @param sce SingleCellExperiment object representing scRNA-seq counts matrix containing celltype field.
#' @param test A parameter of findMarkers function determining which statistical test should be used for selection of DE genes. Default=binom.
#' @param FDR.thresh A False Discovery Rate threshold determining whether given gene would be considered as DE. Default=0.01
#' @param filter_by_ct_specificity A boolean parameter to whether we should perform initial filtering by gene expression being CT specific (q75 == 0 in more than half of celltypes)
#'
#' @return data.frame object where each row represents gene/celltype if the gene is assigned as DE in a given celltype
#' @export
#' @importFrom scran findMarkers
#' @import SingleCellExperiment
#' @examples
#' require(SingleCellExperiment)
#' n_row = 30000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$celltype = as.factor(sample(1:5, n_col, replace=TRUE))
#'
#' markers.binom = get_DE_genes(sce)
#' markers.t = get_DE_genes(sce , test = "t")
#'
get_DE_genes = function(sce , test = "binom", FDR.thresh = 0.01 , filter_by_ct_specificity = TRUE){

  # requires there is a column named celltype and assay named logcounts: abort if the criteria are not met
  if (!.check_counts_matrix_correct(sce)) {
    stop()
  } else{
    if (filter_by_ct_specificity){
      sce = .filter_by_ct_specificity(sce)
    }
    # get relevant genes
    markers <- findMarkers(sce , groups=sce$celltype, direction = "up", pval.type="some", test = test, assay.type = "logcounts")
    # put together in data.frame
    celltypes = names(markers)
    markers = lapply(1:length(celltypes), function(i){
      current.markers = as.data.frame(markers[[i]])
      current.markers = current.markers[!is.na(current.markers$FDR) & current.markers$FDR < FDR.thresh , ]
      if (nrow(current.markers) > 0){
        out = data.frame( celltype = celltypes[i], gene = rownames(current.markers))
        return(out)
      }
    })
    markers = do.call(rbind,markers)
    markers$gene = factor(markers$gene)
    markers$celltype = factor(markers$celltype)
    out = list(genes_de.stat = markers , genes_de = levels(markers$gene))
    return(out)
  }
}

#' @importFrom dplyr summarise group_by
#' @importFrom stats quantile
.filter_by_ct_specificity = function(sce){
  if (!.check_counts_matrix_correct(sce)) {
    stop()
  } else{
    stat = lapply(unique(sce$celltype) , function(celltype){
      logcounts = as.matrix( logcounts(sce[, sce$celltype == celltype]) )
      current.stat = apply(logcounts , 1 , function(x) quantile(x , .75))
      current.stat = data.frame(gene = names(current.stat) , q75 = current.stat , celltype = celltype)
    })
    stat = do.call(rbind , stat)
    stat = as.data.frame( stat %>% group_by(gene) %>% summarise(frac_q75_positive = mean(q75 > 0)))

    # keep only lowly expressed and specific celltypes
    thresh.n_ct = 0.5
    genes = as.character( stat$gene[stat$frac_q75_positive < thresh.n_ct] )
    sce = sce[rownames(sce) %in% genes , ]
    return(sce)
  }
}
