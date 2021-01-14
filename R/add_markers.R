#' Returns markers between selected celltypes
#'
#' @param sce \code{SingleCellExperiment} object representing scRNA-seq counts matrix containing celltype field.
#' @param celltypes Celltypes for which we want to find markers. Default=NULL and in this case all celltypes will be included.
#' @param Top Positive integer representing how many Top markers should be returned.
#'
#' @return A \code{data.frame} object with selected markers.
#' @export
#' @importFrom scran findMarkers
#' @importFrom tibble rownames_to_column
#' @examples
#' require(SingleCellExperiment)
#' n_row = 30000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$celltype = as.factor(sample(1:10, n_col, replace=TRUE))
#' out = add_markers(sce)
add_markers = function(sce, celltypes_distinguish = NULL, celltypes_add_markers = celltypes_distinguish, test = "t", Top=1){
  if (!.check_counts_matrix_correct(sce)) {
    stop()
  } else {
    if (!is.null(celltypes_distinguish)){
      sce = sce[, sce$celltype %in% celltypes_distinguish]
    }
    markers = findMarkers(sce , groups=sce$celltype, direction = "up", test = test, assay.type = "logcounts", pval.type="any")
    markers.df = lapply(celltypes_add_markers, function(celltype){
      current_markers = as.data.frame( markers[[which(names(markers) == celltype)]] )
      current_markers$celltype = celltype
      current_markers = rownames_to_column(current_markers , var = "gene")
      return(current_markers[current_markers$Top <= Top , c("Top", "p.value", "FDR", "summary.logFC", "celltype" , "gene")])
    })
    markers.df = do.call(rbind , markers.df)
    return(markers.df)
  }
}
