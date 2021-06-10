
#' get_DE_genes
#'
#' Essentially a wrapper for scran::findMarkers. Returns a combined data.frame, each row corresponds
#' to cell type and gene which found to be DE in the cell type (using designated FDR as a threshold).
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param test.type String specifying which testing will be performed. Available options are the same as for scran::findMarkers, parameter test.use. Default test.type="binom".
#' @param pval.type String specifying how p-values are combined. Available options are the same as for scran::findMarkers, parameter pval.type Default pval.type="some".
#' @param FDR.thresh Scalar specifying threshold of FDR to be considered as significant. Default FDR=0.01.
#' @param ... Additional arguments
#'
#' @return
#' @export
#' @import scran
#' @import tibble
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 1000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$celltype = as.character(sample.int(5, n_col, replace = TRUE))
#' out = get_DE_genes(sce)
#'
get_DE_genes = function(sce , test.type = "binom", pval.type = "some", FDR.thresh = 0.01, ...){
  args = c(as.list(environment()) , list(...))
  if (!"check_args" %in% names(args)){
    sce = .prepare_sce(sce)
    out = .general_check_arguments(args) & .check_celltype_in_sce(sce)
  }
  else {
    if (args[[which(names(args) == "check_args")]]){
      sce = .prepare_sce(sce)
      out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_celltype_in_sce(sce)
    }
  }

  sce$celltype = as.character(sce$celltype)
  # get potential relevant genes
  markers = findMarkers(sce , groups=sce$celltype, direction = "up", pval.type=pval.type, test = test.type, assay.type = "logcounts")
  markers = lapply(1:length(markers) , function(i) {
    current.markers = as.data.frame(markers[[i]])
    current.markers = rownames_to_column(current.markers , var = "gene")
    current.markers$celltype = names(markers)[i]
    current.markers = current.markers[, c("celltype" , "gene" , "FDR", "summary.logFC")]
    return(current.markers)
  })
  markers = do.call(rbind , markers)
  if (!is.null(FDR.thresh)){
    markers = markers[!is.na(markers$FDR) & markers$FDR <= FDR.thresh , ]
  }
  if (nrow(markers) < 1){
    message("No DE genes discovered with these settings - consider tuning test, type and/or FDR threshold.")
    return(NULL)
  }
  else {
    return(markers[, c("celltype" , "gene", "summary.logFC")])
  }
}

