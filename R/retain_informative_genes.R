
#' retain_informative_genes
#'
#' Function selects variable genes and retains them to represent counts matrix for downstream analysis. This is essentially a wrapper for scran::getTopHVGs +
#' simple check for mt-genes.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param n Scalar (or NULL, if not applied) specifying number of top HVGs to be return. Default=NULL.
#' @param var.thresh Scalar (float) specifying threshold for variation to filter HVGs. Default=0.
#' @param select.hvgs Boolean specifying if we want to discard not HVGs. Default=TRUE.
#' @param discard.mt Boolean specifying if mitochondrial genes should be discarded. Default=TRUE.
#'
#' @return Reduced sce with prefiltered set of genes.
#' @export
#' @import SingleCellExperiment
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 1000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' out = retain_informative_genes(sce)
#'
retain_informative_genes = function(sce, n = NULL, var.thresh = 0, select.hvgs = TRUE, discard.mt = TRUE){
  args = c(as.list(environment()))
  sce = .prepare_sce(sce)
  if (.general_check_arguments(args)){
    if (select.hvgs){
      hvgs = .get_hvgs(sce, n = n, var.thresh = var.thresh)
      sce = sce[hvgs,]
    }
    if (nrow(sce) < 2){
      stop("Less than 2 genes are selected. Consider checking your counts data, increasing n and/or decreasing var.thresh")
    }
    else {
      if (discard.mt){
        rownames.sce = rownames(sce)
        idx = which(!grepl("mt-" , rownames.sce) & !grepl("MT-" , rownames.sce) & !grepl("Mt-" , rownames.sce))
        sce = sce[idx, ]
      }
    }
    cat(paste0(nrow(sce) , " genes retained"))
    if (nrow(sce) < 2){
      stop("Less than 2 genes are selected. Consider checking your counts data, increasing n and/or decreasing var.thresh")
    }
    else {
      return(sce)
    }
  }
}


#' @import scran
.get_hvgs = function(sce, n = NULL, var.thresh = 0){
  out = tryCatch(
    {
      dec.sce = modelGeneVar(sce)
      hvg.genes = getTopHVGs(dec.sce, var.threshold = var.thresh)
      if (!is.null(n)){
        if (n < length(hvg.genes)){
          hvg.genes = getTopHVGs(dec.sce, n = n)
        }
      }
      return(hvg.genes)
    },
    error = function(dump){
      stop("Can not perform modelGeneVar on this counts matrix - check your input data.")
    }
  )
  return(out)
}

