
#' retain_informative_genes
#'
#' Function selects variable genes and retains them to represent counts matrix for downstream analysis. This is essentially a wrapper for scran::getTopHVGs +
#' simple check for mt-genes.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param n Scalar (or NULL, if not applied) specifying number of top HVGs to be return. Default=NULL.
#' @param var.thresh Scalar (float) specifying threshold for variation to filter HVGs. Default=0.
#' @param discard.mt Boolean specifying if mitochondrial genes should be discarded. Default=TRUE.
#'
#' @return Reduced sce with prefiltered set of genes.
#' @export
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 1000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' out = retain_informative_genes(sce, n = 100)
#'
retain_informative_genes = function(sce, n = NULL, var.thresh = 0, discard.mt = T){
  sce = .prepare_sce(sce)
  if (.check_sce(sce)){
    hvgs = .get_hvgs(sce, n = n, var.thresh = var.thresh)
    sce = sce[hvgs,]
    if (discard.mt){
      rownames.sce = rownames(sce)
      idx = which(!grepl("mt-" , rownames.sce) & !grepl("MT-" , rownames.sce) & !grepl("Mt-" , rownames.sce))
      sce = sce[idx, ]
    }
    cat(paste0(nrow(sce) , " genes retained"))
    return(sce)
  }
}


#' @import scran
.get_hvgs = function(sce, n = NULL, var.thresh = 0){
  dec.sce = modelGeneVar(sce)
  hvg.genes = getTopHVGs(dec.sce, var.threshold = var.thresh)
  if (!is.null(n)){
    if (n < length(hvg.genes)){
      hvg.genes = getTopHVGs(dec.sce, n = n)
    }
  }
  return(hvg.genes)
}

