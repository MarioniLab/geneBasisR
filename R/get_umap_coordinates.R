
#' get_umap_coordinates
#'
#' Function returns UMAP-coordinates if onle selected genes are used.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector specifying genes to use for calculating umap-coordinates.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param nPC Scalar specifying number of PCs to use. Default = length(genes)-1
#'
#' @return data.frame, where UMAP-coordinates are assigned as x and y.
#' @export
#' @import SingleCellExperiment
#' @import batchelor
#' @importFrom uwot umap
#' @importFrom irlba prcomp_irlba
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 1000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' out = get_umap_coordinates(sce, genes = rownames(sce)[1:10])
#'
get_umap_coordinates = function(sce, genes, batch = NULL, nPC = length(genes) - 1){
  set.seed(32)
  sce = sce[genes, ]
  counts = as.matrix(logcounts(sce))
  meta = as.data.frame(colData(sce))
  if (is.null(batch)){
    pcs = suppressWarnings( prcomp_irlba(t(counts) , n = min(nPC, (nrow(counts)-1) , (ncol(counts) - 1))) )
    counts = pcs$x
    rownames(counts) = colnames(sce)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    counts = multiBatchPCA(counts , batch = batchFactor , d = nPC)
    counts = do.call(reducedMNN , as.list(counts))
    counts = counts$corrected
  }
  # get umap-coordinates
  umaps = as.data.frame( umap(counts, min_dist = 1) )
  rownames(umaps) = rownames(counts)
  colnames(umaps) = c("x" , "y")
  umaps = rownames_to_column(umaps, var = "cell")
  return(umaps)
}

