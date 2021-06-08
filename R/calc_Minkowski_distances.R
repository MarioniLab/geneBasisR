
#' calc_Minkowski_distances
#'
#' Function returns Minkowski distances for kNN-graphs constructed from counts matrix in sce (stored in assay logcounts) using specified genes.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector containing genes that are used to construct kNN-graph.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Scalar specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC Scalar (or NULL) specifying number of PCs to use for kNN-graph. Default nPC=NULL (no PCA).
#' @param genes.predict Character vector containing genes for which we want to calculate Minkowsky distances. Default genes.predict = rownames(sce).
#' @param p.minkowski Order of Minkowski distance. Default p.minkowski=3.
#' @param ... Additional arguments
#'
#' @return data.frame, field 'gene' to gene from genes.predict; field 'dist' corresponds to calculated Minkowski distance.
#' @export
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 1000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.character(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' genes = rownames(sce)
#' out = calc_Minkowski_distances(sce, genes)
#'
calc_Minkowski_distances = function(sce , genes , batch = NULL , n.neigh = 5 , nPC = NULL , genes.predict = rownames(sce) , p.minkowski = 3, ...){
  args = c(as.list(environment()) , list(...))
  if (!"check_args" %in% names(args)){
    sce = .prepare_sce(sce)
    out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.predict) & .check_genes_in_sce(sce , genes)
  }
  else {
    if (args[[which(names(args) == "check_args")]]){
      sce = .prepare_sce(sce)
      out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.predict) & .check_genes_in_sce(sce , genes)
    }
  }
  if (!is.null(genes)){
    neighs = .get_mapping(sce , genes = genes, batch = batch , n.neigh = n.neigh , nPC = nPC)
    neighs = neighs$cells_mapped
  }
  else {
    neighs = .initiate_random_mapping(sce , batch = batch , n.neigh = n.neigh)
  }
  counts_predict = as.matrix(logcounts(sce[genes.predict , ]))
  stat_predict = lapply(1:ncol(neighs) , function(j){
    cells = neighs[,j]
    current.stat_predict = counts_predict[, cells]
    return(current.stat_predict)
  })
  stat_predict = Reduce("+", stat_predict) / length(stat_predict)
  stat_real = counts_predict[, rownames(neighs)]
  stat = lapply(1:nrow(counts_predict) , function(i){
    out = data.frame(gene = rownames(counts_predict)[i] , dist = as.numeric(dist(rbind(stat_real[i,] , stat_predict[i,]) , method = "minkowski" , p = p.minkowski)))
    return(out)
  })
  stat = do.call(rbind , stat)
  return(stat)
}