
#' get_neighborhood_preservation_scores
#'
#' Calculates cell neighborhood preservation scores by comparing neighbors from True and Selection k-NN graphs.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param neighs.all If not NULL, contains information about True kNN-graph (for each cell - ordered neighbors and distances). Useful to have a priori if cell score will be recalculated multiple times.
#' @param genes.all String specifying genes to be used for construction of True kNN-graph.
#' @param genes.selection String specifying genes to be used for construction of Selection kNN-graph.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Positive integer > 1, specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC.all Scalar specifying number of PCs to use for construction of True kNN-graph. Default nPC.all=50.
#' @param nPC.selection Scalar specifying number of PCs to use for construction of True kNN-graph. Default nPC.selection=NULL (no PCA).
#' @param ... Additional arguments
#'
#' @return data.frame, each row corresponds to cell from counts matrix, contains field cell_score = cell neighborhood preservation score
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
#' genes.selection = sample(rownames(sce) , 20)
#' out = get_neighborhood_preservation_scores(sce, genes.selection = genes.selection)
#'
get_neighborhood_preservation_scores = function(sce, neighs.all = NULL,  genes.all = rownames(sce),
                                  genes.selection, batch = NULL, n.neigh = 5, nPC.all = 50, nPC.selection = NULL, ...){
  args = c(as.list(environment()) , list(...))
  if (!"check_args" %in% names(args)){
    sce = .prepare_sce(sce)
    out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.all) & .check_genes_in_sce(sce , genes.selection)
  }
  else {
    if (args[[which(names(args) == "check_args")]]){
      sce = .prepare_sce(sce)
      out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.all) & .check_genes_in_sce(sce , genes.selection)
    }
  }

  if (is.null(batch)){
    out = .get_neighborhood_preservation_scores_single_batch(sce , neighs.all = neighs.all, genes.all = genes.all,
                                                             genes.selection = genes.selection, n.neigh = n.neigh , nPC.all = nPC.all , nPC.selection = nPC.selection)
    return(out)
  }
  else {
    if (is.null(neighs.all)){
      neighs.all = get_z_scaled_distances(sce , genes.all = genes.all , batch = batch, n.neigh = n.neigh, nPC.all = nPC.all, check_args = FALSE)
    }
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    score = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.sce = sce[, idx]
      current.neighs.all = neighs.all[[which(names(neighs.all) == current.batch)]]
      out = .get_neighborhood_preservation_scores_single_batch(current.sce , neighs.all = current.neighs.all , genes.all = genes.all ,
                                                               genes.selection = genes.selection, n.neigh = n.neigh , nPC.all = nPC.all, nPC.selection = nPC.selection)
      return(out)
    })
    score = do.call(rbind, score)
    return(score)
  }
}


#'
.get_neighborhood_preservation_scores_single_batch = function(sce, neighs.all = NULL,  genes.all = rownames(sce),
                                                              genes.selection, n.neigh = 5, nPC.all = 50, nPC.selection = NULL){
  if (is.null(neighs.all)){
    neighs.all = get_z_scaled_distances(sce, genes.all = genes.all, batch = NULL, n.neigh = n.neigh, nPC.all = nPC.all)
  }
  else {
    if (!class(neighs.all) == "list" | !sum(names(neighs.all) == c("cells_mapped" , "distances")) == 2){
      stop("neighs.all is of the wrong format - should be a list, with elements names 'cells_mapped' and 'distances'.")
      return(F)
    }
  }
  neighs.compare = .get_mapping(sce , genes = genes.selection, batch = NULL, n.neigh = n.neigh, nPC = nPC.selection)
  neighs.compare = neighs.compare$cells_mapped
  neighs.all.cells_mapped = neighs.all$cells_mapped
  neighs.all.distances = neighs.all$distances
  n.cells = ncol(sce)

  score = lapply(1:nrow(neighs.compare) , function(i){
    cells = neighs.all.cells_mapped[i,]
    idx_all = c(1:n.neigh)
    idx_compare = which(cells %in% neighs.compare[i,] )

    dist_all = neighs.all.distances[i, idx_all]
    dist_compare = neighs.all.distances[i, idx_compare]
    current.score = median(-dist_compare)/median(-dist_all)
    out = data.frame(cell_score = current.score )
    return(out)
  })
  score = do.call(rbind , score)
  score$cell = rownames(neighs.compare)
  return(score)
}




#' get_z_scaled_distances
#'
#' For each cell: ranks all other cells based on the transcriptional similarity and returns ordered z-scaled distances (lower distance -- closer cell).
#' It is an intermediate step for get_neighborhood_preservation_scores. It can be handy to calculate this part separately to be recycled multiple times
#' to compare different seelctions.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix.
#' @param genes.all String specifying genes to be used for construction of True kNN-graph.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Positive integer > 1, specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC.all Scalar specifying number of PCs to use for construction of True kNN-graph. Default nPC.all=50.
#' @param ... Additional arguments.
#'
#' @return kNN-graph with assigned neighbors and corresponding z-scored distances.
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
#' out = get_z_scaled_distances(sce)
get_z_scaled_distances = function(sce , genes.all = rownames(sce) , batch = NULL, n.neigh = 5 , nPC.all = 50, ...){
  args = c(as.list(environment()) , list(...))
  if (!"check_args" %in% names(args)){
    sce = .prepare_sce(sce)
    out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.all)
  }
  else {
    if (args[[which(names(args) == "check_args")]]){
      sce = .prepare_sce(sce)
      out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.all)
    }
  }
  if (is.null(batch)){
    out = .get_z_scaled_distances_single_batch(sce , genes.all = genes.all , n.neigh = n.neigh , nPC.all = nPC.all)
    return(out)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    neighs.all = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.sce = sce[, idx]
      out =  .get_z_scaled_distances_single_batch(current.sce , genes.all = genes.all , n.neigh = n.neigh , nPC.all = nPC.all)
    })
    names(neighs.all) = unique(batchFactor)
    return(neighs.all)
  }
}


.get_z_scaled_distances_single_batch = function(sce , genes.all = rownames(sce) , n.neigh = 5 , nPC.all = 50){
  neighs.all = .get_mapping(sce , genes = genes.all, batch = NULL, n.neigh = "all", nPC = nPC.all , get.dist = T)
  distances = neighs.all$distances
  distances_scaled = t( apply(distances , 1 , function(x) scale(x)) )
  # if there are NaNs (i.e. same distance to all neighbors, return 0) -- the final dist will be NaN which is ok for this cases
  distances_scaled[is.na(distances_scaled)] = 0

  rownames(distances_scaled) = rownames(distances)
  neighs.all$distances = distances_scaled
  return(neighs.all)
}
