#' get_neighborhood_preservation_scores_alternative
#'
#' Calculates cell neighborhood preservation scores by comparing neighbors from True and Selection k-NN graphs.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
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
#' out = get_neighborhood_preservation_scores_alternative(sce, genes.selection = genes.selection)
#'
get_neighborhood_preservation_scores_alternative = function(sce, neighs.all_stat = NULL, genes.all = rownames(sce),
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
    score = .get_neighborhood_preservation_scores_single_batch_alternative(sce , neighs.all_stat = neighs.all_stat, genes.all = genes.all,
                                                               genes.selection = genes.selection, n.neigh = n.neigh ,
                                                               nPC.all = nPC.all , nPC.selection = nPC.selection)
    return(score)
  }
  else {
    if (is.null(neighs.all_stat)){
      neighs.all_stat = get_neighs.all_stat(sce , genes.all = genes.all , batch = batch, n.neigh = n.neigh , nPC.all = nPC.all)
    }
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    score = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.sce = sce[, idx]
      current.neighs.all_stat = neighs.all_stat[[which(names(neighs.all_stat) == current.batch)]]
      current.score = .get_neighborhood_preservation_scores_single_batch_alternative(current.sce , neighs.all_stat = current.neighs.all_stat,
                                                               genes.all = genes.all ,
                                                               genes.selection = genes.selection, n.neigh = n.neigh ,
                                                               nPC.all = nPC.all, nPC.selection = nPC.selection)
      return(current.score)
    })
    score = do.call(rbind, score)
    return(score)
  }
}


#' @import Rfast
#' @importFrom irlba prcomp_irlba
.get_neighborhood_preservation_scores_single_batch_alternative = function(sce, neighs.all_stat = NULL, genes.all = rownames(sce),
                                                              genes.selection, n.neigh = 5, nPC.all = 50, nPC.selection = NULL){
  set.seed(32)
  sce = sce[genes.all , ]

  if (is.null(neighs.all_stat)){
    neighs.all_stat = get_neighs.all_stat(sce , genes.all = genes.all , batch = NULL, n.neigh = n.neigh , nPC.all = nPC.all)
  }
  counts = neighs.all_stat$counts
  neighs.all = neighs.all_stat$neighs.all
  neighs.compare = .get_mapping(sce , genes = genes.selection, batch = NULL, n.neigh = n.neigh, nPC = nPC.selection)
  neighs.compare = neighs.compare$cells_mapped
  cells_random = sample(colnames(sce), min(100, ncol(sce)))

  counts = counts[order(rownames(counts)),]
  neighs.all = neighs.all[order(rownames(neighs.all)),]
  neighs.compare = neighs.compare[order(rownames(neighs.compare)),]

  score = lapply(rownames(counts) , function(cell){
    dist_neighs.all_coords = median(dista(t(counts[cell, ]), counts[neighs.all[cell,], ]))
    dist_neighs.compare_coords = median(dista(t(counts[cell, ]), counts[neighs.compare[cell,], ]))
    dist_random_coords = median(dista(t(counts[cell, ]), counts[cells_random, ]))
    current_score = ( dist_random_coords - dist_neighs.compare_coords )/( dist_random_coords - dist_neighs.all_coords )
    out = data.frame(cell_score = current_score )
  })
  score = do.call(rbind , score)
  score$cell = rownames(neighs.compare)
  score = score[, c("cell", "cell_score")]
  return(score)
}


#' get_neighs.all_stat
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
get_neighs.all_stat = function(sce , genes.all = rownames(sce) , batch = NULL, n.neigh = 5 , nPC.all = 50){
  if (is.null(batch)){
    out = .get_neighs.all_stat_single_batch(sce , genes.all = genes.all , n.neigh = n.neigh , nPC.all = nPC.all)
    return(out)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    neighs.all = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.sce = sce[, idx]
      out =  .get_neighs.all_stat_single_batch(current.sce , genes.all = genes.all , n.neigh = n.neigh , nPC.all = nPC.all)
    })
    names(neighs.all) = unique(batchFactor)
    return(neighs.all)
  }
}

.get_neighs.all_stat_single_batch = function(sce , genes.all = rownames(sce) , n.neigh = 5 , nPC.all = 50){

  counts = t(as.matrix(logcounts(sce)))
  pcs = suppressWarnings( prcomp_irlba(counts , n = min(nPC.all, (nrow(counts)-1) , (ncol(counts) - 1))) )
  counts = pcs$x
  rownames(counts) = colnames(sce)

  neighs.all = .get_mapping(sce , genes = genes.all, batch = NULL, n.neigh = n.neigh, nPC = nPC.all , get.dist = F)
  neighs.all = neighs.all$cells_mapped

  out = list(counts = counts, neighs.all = neighs.all)
  return(out)
}


