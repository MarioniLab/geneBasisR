

#' Calculates gene prediction scores by comparing True and Selection k-NN graphs.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes.selection Character vector containing names of selected genes.
#' @param genes.all Character vector containing names of all (potentially relevant/variable) genes.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Scalar specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC.all Scalar (or NULL) specifying number of PCs to use for construction of True kNN-graph. Default nPC.all=50.
#' @param nPC.selection Scalar (or NULL) specifying number of PCs to use for construction of Selection kNN-graph. Default nPC.selection=NULL (no PCA).
#' @param genes.predict Character vector containing names of genes for which we want to calculate gene prediction score.
#' @param method Character specifying method for correlation. Availbale options are c("spearman" , "pearson"). Default method="spearman".
#' @param corr_all.thresh Scalar specifying suitable threshold for correlation to consider (on True graph).
#' @param stat_all If correlation-stat is pre-calculated for True graph, pass it here (default stat_all=NULL). This is useful if this variable will be used re-used multiple times.
#'
#' @return data.frame, each row corresponds to gene, contains field gene_score = gene prediction score.
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
#' out = get_gene_prediction_scores(sce, genes.selection)
#'
get_gene_prediction_scores = function(sce, genes.selection, genes.all = rownames(sce) , batch = NULL, n.neigh = 5, nPC.all = 50 , nPC.selection = NULL,
                                      genes.predict = genes.all, method = "spearman" , corr_all.thresh = 0.25 , stat_all = NULL){
  sce = .prepare_sce(sce)
  if (is.null(stat_all)){
    stat_all = suppressWarnings( .get_gene_corr(sce, genes.all, batch = batch, n.neigh = n.neigh,
                                                nPC = nPC.all, genes.predict = genes.predict, method = method) )
    colnames(stat_all) = c("gene" , "corr_all")
  }
  stat_all = stat_all[stat_all$corr_all > corr_all.thresh , ]
  if (nrow(stat_all) == 0){
    message("No genes retain for this corr_all.thresh. Consider decreasing the threshold.")
    return(NULL)
  }
  else {
    stat_selection = suppressWarnings( .get_gene_corr(sce, genes.selection, batch = batch, n.neigh = n.neigh,
                                                      nPC = nPC.selection, genes.predict = as.character(stat_all$gene), method = method) )
    stat_selection = merge(stat_all, stat_selection)
    stat_selection$gene_score = stat_selection$corr/stat_selection$corr_all
    return(stat_selection)
  }
}


#'
.get_gene_corr = function(sce, genes, batch = NULL, n.neigh = 5, nPC = NULL, genes.predict = rownames(sce), method = "spearman"){
  eps = 0.00001
  neighs = .get_mapping(sce , genes = genes, batch = batch , n.neigh = n.neigh , nPC = nPC)
  neighs = neighs$cells_mapped

  counts_predict = as.matrix(logcounts(sce[genes.predict , ]))

  stat_predict = lapply(1:ncol(neighs) , function(j){
    cells = neighs[,j]
    current.stat_predict = counts_predict[, cells]
    return(current.stat_predict)
  })
  stat_predict = Reduce("+", stat_predict) / length(stat_predict)
  stat_real = counts_predict[, rownames(neighs)]

  stat = lapply(1:nrow(counts_predict) , function(i){
    out = data.frame(gene = rownames(counts_predict)[i] , corr = cor(stat_real[i,] , stat_predict[i,] , method = method))
    return(out)
  })
  stat = do.call(rbind , stat)
  stat$corr[is.na(stat$corr)] = eps
  stat$corr[stat$corr < eps] = eps
  return(stat)
}
