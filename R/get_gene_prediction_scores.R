

#' get_gene_prediction_scores
#'
#' Calculates gene prediction scores by comparing True and Selection k-NN graphs.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes.selection Character vector containing names of selected genes.
#' @param genes.all Character vector containing names of all (potentially relevant/variable) genes.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Positive integer > 1, specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC.all Scalar (or NULL) specifying number of PCs to use for construction of True kNN-graph. Default nPC.all=50.
#' @param nPC.selection Scalar (or NULL) specifying number of PCs to use for construction of Selection kNN-graph. Default nPC.selection=NULL (no PCA).
#' @param genes.predict Character vector containing names of genes for which we want to calculate gene prediction score. Default = genes.all.
#' @param method Character specifying method for correlation. Availbale options are c("spearman", "pearson", "kendall"). Default method="spearman".
#' @param corr_all.thresh Scalar specifying suitable threshold for correlation to consider (on True graph).
#' @param gene_stat_all If not NULL (NULL is default), gene_stat_all is pre-calculated stat for True graph. This is useful if this variable will be used re-used multiple times. Use geneBasisR::get_gene_correlation_scores to calculate this.
#' @param ... Additional arguments
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
                                      genes.predict = rownames(sce), method = "spearman", corr_all.thresh = 0.25 , gene_stat_all = NULL, ...){
  genes.predict = intersect(genes.all, genes.predict)
  args = c(as.list(environment()) , list(...))
  if (!"check_args" %in% names(args)){
    sce = .prepare_sce(sce)
    out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.selection) & .check_genes_in_sce(sce , genes.all) & .check_genes_in_sce(sce, genes.predict)
  }
  else {
    if (args[[which(names(args) == "check_args")]]){
      sce = .prepare_sce(sce)
      out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.selection) & .check_genes_in_sce(sce , genes.all) & .check_genes_in_sce(sce, genes.predict)
    }
  }

  if (is.null(gene_stat_all)){
    gene_stat_all = suppressWarnings( get_gene_correlation_scores(sce, genes.all, batch = batch, n.neigh = n.neigh,
                                                nPC = nPC.all, genes.predict = genes.predict, method = method, check_args = FALSE) )
    colnames(gene_stat_all) = c("gene" , "corr_all")
  }
  else {
    if (!sum(c("gene" , "corr_all" ) %in% colnames(gene_stat_all)) == 2){
      stop("gene_stat_all is of the wrong format - should contain fields gene and corr_all.")
      return(FALSE)
    }
  }

  gene_stat_all = gene_stat_all[gene_stat_all$corr_all > corr_all.thresh & gene_stat_all$gene %in% genes.predict, ]
  if (nrow(gene_stat_all) == 0){
    message("No genes are retained for this corr_all.thresh. Consider decreasing the threshold.")
    return(NULL)
  }
  else {
    stat_selection = suppressWarnings( get_gene_correlation_scores(sce, genes.selection, batch = batch, n.neigh = n.neigh,
                                                      nPC = nPC.selection, genes.predict = as.character(gene_stat_all$gene), method = method, check_args = FALSE) )
    stat_selection = merge(gene_stat_all, stat_selection)
    stat_selection$gene_score = stat_selection$corr/stat_selection$corr_all
    return(stat_selection)
  }
}


#' get_gene_correlation_scores
#'
#' Calculates per gene correlation between measured expression levels and estimated expression levels (average across neighbors).
#' It is an intermediate step for 'get_gene_prediction_scores' function. It can be handy to calculate this part separately to be recycled multiple times
#' to compare different seelctions.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector containing names of selected genes.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Positive integer > 1, specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC Scalar (or NULL) specifying number of PCs to use for construction of kNN-graph. Default nPC=NULL.
#' @param genes.predict Character vector containing names of genes for which we want to calculate gene prediction score. Default = rownames(sce).
#' @param method Character specifying method for correlation. Availbale options are c("spearman", "pearson", "kendall"). Default method="spearman".
#' @param ... Additional arguments.
#'
#' @return
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
#' genes = rownames(sce)
#' out = get_gene_correlation_scores(sce, genes)
#'
get_gene_correlation_scores = function(sce, genes, batch = NULL, n.neigh = 5, nPC = NULL, genes.predict = rownames(sce), method = "spearman", ...){
  args = c(as.list(environment()) , list(...))
  if (!"check_args" %in% names(args)){
    sce = .prepare_sce(sce)
    out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes) & .check_genes_in_sce(sce, genes.predict)
  }
  else {
    if (args[[which(names(args) == "check_args")]]){
      sce = .prepare_sce(sce)
      out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes) & .check_genes_in_sce(sce, genes.predict)
    }
  }
  if (length(genes.predict) < 1){
    stop("Select at least one entry for genes.predict.")
  }
  eps = 0.00001
  neighs = .get_mapping(sce , genes = genes, batch = batch , n.neigh = n.neigh , nPC = nPC)

  res = tryCatch(
    {
      counts_predict = as.matrix(logcounts(sce[genes.predict , ]))
      stat_predict = matrix( 0L, nrow = dim(counts_predict)[1], ncol = dim(counts_predict)[2])
      for (j in c(1:n.neigh)){
        cells = neighs[,j]
        stat_predict = stat_predict + counts_predict[, cells]
      }
      stat_predict = stat_predict / n.neigh
      stat_real = counts_predict[, rownames(neighs)]

      if (nrow(counts_predict) > 1){
        stat = lapply(1:nrow(counts_predict) , function(i){
          out = data.frame(gene = rownames(counts_predict)[i] , corr = cor(stat_real[i,] , stat_predict[i,] , method = method))
          return(out)
        })
        stat = do.call(rbind , stat)
      }
      else if (nrow(counts_predict) == 1){
        stat = data.frame(gene = rownames(counts_predict)[1] , corr = cor(as.numeric(stat_real) , as.numeric(stat_predict) , method = method))
      }
      stat$corr[is.na(stat$corr)] = eps
      stat$corr[stat$corr < eps] = eps
      stat
      #return(stat)
    },
    error = function(dump){
      #message("Count matrix is too big - we will be working with sparse matrices.")
      counts_predict = logcounts(sce[genes.predict , ])
      stat_predict = counts_predict
      stat_predict[!stat_predict == 0] = 0

      for (j in c(1:n.neigh)){
        cells = neighs[,j]
        stat_predict = stat_predict + counts_predict[, cells]
      }
      stat_predict = stat_predict / n.neigh
      stat_real = counts_predict[, rownames(neighs)]
      if (nrow(counts_predict) > 1){
        stat = lapply(1:nrow(counts_predict) , function(i){
          out = data.frame(gene = rownames(counts_predict)[i] , corr = cor(stat_real[i,] , stat_predict[i,] , method = method))
          return(out)
        })
        stat = do.call(rbind , stat)
      }
      else if (nrow(counts_predict) == 1){
        stat = data.frame(gene = rownames(counts_predict)[1] , corr = cor(as.numeric(stat_real) , as.numeric(stat_predict) , method = method))
      }
      stat$corr[is.na(stat$corr)] = eps
      stat$corr[stat$corr < eps] = eps
      return(stat)
    },
    error = function(dump){
      message("Something went wrong: likely memory is exhausted or features you selected can not be used for pca. Try downsampling, smaller n.neigh or smaller nPC.")
      return(NA)
    }
  )
  return(res)
}
