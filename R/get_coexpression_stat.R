
#' get_coexpression_stat
#'
#' For each gene of interest, returns top co-expressed genes from the sce counts matrix (either filtered by number of co-expressed genes or correlation)
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Charcter vector specifying genes to which co-expression stat should be returned.
#' @param n Scalar (or NULL) specifying number of top HVGs to be return. Default n=10. If NULL, returns all genes with correlation higher than threshold.
#' @param corr.thresh Scalar (float) specifying threshold for correlation to use. Default corr.thresh=0.5.
#' @param method Character specifying method for correlation. Availbale options are c("spearman", "pearson", "kendall"). Default method="spearman".
#' @param verbose Boolean identifying whether intermediate print outputs should be returned. Default verbose=TRUE.
#' @return List, each element corresponds to a gene in genes. Each element is data.frame with co-expressed genes.
#' @export
#' @import SingleCellExperiment
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 500
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' out = get_coexpression_stat(sce, c("1","2"))
#'
get_coexpression_stat = function(sce, genes, n = 10, corr.thresh = 0.5, method = "spearman", verbose = TRUE){
  args = c(as.list(environment()))
  sce = .prepare_sce(sce)
  if (.general_check_arguments(args)){
    stat = lapply(genes, function(gene){
      if (verbose){
        cat(paste0("Calculating co-expression stat for gene ", gene, ".\n"))
      }
      current.stat = .get_coexpression_stat_single_gene(sce, gene, n = n, corr.thresh = corr.thresh, method = method)
    })
    names(stat) = genes
    return(stat)
  }
}

.get_coexpression_stat_single_gene = function(sce, gene, n = 10 , corr.thresh = 0.5, method = "spearman"){

  res = tryCatch(
    {
      counts = as.matrix( logcounts(sce))
      counts.gene = counts[gene, ]
      stat = lapply(1:nrow(counts), function(i){
        out = data.frame(corr = cor(counts.gene, counts[i,] , method = method))
      })
      stat = do.call(rbind , stat)
      stat$gene = rownames(counts)
      stat = stat[!stat$gene == gene, ]
      stat = stat[order(stat$corr, decreasing = T),]
      if (is.null(n)){
        stat = stat[stat$corr >= corr.thresh, ]
      }
      else {
        stat = stat[1:n,]
      }
      stat
    },
    error = function(dump){
      counts = logcounts(sce)
      counts.gene = counts[gene, ]
      stat = lapply(1:nrow(counts), function(i){
        out = data.frame(corr = cor(counts.gene, counts[i,] , method = method))
      })
      stat = do.call(rbind , stat)
      stat$gene = rownames(counts)
      stat = stat[!stat$gene == gene, ]
      stat = stat[order(stat$corr, decreasing = T)]
      if (is.null(n)){
        stat = stat[stat$corr >= corr.thresh, ]
      }
      else {
        stat = stat[1:n,]
      }
      return(stat)
    },
    error = function(dump){
      message("Something went wrong: likely memory is exhausted.")
      return(NULL)
    }
  )
  return(res)

}
