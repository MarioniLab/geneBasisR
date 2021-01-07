#' Filtering for highly expressed genes.
#'
#' @param sce SingleCellExperiment object representing scRNA-seq counts matrix containing celltype field.
#' @param method String specifying the approach to take to detect highly expressed genes. If set to "absolute", all genes that have q75 higher than a threshold specified in logcounts.thresh (in at least one celltype) will be discarded. If set to "ribosome_based", the logcounts.thresh is calculated as minimum of upper bounds of q75 across all stably expressed ribosomal genes.
#' @param logcounts.thresh Numeric positive scalar specifying the threshold to detect highly expressed genes.
#'
#' @return Filtered counts matrix (SingleCellExperiment object).
#' @export
#' @importFrom stats quantile
#' @importFrom dplyr summarise group_by
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 30000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$celltype = as.factor(sample(1:5, n_col, replace=TRUE))
#' sce = filter_highly_expressed_genes(sce)
#' sce = filter_highly_expressed_genes(sce, "ribosome_based")
#'
filter_highly_expressed_genes = function(sce, method = "absolute", logcounts.thresh = 20){
  if (!method %in% c("absolute", "ribosome_based")){
    stop("method variable should be either absolute or ribosome.")
  }
  else {
    if (!.check_counts_matrix_correct(sce)) {
      stop()
    } else{
      stat = lapply(unique(sce$celltype) , function(celltype){
        counts = logcounts(sce[, sce$celltype == celltype])
        current.stat = apply(counts , 1 , function(x) quantile(x , 0.75))
        current.stat = data.frame(gene = names(current.stat) , q75 = current.stat , celltype = celltype)
        return(current.stat)
      })
      stat = do.call(rbind, stat)
      stat = as.data.frame(stat %>% group_by(gene) %>% summarise(max.q75 = max(q75)))
      if (method == "ribosome_based"){
        stably_expressed_ribosomal_genes = geneBasisR:::stably_expressed_ribosomal_genes
        message(paste0(length(stably_expressed_ribosomal_genes) , " genes in total are assigned as stably expressed ribosomal genes."))
        stat$gene = toupper(stat$gene)
        stat.stably_expressed_ribosomal_genes = stat[stat$gene %in% stably_expressed_ribosomal_genes, ]
        message(paste0(nrow(stat.stably_expressed_ribosomal_genes) , " stably expressed ribosomal genes are found in counts matrix."))
        logcounts.thresh = min(stat.stably_expressed_ribosomal_genes$max.q75)
      }
      genes = as.character( stat$gene[stat$max.q75 < logcounts.thresh] )
      message(paste0(length(genes) , " genes are retained as not highly expressed."))
      sce.filtered = sce[rownames(sce) %in% genes , ]
      return(sce.filtered)
    }
  }
}

