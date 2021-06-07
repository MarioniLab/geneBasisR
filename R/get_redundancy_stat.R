

#' get_redundancy_stat
#'
#' Functions calculates relevance of each gene (within current selection) to celltype mapping.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector specifying selected library.
#' @param genes_to_assess Character vector specifying gene names for which we want to assess redundancy. Should be a subset of genes.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param ... Additional arguments you can pass to get_celltype_mapping.
#'
#' @return data.frame, each row corresponds to celltype/gene; frac_correctly_mapped_ratio corresponds to ratio between accuracy of the mapping for a given celltype while using genes excluding and including the gene in question.
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
#' sce$celltype = as.character(sample.int(5, n_col, replace = TRUE))
#' genes = rownames(sce)
#' genes_to_assess = sample(rownames(sce),5)
#' out = get_redundancy_stat(sce, genes, genes_to_assess = genes_to_assess)
#'
get_redundancy_stat = function(sce, genes, genes_to_assess = genes, batch = NULL, ...){
  args = c(as.list(environment()) , list(...))
  sce = .prepare_sce(sce)
  out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes) & .check_genes_in_sce(sce , genes_to_assess)

  mapping_all = get_celltype_mapping(sce , genes , batch = batch, return.stat = TRUE, check_args = FALSE, ...)
  stat_all = mapping_all$stat
  colnames(stat_all) = c("celltype", "frac_correctly_mapped_all")

  stat_reduced = lapply(genes_to_assess, function(gene){
    current.mapping = get_celltype_mapping(sce , setdiff(genes,gene) , batch = batch, return.stat = T, check_args = FALSE, ...)
    current.stat = current.mapping$stat
    current.stat$gene = gene
    return(current.stat)
  })
  stat_reduced = do.call(rbind, stat_reduced)
  stat_reduced = merge(stat_reduced, stat_all, all.x = T)
  stat_reduced$frac_correctly_mapped_ratio = stat_reduced$frac_correctly_mapped/stat_reduced$frac_correctly_mapped_all
  return(stat_reduced)
}


