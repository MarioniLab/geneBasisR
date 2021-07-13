

#' get_celltype_mapping
#'
#' For each cell, returns an estimate of its cell type based on cell type labels of neighbors in Selection graph.
#' Requires at least 2 genes for the Selection graph. If cell type mapping for the Selection graph is not possible, returns NULL.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes.selection Character vector specifying genes (at least 2 genes) to be used for the construction of Selection kNN-graph.
#' @param celltype.id Character specifying which field in colData(sce) should be used as celltype. Default celltype.id="celltype".
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Positive integer > 1, specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC.selection Positive integer (or NULL, if no PCA to be applied) specifying number of PCs to use for construction of True kNN-graph. Default nPC.selection=NULL.
#' @param cosine Boolean specifying if cosine normalization should be applied prior to constructing kNN-graph. Default cosine=FALSE.
#' @param return.stat Boolean specifying if stat for the mapping (= diagonal values for the cell type confusion matrix) should be returned alongside the mapping.
#' @param which_genes_to_use String specifying whether cell type mapping should be performed only on differentially expressed (between cell types) genes (= 'DE') or all selected genes (= 'all'). Default which_genes_to_use="all".
#' @param ... Additional arguments (e.g. the ones you can pass to get_DE_genes).
#' @return List, containing field 'mapping' - data.frame with cell IDs (column 'cell'), originally assigned cell type labels (column 'celltype') and estimated cell type labels
#' (column 'mapped_celltype'). If return.stat == TRUE, also returns field 'stat' - data.frame with cell type IDs (column 'celltype') and
#' fraction of cells from this cell type that mapped correctly (column 'frac_correctly_mapped')
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
#' genes.selection = sample( rownames(sce) , 20)
#' out = get_celltype_mapping(sce, genes.selection = genes.selection)
#'
get_celltype_mapping = function(sce , genes.selection , celltype.id = "celltype", batch = NULL, n.neigh = 5 , nPC.selection = NULL, cosine = F, return.stat = T, which_genes_to_use = "all", ...){

  # checks that inputs are eligible
  args = c(as.list(environment()) , list(...))
  if (!"check_args" %in% names(args)){
    sce = .prepare_sce(sce)
    out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.selection)
  }
  else {
    if (args[[which(names(args) == "check_args")]]){
      sce = .prepare_sce(sce)
      out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.selection)
    }
  }
  if (n.neigh > ncol(sce) - 1){
    stop("n.neigh should be less than number of cells. Decrease n.neigh.")
  }
  if (length(genes.selection) < 2){
    message("Less than 2 genes are selected for celltype mapping - celltype mapping is not possible.")
    return(NULL)
  }

  # run mapping
  else {
    sce = .update_sce_w_custom_celltype_id(sce , celltype.id = celltype.id)
    if (which_genes_to_use == "DE"){
      markers = get_DE_genes(sce, celltype.id = celltype.id, check_args = FALSE, ...)
      if (is.null(markers)){
        message("Consider set which_genes_to_use = 'all' or tune settings for DE testings.")
        return(NULL)
      }
      else {
        genes.selection = intersect(as.character(markers$gene) , as.character(genes.selection))
        if (length(genes.selection) < 2){
          message("Less than 2 genes are selected as DE between celltypes - celltype mapping is not possible. Consider option 'all' or change settings for identifying DE genes.")
          return(NULL)
        }
      }
    }
    neighs = suppressWarnings( .get_MNN_corrected_mapping(sce , genes = genes.selection, batch = batch, n.neigh = n.neigh, nPC = nPC.selection , cosine = cosine) )
    if (!is.null(neighs)){
      meta = as.data.frame(colData(sce))
      meta$celltype = as.character(meta$celltype)
      mapping = data.frame(cell = rownames(neighs) ,
                                   celltype = sapply(1:nrow(neighs) , function(i) meta$celltype[meta$cell == rownames(neighs)[i]]) ,
                                   mapped_celltype = sapply(1:nrow(neighs), function(i) .getmode(meta$celltype[match(neighs[i,] , meta$cell)] , 1:n.neigh) ))
      if (return.stat){
        stat = .get_fraction_mapped_correctly(mapping)
        final_stat = list(mapping = mapping, stat = stat)
      }
      else {
        final_stat = list(mapping = mapping)
      }
      return(final_stat)
    }
    else {
      return(NULL)
    }
  }
}


.get_fraction_mapped_correctly = function(mapping){
  tab = table(mapping$celltype , mapping$mapped_celltype)
  tab = sweep(tab, 1, rowSums(tab), "/" )
  tab = as.data.frame(tab)
  colnames(tab) = c("celltype" , "mapped_celltype" , "frac")

  celltypes = as.character(unique(tab$celltype))
  stat = lapply(celltypes, function(celltype){
    current.tab = tab[tab$celltype == celltype , ]
    if (celltype %in% current.tab$mapped_celltype){
      out = data.frame(celltype = celltype , frac_correctly_mapped = current.tab$frac[current.tab$mapped_celltype == celltype])
    }
    else {
      out = data.frame(celltype = celltype , frac_correctly_mapped = 0)
    }
    return(out)
  })
  stat = do.call(rbind, stat)
  return(stat)
}



.getmode <- function(v, dist) {
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}
