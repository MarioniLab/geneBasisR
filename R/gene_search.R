

#' gene_search
#'
#' Main function of the package - returns optimal library of the selected size.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes_base Character vector specifying base genes to construct first Selection graph. Default=NULL in case no genes are supplied.
#' @param n_genes_total Scalar specifying total number of genes to be selected (this includes base genes).
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Positive integer > 1, specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param p.minkowski Order of Minkowski distance. Default p.minkowski=3.
#' @param nPC.selection Scalar specifying number of PCs to use for Selection Graphs. Default nPC=NULL.
#' @param nPC.all Scalar specifying number of PCs to use for True Graph. Default nPC.all=50.
#' @param genes.discard Character vector containing genes to be excluded from candidates (note that they still will be used for graphs construction. If you want to exclude them from graph construction as well, just discard them prior in sce object). Default = NULL and no genes will be discarded.
#' @param genes.discard_prefix Character vector containing prefixes of genes to be excluded (e.g. Rpl for L ribosomal proteins. Note that they still will be used for graphs construction. If you want to exclude them from graph construction as well, just discard them prior in sce object). Default = NULL and no genes will be discarded.
#' @param verbose Boolean identifying whether intermediate print outputs should be returned. Default verbose=TRUE.
#' @param stat_all If True graph and corresponding Minkowski distances have been calculated prior to search, provide this data here.
#' It can be useful if gene_search is desired to be recycled (e.g. for selecting multiple libraries with different inputs such as n_genes_total and genes_base)
#' Ensure that colnames = c("gene", "dist_all"). Default stat_all=NULL - in case this info is not supplied.
#'
#' @return data.frame containing selected genes and corresponding ranks. In case genes_base are supplied, rank among them will be assigned based on the order they are supplied in the corresponding string.
#' @export
#' @importFrom gdata startsWith
#' @import SingleCellExperiment
#' @examples
#' require(SingleCellExperiment)
#' n_row = 1000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' genes = rownames(sce)
#' out = gene_search(sce, n_genes_total = 5)
#'
gene_search = function(sce , genes_base = NULL, n_genes_total , batch = NULL, n.neigh = 5, p.minkowski = 3,
                       nPC.selection = NULL, nPC.all = 50, genes.discard = NULL, genes.discard_prefix = NULL, verbose = TRUE, stat_all = NULL){
  sce = .prepare_sce(sce)
  args = c(as.list(environment()))
  out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes_base)

  # get together genes to discard
  if (!is.null(genes.discard_prefix)){
    rownames.sce = rownames(sce)
    idx = sapply(1:nrow(sce) , function(i) max(startsWith(rownames.sce[i] , genes.discard_prefix)))
    idx = which(idx == 1)
    genes.discard = unique(c(genes.discard , rownames.sce[idx]))
  }

  if (n_genes_total >= nrow(sce)){
    stop("Selected library size should be smaller than number of genes in the counts matrix.")
  }
  else if (n_genes_total >= length(setdiff(rownames(sce) , genes.discard))){
    stop("Selected library size should be smaller than number of non-discarded genes in counts matrix. Reduce n_genes_total or list of genes to be discraded.")
  }
  else {
    if (verbose){
      cat("Constructing the True graph.\n")
    }

    # get baseline stat-all
    if (is.null(stat_all)){
      stat_all = suppressWarnings( calc_Minkowski_distances(sce, genes = rownames(sce), batch = batch, n.neigh = n.neigh, nPC = nPC.all,
                                                      genes.predict = rownames(sce) , p.minkowski = p.minkowski,
                                                      genes.discard = genes.discard, genes.discard_prefix = genes.discard_prefix, check_args = FALSE) )
      colnames(stat_all) = c("gene" , "dist_all")
      if (verbose){
        cat("True graph is constructed.\n")
      }
    }
    else {
      if (!sum(c("gene" , "dist_all" ) %in% colnames(stat_all)) == 2){
        stop("stat_all is of the wrong format - should contain fields gene and dist_all.")
      }
    }
    # add first gene if selection is empty
    if (is.null(genes_base)){
      K = 5
      genes_base = .add_first_gene(sce , stat_all, batch = batch , n.neigh = n.neigh, p.minkowski = p.minkowski,
                                   genes.discard = genes.discard, genes.discard_prefix = genes.discard_prefix, K = K)
      if (verbose){
        cat(paste0("First gene is added: ", genes_base , ". " , n_genes_total - 1, " left.\n"))
      }
    }
    genes_all = genes_base
    while(length(genes_all) < n_genes_total){
      gene = .add_gene_to_current_selection(sce , stat_all, genes = genes_all , batch = batch , n.neigh = n.neigh, nPC = nPC.selection,
                                            p.minkowski = p.minkowski, genes.discard = genes.discard, genes.discard_prefix = genes.discard_prefix)
      if (!is.null(gene)){
        genes_all = c(genes_all , as.character(gene))
        if (verbose){
          cat(paste0("New gene is added: ", as.character(gene) , ". " , n_genes_total - length(genes_all) , " left.\n"))
        }
      }
      else {
        break
      }
    }
    out = data.frame(rank = c(1:length(genes_all)) , gene = genes_all)
    return(out)
  }
}


.add_first_gene = function(sce , stat_all, batch = NULL , n.neigh = 5, p.minkowski = 3 , genes.discard = NULL, genes.discard_prefix = NULL, K = 5){
  first_genes = lapply(1:K , function(dump){
    gene = .add_gene_to_current_selection(sce , stat_all, genes = NULL , batch = batch , n.neigh = n.neigh, nPC = NULL,
                                          p.minkowski = p.minkowski, genes.discard = genes.discard, genes.discard_prefix = genes.discard_prefix)
  })
  first_genes = unlist(first_genes)
  out = as.character(names(sort(table(first_genes),decreasing=TRUE)[1]))
  return(out)
}

.add_gene_to_current_selection = function(sce , stat_all, genes = NULL , batch = NULL , n.neigh = 5, nPC = NULL,
                                          p.minkowski = 3, genes.discard = NULL, genes.discard_prefix = NULL){
  stat_genes = suppressWarnings( calc_Minkowski_distances(sce , genes = genes , batch = batch, n.neigh = n.neigh , nPC = nPC,
                                                  genes.predict = rownames(sce) , p.minkowski = p.minkowski,
                                                  genes.discard = genes.discard, genes.discard_prefix = genes.discard_prefix, check_args = FALSE) )
  stat_genes = stat_genes[!stat_genes$gene %in% genes , ]
  if (nrow(stat_genes) >= 1){
    stat_genes = merge(stat_genes , stat_all)
    stat_genes$dist_diff = stat_genes$dist - stat_genes$dist_all
    idx = which(stat_genes$dist_diff == max(stat_genes$dist_diff))
    gene = stat_genes$gene[idx[1]]
    return(gene)
  }
  else {
    message("No genes are left to be added.")
    return(NULL)
  }
}
