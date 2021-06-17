

#' evaluate_library
#'
#' A wrapper to return the estimates of the library quality (at cell, gene and/or celltype levels) as a function of number of genes.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes.selection Character vector specifying genes to be used for the construction of Selection kNN-graph.
#' @param genes.all Character vector specifying genes to be used for the construction of True kNN-graph.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Positive integer > 1, specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param library.size_type String identifying whether evaluation should be performed only on the whole inserted library (= 'single') or on a series of subsets of the library (= 'series'). Default library.size_type="single".
#' @param n_genes.step In case library.size_type == "series", a scalar identifying the step of the grid for library subsets. Default n_genes.step=10.
#' @param return.cell_score_stat Boolean identifying whether stat on cell neighborhood preservation score should be returned. Default return.cell_score_stat=TRUE.
#' @param return.gene_score_stat Boolean identifying whether stat on gene prediction score should be returned. Default return.gene_score_stat=TRUE.
#' @param return.celltype_stat Boolean identifying whether stat on celltype mapping should be returned. Default return.celltype_stat=TRUE.
#' @param verbose Boolean identifying whether intermediate print outputs should be returned. Default verbose=TRUE.
#' @param neighs.all If not NULL (NULL is default), contains information about True kNN-graph (for each cell - ordered neighbors and distances). Useful to have a priori if cell score will be recalculated multiple times.
#' @param gene_stat_all If not NULL (NULL is default), gene_stat_all is pre-calculated stat for True graph. This is useful if this variable will be used re-used multiple times.
#' @param ... Additional parameters
#'
#' @return
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
#' sce$celltype = as.character(sample.int(5, n_col, replace = TRUE))
#' genes.selection = sample(rownames(sce) , 20)
#' out = evaluate_library(sce, genes.selection)
#'
evaluate_library = function(sce, genes.selection, genes.all = rownames(sce), batch = NULL, n.neigh = 5,
                            library.size_type = "single" , n_genes.step = 10,
                            return.cell_score_stat = T, return.gene_score_stat = T, return.celltype_stat = T, verbose = TRUE,
                            neighs.all = NULL, gene_stat_all = NULL, ...){
  sce = .prepare_sce(sce)
  args = c(as.list(environment()), list(...))
  out = .general_check_arguments(args) & .check_batch(sce , batch) & .check_genes_in_sce(sce , genes.selection) & .check_genes_in_sce(sce, genes.all)
  if (return.celltype_stat){
    out = .check_celltype_in_sce(sce)
  }

  if (library.size_type == "single"){
    n_genes.grid = c(length(genes.selection))
  } else if (library.size_type == "series"){
    n_genes.grid = c( seq(n_genes.step , length(genes.selection) , n_genes.step) , length(genes.selection))
    n_genes.grid = sort(unique(n_genes.grid))
  }

  if (sum(return.cell_score_stat, return.gene_score_stat, return.celltype_stat) == 0){
    stop("Select at least one stat to return.")
    return(F)
  }
  else {
    final_stat = list()
    names(final_stat) = c()


    #### cell score stat
    if (return.cell_score_stat){
      if (verbose){
        cat("Calculating cell neighborhood preservation scores.\n")
      }
      if (!is.null(neighs.all)){
        if (!is.null(batch)){
          out = .check_neighs.all_multipleBatches(sce , batch = batch , neighs.all = neighs.all)
        }
      }
      else {
        neighs.all = suppressWarnings( get_z_scaled_distances(sce , genes.all = genes.all , batch = batch, n.neigh = n.neigh, ...) )
      }
      cell_score_stat = lapply(n_genes.grid, function(n_genes){
        current.stat = suppressWarnings( get_neighborhood_preservation_scores(sce, neighs.all = neighs.all,  genes.all = genes.all,
                                                                       genes.selection = genes.selection[1:n_genes], batch = batch, n.neigh = n.neigh, check_args = FALSE, ...) )
        current.stat$n_genes = n_genes
        if (verbose){
          cat(paste("Finished for the selection of" , n_genes , "genes.\n"))
        }
        return(current.stat)
      })
      cell_score_stat = do.call(rbind , cell_score_stat)
      cell_score_stat$n_genes = factor(cell_score_stat$n_genes, levels = sort(unique(cell_score_stat$n_genes)))
      final_stat[[length(final_stat) + 1]] = cell_score_stat
      names(final_stat)[length(final_stat)] = "cell_score_stat"
      if (verbose){
        cat("Finished calculation of cell neighborhood preservation scores.\n")
      }
    }

    #### gene score stat
    if (return.gene_score_stat){
      if (verbose){
        cat("Calculating gene prediction scores.\n")
      }
      if (is.null(gene_stat_all)){
        gene_stat_all = suppressWarnings( get_gene_correlation_scores(sce, genes.all, batch = batch, n.neigh = n.neigh, ...) )
        colnames(gene_stat_all) = c("gene" , "corr_all")
      }
      gene_score_stat = lapply(n_genes.grid, function(n_genes){
        current.stat = suppressWarnings( get_gene_prediction_scores(sce, genes.selection[1:n_genes], genes.all = genes.all, batch = batch,
                                                                    n.neigh = n.neigh, gene_stat_all = gene_stat_all, check_args = FALSE, ...) )
        if (!is.null(current.stat)){
          current.stat$n_genes = n_genes
        }
        if (verbose){
          cat(paste("Finished for the selection of" , n_genes , "genes.\n"))
        }
        return(current.stat)
      })
      gene_score_stat = do.call(rbind , gene_score_stat)
      if (!is.null(gene_score_stat)){
        gene_score_stat$n_genes = factor(gene_score_stat$n_genes, levels = sort(unique(gene_score_stat$n_genes)))
        final_stat[[length(final_stat) + 1]] = gene_score_stat
        names(final_stat)[length(final_stat)] = "gene_score_stat"
      }
      if (verbose){
        cat("Finished calculation of gene prediction scores.\n")
      }
    }


    ### celltype stat
    if (return.celltype_stat){
      if (verbose){
        cat("Calculating accuracy of cell type mappings.\n")
      }
      celltype_stat = lapply(n_genes.grid, function(n_genes){
        current.stat = get_celltype_mapping(sce , genes.selection[1:n_genes] , batch = batch, n.neigh = n.neigh, return.stat = T, check_args = FALSE, ...)
        if (!is.null(current.stat)){
          current.stat = current.stat$stat
          current.stat$n_genes = n_genes
        }
        if (verbose){
          cat(paste("Finished for the selection of" , n_genes , "genes.\n"))
        }
        return(current.stat)
      })
      celltype_stat = do.call(rbind , celltype_stat)
      if (!is.null(celltype_stat)){
        celltype_stat$n_genes = factor(celltype_stat$n_genes, levels = sort(unique(celltype_stat$n_genes)))
        final_stat[[length(final_stat) + 1]] = celltype_stat
        names(final_stat)[length(final_stat)] = "celltype_stat"
      }
      if (verbose){
        cat("Finished calculation of accuracy of cell type mappings.\n")
      }
    }
    return(final_stat)
  }
}


