# This script contains functions to construct kNN_graphs (designated for the internal use)


#' For each cell, function assigns neighbor cells (in kNN-graph) and corresponding distances.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector containing genes that are used to construct kNN-graph.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Scalar specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC Scalar (or NULL, if no PCA applied) specifying number of PCs to use for kNN-graph. Default nPC=50.
#' @param get.dist Boolean specifying if distances for kNN-graph should be returned. Default get.dist=FALSE.
#' @param cosine Boolean specifying if cosine normalization should be applied prior to constructing kNN-graph. Default cosine=FALSE.
#'
#' @return List containing entry named 'cells_mapped': data.frame, each row corresponds to cellID (in rownames) and columns contain information about first K neighbors;
#' (if get.dist=T) entry named 'distances': corresponding to 'cells_mapped' data.frame, columns contain distances to the first K neighbors.
#' @importFrom paleotree reverseList
#'
.get_mapping = function(sce , genes = rownames(sce), batch = NULL, n.neigh = 5, nPC = 50 , get.dist = F, cosine = F){
  args = c(as.list(environment()))
  if (.check_genes_in_sce(sce , genes) & .general_check_arguments(args) & .check_batch(sce , batch)){
    if (is.null(batch)){
      out = .get_mapping_single_batch(sce , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist, cosine = cosine)
    }
    else {
      meta = as.data.frame(colData(sce))
      batchFactor = factor(meta[, colnames(meta) == batch])
      neighs = lapply(unique(batchFactor) , function(current.batch){
        idx = which(batchFactor == current.batch)
        current.neighs = .get_mapping_single_batch(sce[, idx] , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist, cosine = cosine)
        return(current.neighs)
      })
      if (!get.dist){
        neighs = reverseList(neighs)
        cells_mapped = do.call(rbind , neighs[[1]])

        cells_mapped = cells_mapped[ match(colnames(sce), rownames(cells_mapped)), ]
        out = list(cells_mapped = as.data.frame(cells_mapped))
      }
      else {
        neighs = reverseList(neighs)
        cells_mapped = do.call(rbind , neighs[[1]])
        cells_mapped = cells_mapped[ match(colnames(sce), rownames(cells_mapped)), ]
        distances = do.call(rbind , neighs[[2]])
        distances = distances[ match(colnames(sce), rownames(distances)), ]
        out = list(cells_mapped = as.data.frame(cells_mapped) , distances = as.data.frame(distances))
      }
    }
    return(out)
  }
}


#' @importFrom BiocNeighbors queryKNN
#'
.assign_neighbors = function(counts , reference_cells , query_cells, n.neigh = 5, get.dist = F){
  if (is.numeric(n.neigh) & n.neigh > nrow(counts)-1){
    stop("Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.")
  }
  else {
    knns = queryKNN( counts[reference_cells ,], counts[query_cells ,], k = (n.neigh+1), get.distance = get.dist)
    cells_mapped = t( apply(knns$index, 1, function(x) reference_cells[x[2:(n.neigh+1)]]) )
    rownames(cells_mapped) = query_cells
    if (!get.dist){
      out = list(cells_mapped = cells_mapped)
    }
    if (get.dist){
      distances = knns$distance[, 2:(n.neigh+1)]
      rownames(distances) = query_cells
      out = list(cells_mapped = cells_mapped , distances = distances)
    }
    return(out)
  }
}


#' @importFrom irlba prcomp_irlba
#'
.get_mapping_single_batch = function(sce , genes = rownames(sce), n.neigh = 5, nPC = 50 , get.dist = F , cosine = F){
  if (is.numeric(n.neigh) & n.neigh > ncol(sce)-1){
    stop("Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.")
  }
  else {
    set.seed(32)
    sce = sce[genes , ]
    if (cosine){
      logcounts(sce) = cosineNorm(logcounts(sce))
    }
    counts = as.matrix( logcounts(sce) )
    meta = as.data.frame(colData(sce))
    res = tryCatch(
      {
        if (!is.null(nPC)){
          pcs = prcomp_irlba(t(counts) , n = min(nPC, (nrow(counts)-1) , (ncol(counts) - 1)))
          counts = pcs$x
          rownames(counts) = colnames(sce)
        } else {
          counts = t(counts)
        }
        reference_cells = colnames(sce)
        query_cells = colnames(sce)
        if (n.neigh == "all"){
          n.neigh = length(reference_cells) - 1
        }
        out = .assign_neighbors(counts , reference_cells, query_cells, n.neigh = n.neigh, get.dist = get.dist)
        return(out)
      },
      error = function(dump){
        message("Features you selected can not be used for pca")
        return(NULL)
      }
    )
    return(res)
  }
}


#' For each cell, function assigns neighbor cells (in MNN-correced kNN-graph) and corresponding distances.
#'
#' @param sce SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).
#' @param genes Character vector containing genes that are used to construct kNN-graph.
#' @param batch Name of the field in colData(sce) to specify batch. Default batch=NULL if no batch is applied.
#' @param n.neigh Scalar specifying number of neighbors to use for kNN-graph. Default n.neigh=5.
#' @param nPC Scalar specifying number of PCs to use for kNN-graph. Default nPC=50.
#' @param cosine Boolean specifying if cosine normalization should be applied prior to constructing kNN-graph. Default cosine=FALSE.
#'
#' @import batchelor
#'
.get_MNN_corrected_mapping = function(sce , genes = rownames(sce), batch = NULL, n.neigh = 5, nPC = 50, cosine = F){
  sce = .prepare_sce(sce)
  args = c(as.list(environment()))
  if (.check_genes_in_sce(sce , genes) & .general_check_arguments(args) & .check_batch(sce , batch)){
    if (is.null(batch)){
      out = .get_mapping(sce , genes = genes, batch = NULL, n.neigh = n.neigh, nPC = nPC , cosine = cosine)
    }
    else {
      sce = sce[genes , ]
      if (cosine){
        logcounts(sce) = cosineNorm(logcounts(sce))
      }
      counts = as.matrix( logcounts(sce))
      meta = as.data.frame(colData(sce))
      batchFactor = factor(meta[, colnames(meta) == batch])
      res = tryCatch(
        {
          if (!is.null(nPC)){
            counts = multiBatchPCA(counts , batch = batchFactor , d = nPC)
            counts = do.call(reducedMNN , counts)
            counts = counts$corrected
          } else {
            counts = fastMNN(counts , batch = batchFactor , d = NA)
            counts = reducedDim(counts , "corrected")
          }
          reference_cells = colnames(sce)
          query_cells = colnames(sce)
          if (n.neigh == "all"){
            n.neigh = length(reference_cells) - 1
          }
          out = .assign_neighbors(counts , reference_cells , query_cells, n.neigh = n.neigh)
          return(out)
        },
        error = function(dump){
          message("Features you use are insufficient to perform batch correction via fastMNN")
          return(NULL)
        }
      )
      return(res)
    }
  }
}


