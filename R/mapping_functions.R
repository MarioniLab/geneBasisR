# This script contains functions to construct kNN_graphs (designated for the internal use)

#'
.get_mapping = function(sce , genes = rownames(sce), batch = NULL, n.neigh = 5, nPC = 50 , cosine = F){
  if (is.null(batch)){
    out = .get_mapping_single_batch(sce , genes = genes, n.neigh = n.neigh, nPC = nPC , cosine = cosine)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    neighs = lapply(unique(batchFactor) , function(current.batch){
      #print(current.batch)
      idx = which(batchFactor == current.batch)
      current.neighs = .get_mapping_single_batch(sce[, idx] , genes = genes, n.neigh = n.neigh, nPC = nPC , cosine = cosine)
      return(current.neighs)
    })
    cells_mapped = do.call(rbind , neighs)
    cells_mapped = cells_mapped[ match(colnames(sce), rownames(cells_mapped)), ]
    out = cells_mapped
  }
  return(out)
}


#' @importFrom BiocNeighbors queryKNN
.assign_neighbors = function(counts , reference_cells , query_cells, n.neigh = 5){
  set.seed(32)
  if (is.numeric(n.neigh) & n.neigh > nrow(counts)-1){
    stop("Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.")
  }
  else {
    knns = suppressWarnings( queryKNN( counts[reference_cells ,], counts[query_cells ,], k = (n.neigh+1), get.distance = FALSE) )
    cells_mapped = lapply(1:length(query_cells), function(i){
      current.neighs = knns$index[i, ]
      current.neighs = current.neighs[!current.neighs == i]
      current.neighs = current.neighs[1:n.neigh]
      return(current.neighs)
    })
    cells_mapped = do.call(rbind, cells_mapped)
    cells_mapped = t( apply(cells_mapped, 1, function(x) reference_cells[x]) )
    rownames(cells_mapped) = query_cells
    return(cells_mapped)
  }
}




#' @importFrom irlba prcomp_irlba
#' @importFrom batchelor cosineNorm
#' @importFrom SingleCellExperiment logcounts colData
#' @import Matrix
#'
.get_mapping_single_batch = function(sce , genes = rownames(sce), n.neigh = 5, nPC = 50 , cosine = F){
  if (is.numeric(n.neigh) & n.neigh > ncol(sce)-1){
    stop("Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.")
  }
  else {
    set.seed(32)
    sce = sce[genes , ]
    if (cosine){
      logcounts(sce) = cosineNorm(logcounts(sce))
    }
    meta = as.data.frame(colData(sce))
    res = tryCatch(
      {
        counts = t(as.matrix(logcounts(sce)))
        if (!is.null(nPC)){
          n = min(nPC, (nrow(counts)-1) , (ncol(counts) - 1))
          if (n > 2){
            pcs = suppressWarnings( prcomp_irlba(counts , n = n ) )
            counts = pcs$x
          }
        }
        rownames(counts) = colnames(sce)
        out = .assign_neighbors(counts , reference_cells = colnames(sce), query_cells = colnames(sce), n.neigh = n.neigh)
        out
        #return(out)
      },
      error = function(dummy){
        #message("Count matrix is too big - we will be working with sparse matrices.")
        counts = Matrix::t(logcounts(sce))
        if (!is.null(nPC)){
          n = min(nPC, (nrow(counts)-1) , (ncol(counts) - 1))
          if (n > 2){
            pcs = suppressWarnings( prcomp_irlba(counts , n = min(nPC, (nrow(counts)-1) , (ncol(counts) - 1))) )
            counts = pcs$x
          }
        }
        rownames(counts) = colnames(sce)
        out = .assign_neighbors(counts , reference_cells = colnames(sce), query_cells = colnames(sce), n.neigh = n.neigh)
        return(out)
      },
      error = function(dump){
        message("Either memory is exhausted or features you selected can not be used for pca. Try downsampling, smaller n.neigh or smaller nPC.all.")
        return(NULL)
      }
    )
    return(res)
  }
}


#' @importFrom batchelor cosineNorm multiBatchPCA reducedMNN fastMNN
#' @importFrom SingleCellExperiment logcounts reducedDim
#'
.get_MNN_corrected_mapping = function(sce , genes = rownames(sce), batch = NULL, n.neigh = 5, nPC = 50, cosine = F){
  set.seed(32)
  if (is.null(batch)){
    out = .get_mapping(sce , genes = genes, batch = NULL, n.neigh = n.neigh, nPC = nPC , cosine = cosine)
  }
  else {
    sce = sce[genes , ]
    if (cosine){
      logcounts(sce) = cosineNorm(logcounts(sce))
    }
    res = tryCatch(
      {
        counts = as.matrix( logcounts(sce))
        meta = as.data.frame(colData(sce))
        batchFactor = factor(meta[, colnames(meta) == batch])
        if (!is.null(nPC)){

          counts = multiBatchPCA(counts , batch = batchFactor , d = nPC)

          counts = do.call(reducedMNN , counts)

          counts = counts$corrected
        } else {
          counts = fastMNN(counts , batch = batchFactor , d = NA)
          counts = reducedDim(counts , "corrected")
        }
        out = .assign_neighbors(counts , reference_cells = colnames(sce), query_cells = colnames(sce), n.neigh = n.neigh)
        out
        #return(out)
      },
      error = function(dump){
        counts = logcounts(sce)
        meta = as.data.frame(colData(sce))
        batchFactor = factor(meta[, colnames(meta) == batch])
        if (!is.null(nPC)){
          counts = multiBatchPCA(counts , batch = batchFactor , d = nPC)
          counts = do.call(reducedMNN , as.list(counts))
          counts = counts$corrected
        } else {
          counts = fastMNN(counts , batch = batchFactor , d = NA)
          counts = reducedDim(counts , "corrected")
        }
        out = .assign_neighbors(counts , reference_cells = colnames(sce), query_cells = colnames(sce), n.neigh = n.neigh)
        return(out)
      },
      error = function(dump){
        message("Something went wrong: likely memory is exhausted or features you selected can not be used for pca. Try downsampling, smaller n.neigh or smaller nPC.all.")
        return(NULL)
      }
    )
    return(res)
  }
}


#' @importFrom BiocNeighbors queryKNN
#' @importFrom SingleCellExperiment colData
.initiate_random_mapping = function(sce , batch = NULL, n.neigh = 5){
  if (is.null(batch)){
    batchFactor = factor(rep(1 , ncol(sce)))
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
  }
  initial_random_mtrx = suppressWarnings( abs(matrix(rnorm(10),2,ncol(sce))) )
  colnames(initial_random_mtrx) = colnames(sce)
  neighs = lapply(unique(batchFactor) , function(current.batch){
    idx = which(batchFactor == current.batch)
    counts = t( initial_random_mtrx[, idx] )
    reference_cells = colnames(sce[,idx])
    query_cells = colnames(sce[,idx])
    out = .assign_neighbors(counts , reference_cells, query_cells, n.neigh = n.neigh)
    return(out)
  })
  neighs = do.call(rbind , neighs)
  neighs = neighs[ match(colnames(sce), rownames(neighs)), ]
  return(neighs)
}



