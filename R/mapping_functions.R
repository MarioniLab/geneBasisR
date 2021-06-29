# This script contains functions to construct kNN_graphs (designated for the internal use)

#' @importFrom paleotree reverseList
#'
.get_mapping = function(sce , genes = rownames(sce), batch = NULL, n.neigh = 5, nPC = 50 , get.dist = F, cosine = F){
  if (is.null(batch)){
    out = .get_mapping_single_batch(sce , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist, cosine = cosine)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    neighs = lapply(unique(batchFactor) , function(current.batch){
      #print(current.batch)
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


#' @importFrom BiocNeighbors queryKNN
#' @importFrom paleotree reverseList
.assign_neighbors = function(counts , reference_cells , query_cells, n.neigh = 5, get.dist = F ){
  set.seed(32)
  if (is.numeric(n.neigh) & n.neigh > nrow(counts)-1){
    stop("Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.")
  }
  else {
    knns = suppressWarnings( queryKNN( counts[reference_cells ,], counts[query_cells ,], k = (n.neigh+1), get.distance = get.dist) )
    if (!get.dist){
      cells_mapped = lapply(1:length(query_cells), function(i){
        current.neighs = knns$index[i, ]
        current.neighs = current.neighs[!current.neighs == i]
        current.neighs = current.neighs[1:n.neigh]
        return(current.neighs)
      })
      cells_mapped = do.call(rbind, cells_mapped)
      cells_mapped = t( apply(cells_mapped, 1, function(x) reference_cells[x]) )
      rownames(cells_mapped) = query_cells
      out = list(cells_mapped = cells_mapped)
    }
    else {
      cells_mapped_w_dist = lapply(1:nrow(knns$index), function(i){
        current.neighs = knns$index[i, ]
        idx = which(!current.neighs == i)
        current.neighs = current.neighs[idx]
        current.neighs = current.neighs[1:n.neigh]
        current.dist = knns$distance[i, idx]
        current.dist = current.dist[1:n.neigh]
        return(list(cells_mapped = current.neighs , distances = current.dist))
      })
      cells_mapped_w_dist = reverseList(cells_mapped_w_dist)
      cells_mapped = do.call(rbind , cells_mapped_w_dist$cells_mapped)
      cells_mapped = t( apply(cells_mapped, 1, function(x) reference_cells[x]) )
      rownames(cells_mapped) = query_cells
      distances = do.call(rbind , cells_mapped_w_dist$distances)
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

    meta = as.data.frame(colData(sce))
    res = tryCatch(
      {
        counts = t(as.matrix(logcounts(sce)))
        if (!is.null(nPC)){
          pcs = suppressWarnings( prcomp_irlba(counts , n = min(nPC, (nrow(counts)-1) , (ncol(counts) - 1))) )
          counts = pcs$x
        }
        rownames(counts) = colnames(sce)
        reference_cells = colnames(sce)
        query_cells = colnames(sce)
        if (n.neigh == "all"){
          n.neigh = length(reference_cells) - 1
        }
        out = .assign_neighbors(counts , reference_cells, query_cells, n.neigh = n.neigh, get.dist = get.dist)
        print("neighs done")
        return(out)
      },
      error = function(dummy){
        message("Count matrix is too big - we will be working with sparse matrices.")
        counts = t(logcounts(sce))
        if (!is.null(nPC)){
          pcs = suppressWarnings( prcomp_irlba(counts , n = min(nPC, (nrow(counts)-1) , (ncol(counts) - 1))) )
          counts = pcs$x
        }
        rownames(counts) = colnames(sce)
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


#' @import batchelor
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



#' @importFrom BiocNeighbors queryKNN
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
    out = .assign_neighbors(counts , reference_cells, query_cells, n.neigh = n.neigh, get.dist = F)
    return(out$cells_mapped)
  })
  neighs = do.call(rbind , neighs)
  neighs = neighs[ match(colnames(sce), rownames(neighs)), ]
  return(neighs)
}

