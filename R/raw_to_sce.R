


#' raw_to_sce
#'
#' Simply a parser raw .txt files --> SingleCellExperiment object + adding logcounts
#'
#' @param counts_dir String specifying the directory for counts matrix (assuming counts where already calculated)
#' @param counts_type String specifying whether raw data is stored as counts or log-counts. For geneBasis we recommend to work with log-counts.
#' @param transform_counts_to_logcounts In case, raw data are counts, Boolean specifying whether we should perform log-normalization.
#' @param header Boolean specifying if logcounts_dir file has cell IDs stored in colnames.
#' @param sep the field separator string. Note it should be the same for logcounts_dir and meta_dir (if exists)
#' @param meta_dir If not NULL, a string specifying the directory for meta-data (i.e. celltype, batch, UMAP-coordinates).
#' @param batch If not NULL (no batch applied), string specifying a column in meta file that will be used as batchID. Please check that named column exists in meta-file.
#' Note that UMAP coordinates should be named as x and y. Also, if meta contains field cell - this field will be used for cell IDs.
#'
#' @return SingleCellExperiment object with gene counts/logcounts and meta-data (if applicable)
#' @export
#' @import SingleCellExperiment
raw_to_sce = function(counts_dir, counts_type = "counts", transform_counts_to_logcounts = TRUE, header = TRUE, sep = "\t" , meta_dir = NULL, batch = NULL,...){
  if (!file.exists(counts_dir)){
    stop("Counts file does not exist")
  }
  else if (!is.null(meta_dir)) {
    if (!file.exists(meta_dir)) {
      stop("Meta file does not exist")
    }
  }
  else {
    counts = read.table(counts_dir, header = header, sep = sep)
    if (!header){
      colnames(counts) = c(1:ncol(sce))
    }
    if (!is.null(meta_dir)){
      meta = read.table(meta_dir, header = TRUE, sep = sep)
      if (!"cell" %in% colnames(meta)){
        meta$cell = c(1:nrow(meta))
      }
      meta = meta[order(meta$cell) , ]
      counts = counts[, order(colnames(counts))]
      if (nrow(meta) != ncol(counts)){
        stop("Mismatch in number of cells between counts matrix and meta-file.")
        return(FALSE)
      } else if (mean(meta$cell == colnames(counts)) < 1){
        stop("Mismatch in cell IDs between counts matrix and meta-file.")
        return(FALSE)
      } else {
        sce = SingleCellExperiment(assay = counts, colData = meta)
        names(assays(sce)) = counts_type
      }
    }
    else {
      sce = SingleCellExperiment(assay = counts)
      names(assays(sce)) = counts_type
    }

    if (counts_type == "counts" & transform_counts_to_logcounts){
      sce = .log_normalise(sce, batch,...)
    }
    return(sce)
  }
}


#' @import SingleCellExperiment
#' @importFrom scuttle logNormCounts
.log_normalise = function(sce, batch, ...){
  if (is.null(batch)){
    sce = .get_size_factors_single_batch(sce,...)
    sce = logNormCounts(sce)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    sce = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.sce = .get_size_factors_single_batch(sce[,idx],...)
      return(current.sce)
    })
    sce = do.call(multiBatchNorm , sce )
    sce = do.call(cbind, sce)
  }
  return(sce)
}

#' @importFrom scran quickCluster computeSumFactors
.get_size_factors_single_batch = function(sce, d=30, min.mean=0.1){
  clusters = quickCluster(sce, method="igraph", use.ranks=TRUE, d=d, min.mean=min.mean)
  sce = computeSumFactors(sce, clusters=clusters)
  return(sce)
}




