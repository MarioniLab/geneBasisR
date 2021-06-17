


#' raw_to_sce
#'
#' Simply a parser raw .txt files --> SingleCellExperiment object + adding logcounts
#'
#' @param counts_dir String specifying the directory for counts matrix (assuming counts where already calculated)
#' @param counts_type String specifying whether raw data is stored as counts or log-counts (= 'counts' and 'logcounts' respectively).
#' For geneBasis we recommend to work with log-counts. If you do not have log-counts precomputed, they can be computed within this function.
#' @param transform_counts_to_logcounts In case, raw data are counts (as opposed to log-counts), Boolean specifying whether we should perform log-normalization.
#' @param header Boolean specifying if logcounts_dir file has cell IDs stored in colnames.
#' @param sep the field separator string. Note it should be the same for logcounts_dir and meta_dir (if latter exists).
#' @param meta_dir If not NULL, a string specifying the directory for meta-data (i.e. celltype, batch, UMAP-coordinates).
#' Store UMAP-coordinates as 'x' and 'y' (relevant for plotting functions).
#' Also, if meta contains field cell - this field will be used for cell IDs (so ensure the values are unique).
#' @param batch If not NULL (i.e. no batch), string specifying a column in meta file that will be used as batchID. Please check that specified batch name exists in meta-file.
#' @param ... Additional arguments. This includes d and min.mean for scran::quickCluster - used to calculate size factors to compute normalized log-counts.
#'
#' @return SingleCellExperiment object with gene counts/logcounts and meta-data (if supplied) stored in colData.
#' @export
#' @import SingleCellExperiment
#' @examples
#' require(SingleCellExperiment)
#' counts_dir = system.file("extdata", "raw_spleen.txt", package = "geneBasisR")
#' meta_dir = system.file("extdata", "raw_spleen_meta.txt", package = "geneBasisR")
#' out = raw_to_sce(counts_dir, counts_type = "logcounts", transform_counts_to_logcounts = FALSE, header = TRUE, sep = "\t" , meta_dir = meta_dir, batch = NULL)
#'
raw_to_sce = function(counts_dir, counts_type = "counts", transform_counts_to_logcounts = TRUE, header = TRUE, sep = "\t" , meta_dir = NULL, batch = NULL,...){
  if (!file.exists(counts_dir)){
    stop("Counts file does not exist")
  }
  else if (!is.null(meta_dir) & !file.exists(meta_dir)) {
    stop("Meta file is supplied, but does not exist")
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
      } else if (mean(meta$cell == colnames(counts)) < 1){
        stop("Mismatch in cell IDs between counts matrix and meta-file.")
      } else {
        sce = SingleCellExperiment(list(counts = counts), colData = meta)
        names(assays(sce)) = counts_type
      }
    }
    else {
      sce = SingleCellExperiment(list(counts = counts))
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




