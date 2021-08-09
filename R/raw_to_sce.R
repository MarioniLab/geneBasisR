


#' raw_to_sce
#'
#' Essentially a parser: raw counts (or log-normalized counts) stored in a file (e.g. .txt) -> SingleCellExperiment object of the right format.
#' If raw counts are passed, log-normalization is performed (optional but recommended) and logcounts will be used downstream.
#' If no batch is supplied, size factors are calculated using scran::quickCluster (arguments d and min.mean regulate the coarsity of clustering) and then caluclated
#' size factors are supplied to scuttle::logNormCounts to calulcate log-normalized counts.
#' If batch is supplied, size factors are calculated per within batch and then size factors are supplied to batchelor::multiBatchNorm.
#' @param counts_dir String specifying the directory for counts matrix (assuming counts where already calculated)
#' @param counts_type String specifying whether raw data is stored as counts or log-counts (= 'counts' and 'logcounts' respectively).
#' For geneBasis we recommend to work with log-counts. If you do not have log-counts precomputed, they can be computed within this function.
#' @param transform_counts_to_logcounts In case, raw data are counts (as opposed to log-counts), Boolean specifying whether we should perform log-normalization.
#' @param header Boolean specifying if logcounts_dir file has cell IDs stored in colnames.
#' @param sep the field separator string. Note it should be the same for logcounts_dir and meta_dir (if latter exists).
#' @param meta_dir If not NULL (NULL is default), a string specifying the directory for meta-data (i.e. celltype, batch, UMAP-coordinates).
#' Store UMAP-coordinates as 'x' and 'y' (relevant for plotting functions).
#' Also, if meta contains field cell - this field will be used for cell IDs (so ensure the values are unique).
#' @param batch If not NULL (i.e. no batch, NULL is default), string specifying a column in meta file that will be used as batchID. Please check that specified batch name exists in meta-file.
#' @param verbose Boolean identifying whether intermediate print outputs should be returned. Default verbose=TRUE.
#' @param d Only used for log-normalization: an integer scalar specifying the number of principal components to retain.
#' @param min.mean Only used for log-normalization: a numeric scalar specifying the filter to be applied on the average count for each filter prior to computing ranks.
#' @param ... Additional arguments. This includes d and min.mean for scran::quickCluster - used to calculate size factors to compute normalized log-counts.
#'
#' @return SingleCellExperiment object with gene counts/logcounts and meta-data (if supplied) stored in colData.
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment colData
#' @importFrom SummarizedExperiment assays
#' @importFrom scuttle readSparseCounts
#' @examples
#' require(SingleCellExperiment)
#' counts_dir = system.file("extdata", "raw_spleen.txt", package = "geneBasisR")
#' meta_dir = system.file("extdata", "raw_spleen_meta.txt", package = "geneBasisR")
#' out = raw_to_sce(counts_dir, counts_type = "logcounts", transform_counts_to_logcounts = FALSE, header = TRUE, sep = "\t" , meta_dir = meta_dir, batch = NULL)
#'
raw_to_sce = function(counts_dir, counts_type = "counts", transform_counts_to_logcounts = TRUE, header = TRUE, sep = "\t" ,
                      meta_dir = NULL, batch = NULL, verbose = TRUE, d = 50, min.mean = 0.1, ...){
  args = c(as.list(environment()), list(...))
  out = .general_check_arguments(args)

  if (!file.exists(counts_dir)){
    stop("Counts file does not exist")
  }
  else {
    if (verbose){
      cat("Counts matrix is being processed.\n")
    }
    counts = readSparseCounts(counts_dir, sep = sep, row.names = TRUE, col.names = header)
    if (!header){
      colnames(counts) = c(1:ncol(sce))
    }
    if (is.null(meta_dir)){
      sce = SingleCellExperiment(list(counts = counts))
      names(assays(sce)) = counts_type
    }
    else {
      if(!file.exists(meta_dir) ){
        stop("Meta file is supplied, but does not exist")
      }
      else {
        if (verbose){
          cat("Meta file is being processed.\n")
        }
        meta = read.table(meta_dir, header = TRUE, sep = sep)
        if (!"cell" %in% colnames(meta)){
          meta$cell = c(1:nrow(meta))
        }
        if (!is.null(batch)){
          if (!batch %in% colnames(meta)){
            stop("Batch should be one the colnames in meta.")
          }
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
    }
    if (counts_type == "counts" & transform_counts_to_logcounts){
      sce = .log_normalise(sce, batch = batch, verbose = verbose, d = d, min.mean = min.mean)
    }
    return(sce)
  }
}


#' @importFrom SingleCellExperiment colData
#' @importFrom scuttle logNormCounts
#' @importFrom batchelor multiBatchNorm
.log_normalise = function(sce, batch = NULL, verbose = TRUE, d = 50, min.mean = 0.1){
  if (verbose){
    cat("Calculating size factors (to perform log-normalization).\n")
  }
  if (is.null(batch)){
    res = tryCatch(
      {
        sce_corrected = .get_size_factors_single_batch(sce, d = d, min.mean = min.mean)
        sce_corrected = logNormCounts(sce)
      },
      error = function(dump){
        message("Can not calculate size factors. Either calculate log-counts separately or tune d / min.mean.")
        return(sce)
      }
    )
    return(res)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    res = tryCatch(
      {
        sce_corrected = lapply(unique(batchFactor) , function(current.batch){
          if (verbose){
            cat(paste0("Calculating size factors for batch " , current.batch , ".\n"))
          }
          idx = which(batchFactor == current.batch)
          current.sce = .get_size_factors_single_batch(sce[,idx], d = d, min.mean = min.mean)
          return(current.sce)
        })
        sce_corrected = do.call(multiBatchNorm , sce_corrected )
        sce_corrected = do.call(cbind, sce_corrected)
        return(sce_corrected)
      },
      error = function(dump){
        message("Can not calculate size factors for at least one batch. Either calculate log-counts separately, tune d / min.mean or consider log normalization for all data together (no batch specified, set batch = NULL).")
        return(sce)
      }
    )
    return(res)
  }
}

#' @importFrom scran quickCluster computeSumFactors
.get_size_factors_single_batch = function(sce, d=50, min.mean=0.1){
  clusters = quickCluster(sce, method="igraph", use.ranks=TRUE, d=d, min.mean=min.mean)
  sce = computeSumFactors(sce, clusters=clusters)
  return(sce)
}


