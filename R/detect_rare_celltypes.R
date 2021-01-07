#' Detecting and filtering rare celltypes from count matrix. Additionally, assigning markers for rare celltypes from the original count matrix.
#'
#' @param sce SingleCellExperiment object representing scRNA-seq counts matrix containing celltype field.
#' @param n_cells.thresh Threshold of how many cells should be used to detect rare celltypes. Default=10.
#' @param Top Positive integer representing how many Top markers should be returned.
#' @param filter_rare_celltypes Boolean representing whether filtered counts matrix should be returned.
#' @param get_markers_rare_celltypes Boolean representing whether markers for rare celltypes should be returned.
#'
#' @return Filtered counts matrix and/or markers for rare celltypes.
#' @export
#'
#' @examples
detect_rare_celltypes = function(sce , n_cells.thresh = 10 , Top = 5, filter_rare_celltypes = TRUE, get_markers_rare_celltypes = TRUE){
  if(!is.numeric(n_cells.thresh) | n_cells.thresh < 0){
    stop("n_cells.thresh must be a non-negative number.")
  }
  if (!.check_counts_matrix_correct(sce)) {
    stop()
  } else {
    celltypes_stat = table(sce$celltype)
    celltypes_2discard = names(celltypes_stat)[celltypes_stat < n_cells.thresh]
    if (isEmpty(celltypes_2discard)){
      message("No rare celltypes found.")
      out = sce
    } else {
      message(paste0("Detected next rare celltypes: " , celltypes_2discard))
      if (get_markers_rare_celltypes & filter_rare_celltypes){
        markers.rare_celltypes = .get_markers_rare_celltypes(sce, celltypes_2discard)
        sce.filtered = .filter_rare_celltypes(sce, celltypes_2discard)
        out = list(filtered_sce = sce.filtered, markers = markers.rare_celltypes)
      } else if (filter_rare_celltypes){
        sce.filtered = .filter_rare_celltypes(sce, celltypes_2discard)
        out = sce.filtered
      } else if (get_markers_rare_celltypes){
        markers.rare_celltypes = .get_markers_rare_celltypes(sce, celltypes_2discard)
        out = markers.rare_celltypes
      }
    }
    return(out)
  }
}

#' @importFrom scran findMarkers
.get_markers_rare_celltypes = function(sce, celltypes){
  markers <- findMarkers(sce , groups=sce$celltype, direction = "up", test = "t", assay.type = "logcounts", pval.type="any")
  markers.rare_celltypes = lapply(celltypes, function(celltype){
    current_sce = sce
    current_sce$celltype.bool = current_sce$celltype == celltype
    current_markers = as.data.frame( markers[[which(names(markers) == celltype)]] )
    return(current_markers[current_markers$Top <= Top ,])
  })
  names(markers.rare_celltypes) = celltypes
  return(markers.rare_celltypes)
}

.filter_rare_celltypes = function(sce, celltypes){
  sce.filtered = sce[ , !sce$celltype %in% celltypes]
  return(sce.filtered)
}




