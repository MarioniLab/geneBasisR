#' Filters rare celltypes
#'
#' @param sce - a SingleCellExperiment object that contains celltype field
#' @param n_cells.thresh - a threshold of how many cells should be used to detect rare celltypes
#'
#' @return filtered sce
#' @export
#'
#' @examples
#'
get_rid_of_rare_celltypes = function(sce , n_cells.thresh = 10){

  if(!is.numeric(n_cells.thresh) | n_cells.thresh < 0){
    stop("n_cells.thresh must be a non-negative number.")
  }
  if (!.check_counts_matrix_correct(sce)) {
    stop()
  } else{
    celltypes_stat = table(sce$celltype)
    celltypes_2discard = names(celltypes_stat)[celltypes_stat < n_cells.thresh]
    sce = sce[ , !sce$celltype %in% celltypes_2discard]
    if (!isEmpty(celltypes_2discard)){
      message(paste0("Discraded these rare celltypes: " , celltypes_2discard))
    }
    else {
      message("No rare celltypes found")
    }
    return(sce)
  }
}
