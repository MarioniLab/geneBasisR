#' MouseEmbryo, E8.5
#' Counts matrix for mouse embryo, E8.5 (data extracted from Pijuan-Sala et al., 2019).
#'
#' SCE object with log-normalised counts matrix (pre-calculated with custom script) strored in assay logcounts.
#' colData contains information about sample (i.e. batch), celltype and provided color-code for celltypes.
#' reducedDim contains UMAP coordinates.
#' Counts matrix is already pre-processed and only non-mitochondrial genes with variance > 0 are retained (top 3000 to reduce memory usage and computational time).
#'
#' @docType data
#' @usage data(sce_mouseEmbryo)
#'
#' @format SCE object
#'
#' @name sce_mouseEmbryo
#' @source MouseGastrulationData package
NULL
