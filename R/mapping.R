#' MNN-based mapping of celltypes
#'
#' @param sce_reference SingleCellExperiment object representing reference scRNA-seq counts matrix containing celltype field.
#' @param sce_query SingleCellExperiment object representing query scRNA-seq counts matrix containing celltype field. Cells and accordingly celltype labels from this ds would be assigned neighbors from sce_reference.
#' @param genes Genes that are used for integration between sce_reference and sce_query
#' @param batch A field from colData of sce object that indicates which field is used as batch. If specified, PCA would be calculated using multiBatchPCA and pcs would be batch-corrected. Default=NULL that assumes no batch.
#' @param cosineNorm A boolean specifying whether cosine normalization should be applied to logcounts. Default=TRUE.
#' @param n.neigh Number of nearets neighbrs used for the integration. Default=10.
#' @param nPC Number of PCs to be used. Default=NULL, and in this case all significant PCs (reflecting dependence between celltype and PC from ANOVA) will be included.
#' @param p.thresh A P-value threshold to determine significance of each PC. Default=0.05.
#'
#' @return data.frame object, where every row corresponds to cell from sce_query. celltype_mapped field contains assigned from sce_reference celltype.
#' @export
#' @importFrom BiocNeighbors queryKNN
#' @examples
#'
mapping = function(sce_reference , sce_query , genes , batch = NULL , cosineNorm = TRUE , n.neigh = 10 , nPC = NULL , p.thresh = 0.05){
  if (!.check_counts_matrix_correct(sce_reference) | !.check_counts_matrix_correct(sce_query)) {
    stop()
  } else {
    if ( (sum(genes %in% rownames(sce_reference)) < length(genes)) | (sum(genes %in% rownames(sce_query)) < length(genes))){
      stop("Either reference or query ds do not have some genes in their rownames.")
    } else {
      sce_reference = sce_reference[rownames(sce_reference) %in% genes , ]
      sce_reference = sce_reference[order(rownames(sce_reference)),]
      sce_query = sce_query[rownames(sce_query) %in% genes , ]
      sce_query = sce_query[order(rownames(sce_query)),]
      sce = cbind(sce_reference , sce_query)
      meta = as.data.frame(colData(sce))

      # if only one celltype is present in the reference ds - assign it immediately.
      if (length(unique(sce_reference$celltype)) == 1){
        mapping = data.frame(cell = colnames(sce_query) , celltype_mapped = unique(sce_reference$celltype) )
        mapping$celltype_mapped = as.character(mapping$celltype_mapped)
      }
      else {
        pcs_corrected = .get_corrected_pcs(sce , genes , batch , cosineNorm)
        relevant_pcs = .get_relevant_for_celltypes_pcs(pcs_corrected , sce_reference, nPC, p.thresh)
        corrected_pcs = corrected_pcs[, relevant_pcs]

        reference_cells = colnames(sce_reference)
        query_cells = colnames(sce_query)
        knns = suppressWarnings( queryKNN( corrected_pcs[reference_cells ,], corrected_pcs[query_cells ,], k = n.neigh, get.index = TRUE) )
        cells_mapped = t( apply(knns$index, 1, function(x) reference_cells[x]) )
        celltypes = t(apply(cells_mapped, 1, function(x) meta$celltype[match(x, meta$cell)]))
        celltype_mapped = apply(celltypes, 1, function(x) .getmode(x, 1:length(x)))
        mapping = data.frame(cell = query_cells , celltype_mapped = celltype_mapped)
        mapping$celltype_mapped = as.character(mapping$celltype_mapped)
      }
      return(mapping)
    }
  }
}
