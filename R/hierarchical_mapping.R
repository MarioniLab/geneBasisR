
#' MNN-based mapping in hierarchical fashion.
#'
#' @param sce_reference SingleCellExperiment object representing reference scRNA-seq counts matrix containing celltype field.
#' @param sce_query SingleCellExperiment object representing query scRNA-seq counts matrix containing celltype field. Cells and accordingly celltype labels from this ds would be assigned neighbors from sce_reference.
#' @param genes Genes that are used for integration between sce_reference and sce_query
#' @param ct_hierarchy Nested lists object representing celltype hierarchy
#' @param batch A field from colData of sce object that indicates which field is used as batch. If specified, PCA would be calculated using multiBatchPCA and pcs would be batch-corrected. Default=NULL that assumes no batch.
#' @param cosineNorm A boolean specifying whether cosine normalization should be applied to logcounts. Default=TRUE.
#' @param n.neigh Number of nearets neighbrs used for the integration. Default=10.
#' @param nPC Number of PCs to be used. Default=NULL, and in this case all significant PCs (reflecting dependence between celltype and PC from ANOVA) will be included.
#' @param p.thresh A P-value threshold to determine significance of each PC. Default=0.05.

#' @return data.frame object, where every row corresponds to cell from sce_query. celltype_mapped field contains assigned from sce_reference celltype.
#' @export
#'
#' @examples
get_hierarchical_mapping = function( sce_reference , sce_query , genes, ct_hierarchy , batch = NULL , cosineNorm = TRUE , n.neigh = 10 , nPC = NULL , p.thresh = 0.05){
  mapping = .get_mapping_single_split(sce_reference, sce_query, genes, ct_hierarchy, batch, cosineNorm, n.neigh, nPC, p.thresh)
  mappings = lapply(1:2, function(i){
    if (!is(ct_hierarchy[[i]] , "list")){
      current_mapping = mapping[mapping$celltype_mapped == ct_hierarchy[[i]] , ]
      return(current_mapping)
    }
    else {
      current_sce_reference = sce_reference[, sce_reference$celltype %in% unlist( ct_hierarchy[[i]] ) ]
      current_sce_query = sce_query[, colnames(sce_query) %in% mapping$cell[mapping$celltype_mapped == i]]
      current_ct_hierarchy = ct_hierarchy[[i]]
      return(get_hierarchical_mapping(current_sce_reference , current_sce_query , genes , current_ct_hierarchy, batch, cosineNorm, n.neigh, nPC, p.thresh))
    }
  })
  mappings = do.call(rbind , mappings)
  return(mappings)
}


.get_mapping_single_split = function(sce_reference , sce_query , genes , ct_hierarchy , batch = NULL , cosineNorm = TRUE , n.neigh = 10, nPC = NULL , p.thresh){
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

      sce_reference = sce_reference[, sce_reference$celltype %in% unlist(ct_hierarchy)]
      celltypes.1 = unlist(ct_hierarchy[[1]])
      celltypes.2 = unlist(ct_hierarchy[[2]])
      sce_reference$celltype = sapply(sce_reference$celltype , function(celltype) ifelse(celltype %in% celltypes.1 , 1 , 2) )
      sce_query$celltype = sapply(sce_query$celltype , function(celltype){
        if (celltype %in% celltypes.1){
          return(1)
        }
        else if (celltype %in% celltypes.2){
          return(2)
        } else {
          return(0)
        }
      })
      mapping = mapping(sce_reference, sce_query, genes, batch, cosineNorm, n.neigh, nPC, p.thresh)
      if (!is.list(ct_hierarchy[[1]])){
        mapping$celltype_mapped[mapping$celltype_mapped == 1] = celltypes.1
      }
      if (!is.list(ct_hierarchy[[2]])){
        mapping$celltype_mapped[mapping$celltype_mapped == 2] = celltypes.2
      }
      return(mapping)
    }
  }
}


