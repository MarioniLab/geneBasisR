
#' MNN-based mapping in hierarchical fashion
#'
#' @param sce_reference SingleCellExperiment object representing reference scRNA-seq counts matrix containing celltype field.
#' @param sce_query SingleCellExperiment object representing query scRNA-seq counts matrix containing celltype field. Cells and accordingly celltype labels from this ds would be assigned neighbors from sce_reference.
#' @param genes Genes that are used for integration between sce_reference and sce_query.
#' @param batch A field from colData of sce object that indicates which field is used as batch. If specified, PCA would be calculated using multiBatchPCA and pcs would be batch-corrected. Default=NULL that assumes no batch.
#' @param cosineNorm A boolean specifying whether cosine normalization should be applied to logcounts. Default=TRUE.
#' @param n.neigh Positive scalar specifying number of nearest neighbors to be used in the integration. Default=10.
#' @param nPC Number of PCs to be used for the integration between sce_reference and sce_query. Default=NULL, and in this case all significant PCs (reflecting dependence between celltype and PC from ANOVA) will be included.
#' @param nPC.ct_hierarchy Number of PCs to be used in building celltype hierarchy. Default=NULL, and in this case all significant PCs (reflecting dependence between celltype and PC from ANOVA) will be included.
#' @param p.thresh A P-value threshold to determine significance of each PC (for the integration between sce_reference and sce_query). Default=0.05.
#' @param p.thresh.ct_hierarchy A P-value threshold to determine significance of each PC (for building the celltype hierarchy). Default=0.05.
#'
#' @return mapping = data.frame object, where every row corresponds to cell from sce_query. 'celltype' field contains actual celltype labels. 'celltype_mapped' field contains assigned from sce_reference celltype. mapping_sensitivity = list containing objects: phylo object 'tree' - alternative representation of ct_hierarchy; vector 'sensitivity_score' - a sensitivity of mapping, each element corresponds to node (i.e. clade or set of celltypes) in 'tree'; list 'celltypes.poorly_mapped' - each element represents clade with poor mapping sensitivity such as no ancestor clade has poor mapping sensitivity.
#' @export
#'
#' @examples
#' require(SingleCellExperiment)
#' n_row = 30000
#' n_col = 100
#' cells_reference = c(1:80)
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$celltype = as.factor(sample(1:5, n_col, replace=TRUE))
#' genes = c(1:20)
#' sce_reference = sce[ , colnames(sce) %in% cells_reference]
#' sce_query = sce[ , !colnames(sce) %in% cells_reference]
#' out = hierarchical_mapping(sce_reference , sce_query , genes, nPC = 5, nPC.ct_hierarchy=5)
#'
hierarchical_mapping = function( sce_reference , sce_query , genes , batch = NULL ,
                                 cosineNorm = TRUE , n.neigh = 10 , nPC = NULL , nPC.ct_hierarchy = NULL,
                                 p.thresh = 0.05, p.thresh.ct_hierarchy = 0.05,
                                 get.sensitivity = TRUE, sensitivity.thresh = 0.75){
  if (!.check_counts_matrix_correct(sce_reference) | !.check_counts_matrix_correct(sce_query)) {
    stop()
  } else {
    if ( (sum(genes %in% rownames(sce_reference)) < length(genes)) | (sum(genes %in% rownames(sce_query)) < length(genes)) ){
      stop("Either reference or query ds do not have some genes in their rownames.")
    } else {
      ct_hierarchy = get_ct_hierarchy(sce_reference , genes = genes , batch = batch,
                                      cosineNorm = cosineNorm , nPC = nPC.ct_hierarchy,
                                      p.thresh = p.thresh.ct_hierarchy, option = "both")

      # select genes entries and order
      sce_reference = sce_reference[rownames(sce_reference) %in% genes , ]
      sce_reference = sce_reference[order(rownames(sce_reference)),]
      sce_query = sce_query[rownames(sce_query) %in% genes , ]
      sce_query = sce_query[order(rownames(sce_query)),]

      mapping = .hierarchical_mapping_intermediate(sce_reference, sce_query, ct_hierarchy$list, batch, cosineNorm, n.neigh, nPC, p.thresh)
      meta = as.data.frame(colData(sce_query))
      mapping = merge(mapping , meta[, c("cell" , "celltype")])

      if (get.sensitivity){
        mapping_sensitivity = .sensitivity_of_hierarchical_mapping(mapping , ct_hierarchy$hclust, sensitivity.thresh)
        out = list(mapping = mapping , mapping_sensitivity = mapping_sensitivity)
      } else {
        out = list(mapping = mapping)
      }
      return(out)
    }
  }
}


.hierarchical_mapping_intermediate = function( sce_reference , sce_query, ct_hierarchy , batch = NULL , cosineNorm = TRUE , n.neigh = 10 , nPC = NULL , p.thresh = 0.05){
  mapping = .mapping_single_split(sce_reference, sce_query, ct_hierarchy, batch, cosineNorm, n.neigh, nPC, p.thresh)
  mappings = lapply(1:2, function(i){
    if (!is(ct_hierarchy[[i]] , "list")){
      current_mapping = mapping[mapping$celltype_mapped == ct_hierarchy[[i]] , ]
      return(current_mapping)
    } else {
      current_sce_reference = sce_reference[, sce_reference$celltype %in% unlist( ct_hierarchy[[i]] ) ]
      current_sce_query = sce_query[, colnames(sce_query) %in% mapping$cell[mapping$celltype_mapped == i]]
      current_ct_hierarchy = ct_hierarchy[[i]]
      if (ncol(current_sce_query) > 0){
        return(.hierarchical_mapping_intermediate(current_sce_reference , current_sce_query, current_ct_hierarchy, batch, cosineNorm, n.neigh, nPC, p.thresh))
      }
    }
  })
  mappings = do.call(rbind , mappings)
  return(mappings)
}


.mapping_single_split = function(sce_reference, sce_query, ct_hierarchy , batch = NULL , cosineNorm = TRUE , n.neigh = 10, nPC = NULL , p.thresh){
  if (!.check_counts_matrix_correct(sce_reference) | !.check_counts_matrix_correct(sce_query)) {
    stop()
  } else {
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
    mapping = .mapping_intermediate(sce_reference, sce_query, batch = batch, cosineNorm = cosineNorm, n.neigh = n.neigh, nPC = nPC, p.thresh = p.thresh)
    if (!is.list(ct_hierarchy[[1]])){
      mapping$celltype_mapped[mapping$celltype_mapped == 1] = celltypes.1
    }
    if (!is.list(ct_hierarchy[[2]])){
      mapping$celltype_mapped[mapping$celltype_mapped == 2] = celltypes.2
    }
    return(mapping)
  }
}


