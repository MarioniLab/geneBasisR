

#' Produces binary tree representing celltype relationships (based on input genes)
#'
#' @param sce SingleCellExperiment object representing scRNA-seq counts matrix containing celltype field.
#' @param genes Genes that are used to assess celltypes similarity. Default=NULL, and in this case all rownames of sce would be assigned as genes.
#' @param batch A field from colData of sce object that indicates which field is used as batch. If specified, PCA would be calculated using multiBatchPCA and pcs would be batch-corrected. Default=NULL that assumes no batch.
#' @param cosineNorm A boolean specifying whether cosine normalization should be applied to logcounts. Default=TRUE.
#' @param nPC Number of PCs to be used. Default=NULL, and in this case all significant PCs (reflecting dependence between celltype and PC from ANOVA) will be included.
#' @param p.thresh A P-value threshold to determine significance of each PC. Default=0.05.
#'
#' @return hclust object representing relationship of celltypes.
#' @export
#' @importFrom stats dist hclust
#' @examples
#' require(SingleCellExperiment)
#' n_row = 30000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$celltype = as.factor(sample(1:5, n_col, replace=TRUE))
#' ct_hierarchy = get_ct_hierarchy(sce , genes = c(1:20), nPC = 10)
#'
get_ct_hierarchy = function(sce , genes = NULL , batch = NULL, cosineNorm = TRUE , nPC = NULL, p.thresh = 0.05, option = "hclust"){
  if (!.check_counts_matrix_correct(sce)) {
    stop()
  } else {
    if (is.null(genes)){
      genes = rownames(sce)
    }
    pcs_corrected = .get_corrected_pcs(sce , genes , batch , cosineNorm)

    # get only those w relevant pcs
    relevant_pcs = .get_relevant_for_celltypes_pcs(pcs_corrected , sce , nPC , p.thresh)
    pcs_corrected = pcs_corrected[, relevant_pcs]

    meta = as.data.frame(colData(sce))
    # get avg pc-value per ct
    celltypes = unique(sce$celltype)
    avg_pc_per_ct.stat = lapply(celltypes , function(celltype){
      cells = rownames(pcs_corrected) %in% meta$cell[meta$celltype == celltype]
      return( apply( pcs_corrected[cells,] , 2 , mean))
    })
    avg_pc_per_ct.stat = do.call(cbind , avg_pc_per_ct.stat)
    colnames(avg_pc_per_ct.stat) = celltypes

    # build hierarchy
    dist = dist(t(avg_pc_per_ct.stat) , diag=TRUE)
    ct_hierarchy = hclust(dist ,  method = "average")
    if (option == "list"){
      ct_hierarchy = .get_listed_ct_hierarchy(ct_hierarchy)
    }
    return(ct_hierarchy)
  }
}

# Transforms celltype hierarchy representation to the nested list format.
.get_listed_ct_hierarchy = function(ct_hierarchy){
  if (!is(ct_hierarchy , "hclust")){
    stop("input should be of hclust class")
  } else {
    ct_hierarchy.merge = ct_hierarchy$merge
    ct_hierarchy.labels = ct_hierarchy$labels

    out = list()
    for (i in 1:nrow(ct_hierarchy.merge)){
      if (ct_hierarchy.merge[i,1] < 0){
        ct_hierarchy.merge.1 = ct_hierarchy.labels[ -1*ct_hierarchy.merge[i,1] ]
      }
      else {
        ct_hierarchy.merge.1 = out[[ct_hierarchy.merge[i,1]]]
      }
      if (ct_hierarchy.merge[i,2] < 0){
        ct_hierarchy.merge.2 = ct_hierarchy.labels[ -1*ct_hierarchy.merge[i,2] ]
      }
      else {
        ct_hierarchy.merge.2 = out[[ct_hierarchy.merge[i,2]]]
      }
      out[[i]] = list(ct_hierarchy.merge.1 , ct_hierarchy.merge.2)
    }
    return(out[[nrow(ct_hierarchy.merge)]])
  }
}
