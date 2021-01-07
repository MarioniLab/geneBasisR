
#' Returns ordered gene selection (by their relevance to PC-space)
#'
#' @param sce SingleCellExperiment object representing scRNA-seq counts matrix
#' @param genes Initially pre-selected DE (per celltype) genes that we assign as relevant for generation of PC-space
#' @param type Setting to assign whether covariance matrix in PCA should be reweighted so each celltype has the same weight (celltype). Otherwise setting should be set to none and each cell will be accounted once.
#'
#' @return list containing stat per PC and order gene vector
#' @export
#' @import forcats
#' @examples
#' require(SingleCellExperiment)
#' n_row = 30000
#' n_col = 100
#' sce = SingleCellExperiment(assays = list(logcounts = matrix(rnorm(n_row*n_col), ncol=n_col)))
#' rownames(sce) = as.factor(1:n_row)
#' colnames(sce) = c(1:n_col)
#' sce$cell = colnames(sce)
#' sce$celltype = as.factor(sample(1:5, n_col, replace=TRUE))
#' genes_ranked = get_gene_ranking(sce , genes = c(1:20))

get_gene_ranking = function(sce , genes, type = "none"){

  if (!(type %in% c("celltype" , "none"))){
    stop("type should be assigned to either celltype or none - choose one.")
  }
  loadings = .get_loadings(sce , genes , type)
  stat = lapply(1:ncol(loadings), function(nPC){
    current_loadings = loadings[, nPC]
    out.max = data.frame(nPC = nPC , sign = "pos" , gene = names(which.max(current_loadings)) , loading = max(current_loadings) )
    out.min = data.frame(nPC = nPC , sign = "neg" , gene = names(which.min(current_loadings)) , loading = min(current_loadings))
    out = rbind(out.max , out.min)
  })
  stat = do.call(rbind , stat)
  stat$gene = factor(stat$gene)
  genes = levels(fct_inorder(stat$gene))
  out = list(stat = stat , genes = genes)
  return(out)
}


#' @importFrom stats prcomp
#' @importFrom batchelor multiBatchPCA
#'
.get_loadings = function(sce , genes , type){
  set.seed(32)
  if (!.check_counts_matrix_correct(sce)) {
    stop()
  } else{
    if (!(type %in% c("celltype" , "none"))){
      stop("type should be assigned to either celltype or none - choose one.")
    }
    sce = sce[rownames(sce) %in% as.character(genes) , ]
    counts = as.matrix(logcounts(sce))
    if (type == "celltype"){
      mbpca = multiBatchPCA(counts , batch = factor( sce$celltype ) , get.variance = T , d = nrow(counts))
      loadings = metadata(mbpca)
      mbpca = do.call(rbind , mbpca)
      out = loadings$rotation
    }
    else if (type == "none") {
      pca = prcomp(t(counts) , rank. = nrow(counts))
      out = pca$rotation
    }
    return(out)
  }
}

