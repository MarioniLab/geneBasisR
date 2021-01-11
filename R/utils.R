## These are utility functions that are called internally in various final functions of the package


# The function does a check whether SingleCellExperiment objest is of the correct format. To be inserted to almost every function prior to everything else.
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
.check_counts_matrix_correct = function(sce){
  if (!is(sce , "SingleCellExperiment")){
    stop("sce should be SingleCellExperiment object")
    return(F)
  } else if (!("logcounts" %in% names(assays(sce)))){
    stop("sce should contain assay slot logcounts")
    return(F)
  } else if (!("celltype" %in% colnames(colData(sce)))){
    stop("colData of sce should contain celltype slot")
    return(F)
  } else {
    return(T)
  }
}

# The function does a check for correct mapping data.frame format
.valid_mapping_df = function(mapping){
  if (!is(mapping , "data.frame")){
    stop("mapping data.frame should be data.frame object")
    return(F)
  } else if (!"celltype" %in% colnames(mapping) | !"celltype_mapped" %in% colnames(mapping)){
    stop("mapping data.frame should contain columns 'celltype' and 'celltype_mapped'.")
    return(F)
  } else {
    return(T)
  }
}

# The function returns corrected PCs on the selected genes
#' @importFrom batchelor cosineNorm multiBatchPCA reducedMNN
#' @import SingleCellExperiment
.get_corrected_pcs = function(sce , genes , batch = NULL, cosineNorm = TRUE){
  set.seed(32)
  if (!.check_counts_matrix_correct(sce)) {
    stop()
  } else {
    # assign whether we work with cosine normalized logcounts
    if (cosineNorm){
      counts = cosineNorm( logcounts(sce[rownames(sce) %in% genes , ]) )
    } else {
      counts = logcounts(sce[rownames(sce) %in% genes , ])
    }
    # perform pca and if multiple batches exist - perform correction
    if (!is.null(batch)){
      meta = colData(sce)
      if (!(batch %in% colnames(meta))){
        stop("batch specified variable should exist in SingleCellExperiment object - are you sure there is no typo?")
      } else {
        batchFactor = factor(meta[, colnames(meta) == batch])
        pcs = suppressWarnings( multiBatchPCA(counts, batch = batchFactor, d = length(genes)) )
        if (length(unique(batchFactor)) > 1){
          pcs_corrected = suppressWarnings( do.call(reducedMNN, as.list(pcs)) )
          pcs_corrected = pcs_corrected$corrected
        } else {
          pcs_corrected = pcs[[1]]
        }
      }
    } else {
      batchFactor = factor(rep(1,1,ncol(counts)))
      pcs_corrected = suppressWarnings( multiBatchPCA(counts, batch = batchFactor, d = length(genes)))
      pcs_corrected = pcs_corrected[[1]]
    }
    return(pcs_corrected)
  }
}


# The function returns PCs that show the most significant dependance between PCs and celltypes
#' @importFrom tibble rownames_to_column
#' @importFrom stats aov p.adjust
#' @importFrom stringr str_remove
.get_relevant_for_celltypes_pcs = function(pcs_corrected , sce_reference , nPC = NULL, p.thresh = 0.05){
  if (!.check_counts_matrix_correct(sce_reference)) {
    stop()
  } else {
    pcs_corrected = as.data.frame( pcs_corrected[rownames(pcs_corrected) %in% colnames(sce_reference) ,] )
    pc_cols = paste0("pc_", c(1:ncol(pcs_corrected)))
    colnames(pcs_corrected) = pc_cols

    pcs_corrected = rownames_to_column(pcs_corrected , var = "cell")
    meta = as.data.frame(colData(sce_reference))
    pcs_corrected = merge(pcs_corrected , meta)

    anova_stat = lapply(pc_cols , function(pc_col){
      fit = aov(pcs_corrected[, pc_col] ~ pcs_corrected$celltype)
      out = data.frame(pc = str_remove(pc_col , "pc_") , p = summary(fit)[[1]][["Pr(>F)"]][1])
      return(out)
    })
    anova_stat = do.call(rbind , anova_stat)
    anova_stat = anova_stat[order(anova_stat$p , decreasing = F) , ]
    anova_stat$p = p.adjust(anova_stat$p, method = "BH", n = length(anova_stat$p))

    if ( !is.null(nPC) ){
      out = as.numeric(  anova_stat$pc[1:min(nPC, nrow(anova_stat))] )
      return(out)
    } else {
      out = as.numeric( anova_stat$pc[anova_stat$p < p.thresh] )
      if (isEmpty(out)){
        stop("Number of components is not specified and none of the components have p-value less than threshold. Either specify number of components or increase p-value threshold.")
      } else {
        return(out)
      }
    }
  }
}


# returns most popular element in vector. If tied, returns the one, which appears first.
.getmode = function(v, dist) {
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}
