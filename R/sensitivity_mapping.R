
#' Returns celltype corresponding vector representing percentage of cells from the celltype that is mapped correctly
#'
#' @param mapping A \code{data.frame} object representing mapping.
#'
#' @return A \code{data.frame} object with fields 'celltype' and 'frac_mapped_correctly'.
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
#' mapping = out$mapping
#' sensitivity_stat = get_sensitivity_mapping(mapping)
get_sensitivity_mapping = function(mapping){
  if (!.valid_mapping_df(mapping)) {
    stop()
  } else {
    tab = table(mapping$celltype , mapping$celltype_mapped)
    tab = sweep(tab, 1, rowSums(tab), "/" )
    tab = as.data.frame(tab)
    colnames(tab) = c("celltype" , "celltype_mapped" , "frac")

    celltypes = as.character(unique(tab$celltype))
    stat = lapply(celltypes, function(celltype){
      current.tab = tab[tab$celltype == celltype , ]
      if (celltype %in% current.tab$celltype_mapped){
        out = data.frame(celltype = celltype , frac_mapped_correctly = current.tab$frac[current.tab$celltype_mapped == celltype])
      }
      else {
        out = data.frame(celltype = celltype , frac_mapped_correctly = 0)
      }
      return(out)
    })
    stat = do.call(rbind, stat)
    return(stat)
  }
}


#' @importFrom stats as.dendrogram
#' @importFrom ape as.phylo extract.clade
#' @import dendextend
.sensitivity_of_hierarchical_mapping = function(mapping , ct_hierarchy, sensitivity.thresh = 0.75){
  if (!is(ct_hierarchy, "hclust")) {
    stop("ct_hierarchy variable should be hclust object")
  } else{
    dend = as.dendrogram(ct_hierarchy)
    tree = as.phylo(dend)

    edges = tree$edge
    nodes = sort(unique(as.vector(edges)))
    tip.nodes = 1:length(tree$tip.label)
    tip.labels = tree$tip.label
    sensitivity_score = sapply(nodes , function(node){
      if (node %in% tip.nodes){
        celltypes = tip.labels[node]
        out = .get_frac_correct_within_split(mapping , celltypes)
        out = round(out,2)
        return(out)
      }
      else {
        current.tree = extract.clade(tree , node)
        celltypes = current.tree$tip.label
        out = .get_frac_correct_within_split(mapping , celltypes)
        out = round(out,2)
        return(out)
      }
    })
    celltypes.poorly_mapped = .get_celltypes_poorly_mapped(tree , sensitivity_score , sensitivity.thresh)
    out = list(tree = tree , sensitivity_score = sensitivity_score, celltypes.poorly_mapped = celltypes.poorly_mapped, sensitivity.thresh = sensitivity.thresh)
    return(out)
  }
}

#
.get_frac_correct_within_split = function(mapping , celltypes){
  frac_mapped_correct = sapply(celltypes , function(celltype){
    out = sum(mapping$celltype_mapped %in% celltypes & mapping$celltype == celltype)/sum(mapping$celltype == celltype)
    return(out)
  })
  return(min(frac_mapped_correct))
}


#
#' @importFrom phangorn Ancestors
.get_celltypes_poorly_mapped = function(tree , sensitivity_score , sensitivity.thresh){
  edges = tree$edge
  nodes = sort(unique(as.vector(edges)))
  tip.nodes = 1:length(tree$tip.label)
  tip.labels = tree$tip.label

  nodes.poorly_mapped = nodes[sensitivity_score < sensitivity.thresh]
  tip.nodes.poorly_mapped = intersect(tip.nodes,nodes.poorly_mapped)

  celltypes.poorly_mapped = lapply(tip.nodes.poorly_mapped, function(tip.node){
    celltypes_add_markers = tip.labels[tip.node]
    current.ancestors = Ancestors(tree, tip.node, type = "all")
    current.ancestors.well_mapped = current.ancestors[!current.ancestors %in% nodes.poorly_mapped]
    closest.ancestor.well_mapped = max(current.ancestors.well_mapped)
    current.tree = extract.clade(tree , closest.ancestor.well_mapped)
    celltypes_distinguish = current.tree$tip.label
    out = list(celltypes_add_markers = celltypes_add_markers , celltypes_distinguish = celltypes_distinguish)
    return(out)
  })
  return(celltypes.poorly_mapped)
}





