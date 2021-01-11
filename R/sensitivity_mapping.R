
#' Returns celltype corresponding vector representing percentage of cells from the celltype that is mapped correctly
#'
#' @param mapping A \code{data.frame} object representing mapping.
#'
#' @return A \code{data.frame} object with fields 'celltype' and 'frac_mapped_correctly'.
#' @export
#'
#' @examples
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
    celltypes.poorly_mapped = .get_clades_poorly_mapped(tree , sensitivity_score , sensitivity.thresh)
    out = list(tree = tree , sensitivity_score = sensitivity_score, celltypes.poorly_mapped = celltypes.poorly_mapped)
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
.get_clades_poorly_mapped = function(tree , sensitivity_score , sensitivity.thresh){
  edges = tree$edge
  nodes = sort(unique(as.vector(edges)))
  tip.nodes = 1:length(tree$tip.label)
  tip.labels = tree$tip.label

  nodes.poorly_mapped = nodes[sensitivity_score < sensitivity.thresh]
  ancestors.poorly_mapped = sapply(nodes.poorly_mapped , function(node){
    current.ancestors = Ancestors(tree, node, type = "all")
    if (sum(nodes.poorly_mapped %in% current.ancestors) == 0){
      return(F)
    } else {
      return(T)
    }
  })
  nodes.poorly_mapped = nodes.poorly_mapped[ancestors.poorly_mapped == F]
  celltypes.poorly_mapped = lapply(nodes.poorly_mapped , function(node){
    if (node %in% tip.nodes){
      out = tip.labels[node]
    } else {
      current.tree = extract.clade(tree , node)
      out = current.tree$tip.label
    }
    return(out)
  })
  return(celltypes.poorly_mapped)
}





