# Contains various check functions to examine whether variables are of the right format


#' @import SingleCellExperiment
.prepare_sce_counts = function(sce){
  if (!is(sce , "SingleCellExperiment")){
    stop("SCE should be a SingleCellExperiment object.")
    return(F)
  } else if ("logcounts" %in% names(assays(sce))){
    return(sce)
  } else if (!"counts" %in% names(assays(sce))){
    stop("SCE should contain either counts or logcounts.")
    return(F)
  } else {
    message("Logcounts assay is not found. Counts will be used instead.")
    logcounts(sce) = counts(sce)
    return(sce)
  }
}

#' @import SingleCellExperiment
.prepare_sce_cell_ids = function(sce){
  if (!is(sce , "SingleCellExperiment")){
    stop("SCE should be a SingleCellExperiment object.")
    return(F)
  }
  else {
    meta = colData(sce)
    if (!"cell" %in% colnames(meta)){
      if (is.null(colnames(sce))){
        colnames(sce) = c(1:ncol(sce))
      }
      #message("Using SCE colnames as cell IDs")
      sce$cell = colnames(sce)
      meta$cell = colnames(sce)
    }
    if (length(unique(meta$cell)) < nrow(meta)){
      stop("Cell entries should contain unique IDs.")
      return(F)
    }
    else {
      return(sce)
    }
  }
}

.prepare_sce = function(sce){
  sce = .prepare_sce_counts(sce)
  sce = .prepare_sce_cell_ids(sce)
  rownames(sce) = as.character(rownames(sce))
  if (!is(sce , "SingleCellExperiment")){
    stop("SCE variable did not suffice must have properties.")
    return(F)
  }
  else {
    return(sce)
  }
}


#' @import SingleCellExperiment
.check_sce = function(sce){
  if (!is(sce , "SingleCellExperiment")){
    stop("SCE should be a SingleCellExperiment object.")
    return(F)
  } else if (!("logcounts" %in% names(assays(sce)))){
    stop("SCE should contain logcounts assay.")
    return(F)
  } else if (length(unique(rownames(sce))) < nrow(sce) ){
    stop("SCE should have unique rownames.")
    return(F)
  } else {
    return(T)
  }
}


.check_celltype_in_sce = function(sce){
  if (.check_sce(sce)){
    if (!("celltype" %in% colnames(colData(sce)))){
      stop("'celltype' field should be in colData(sce)")
      return(F)
    }
    else {
      return(T)
    }
  }
}


.check_genes_in_sce = function(sce, genes){
  if (.check_sce(sce)){
    if (!is.null(genes)){
      out = mean(genes %in% rownames(sce))
      if (out < 1){
        stop("Some gene names are missing from SCE.")
        return(F)
      }
      else {
        return(T)
      }
    }
    else {
      return(T)
    }
  }
}

.check_batch = function(sce, batch){
  if (.check_sce(sce)){
    if (!is.null(batch)){
      meta = as.data.frame(colData(sce))
      if (!batch %in% colnames(meta)){
        stop("Batch should be one the colnames in colData(sce).")
        return(F)
      }
      else {
        return(T)
      }
    }
    else {
      return(T)
    }
  }
}

.check_celltype_mapping_correct = function(mapping){
  if (!class(mapping) == data.frame){
    stop("Celltype mapping should be stored as data.frame")
    return(F)
  } else if (sum(c("celltype" , "mapped_celltype") %in% colnames(mapping)) < 2){
    stop("Either 'celltype' and/or 'mapped_celltype' are not in colnames(mapping)")
    return(F)
  } else {
    return(T)
  }
}

.check_argument_correct = function(dots, arg_name , fun , message){
  if (arg_name %in% names(dots)){
    arg = dots[[which(names(dots) == arg_name)]]
    out = fun(arg)
    if (!out){
      stop(message)
    }
    return(out)
  }
  else {
    return(TRUE)
  }
}


.general_check_arguments = function(dots){
  out = TRUE
  out = .check_argument_correct(dots, "sce", .check_sce, "Check sce - something is wrong (gene names unique? logcounts assay exists?)")
  out = .check_argument_correct(dots, "n.neigh", .check_n.neigh, "Check n.neigh - should be positive integer > 1")
  out = .check_argument_correct(dots, "nPC", .check_positive_integer_or_null, "Check nPC - should be NULL or positive integer")
  out = .check_argument_correct(dots, "nPC.all", .check_positive_integer_or_null, "Check nPC.all - should be NULL or positive integer")
  out = .check_argument_correct(dots, "nPC.selection", .check_positive_integer_or_null, "Check nPC.selection - should be NULL or positive integer")
  out = .check_argument_correct(dots, "batch", .check_string_or_null, "Check batch - should be NULL or string")
  out = .check_argument_correct(dots, "genes", .check_string_or_null, "Check genes - should be NULL or character vector")
  out = .check_argument_correct(dots, "genes_base", .check_string_or_null, "Check genes_base - should be NULL or character vector")
  out = .check_argument_correct(dots, "genes_to_assess", is.character, "Check genes_to_assess - should be NULL or character vector")
  out = .check_argument_correct(dots, "genes.selection", is.character, "Check genes.selection - should be character vector")
  out = .check_argument_correct(dots, "genes.all", is.character, "Check genes.all - should be character vector")
  out = .check_argument_correct(dots, "genes.predict", is.character, "Check genes.predict - should be character vector")
  out = .check_argument_correct(dots, "cosine", .check_boolean, "Check cosine - should be boolean")
  out = .check_argument_correct(dots, "verbose", .check_boolean, "Check verbose - should be boolean")
  out = .check_argument_correct(dots, "get.dist", .check_boolean, "Check get.dist - should be boolean")
  out = .check_argument_correct(dots, "discard.mt", .check_boolean, "Check discard.mt - should be boolean")
  out = .check_argument_correct(dots, "select.hvgs", .check_boolean, "Check select.hvgs - should be boolean")
  out = .check_argument_correct(dots, "return.cell_score_stat", .check_boolean, "Check return.cell_score_stat - should be boolean")
  out = .check_argument_correct(dots, "return.gene_score_stat", .check_boolean, "Check return.gene_score_stat - should be boolean")
  out = .check_argument_correct(dots, "return.celltype_stat", .check_boolean, "Check return.celltype_stat - should be boolean")
  out = .check_argument_correct(dots, "return.stat", .check_boolean, "Check return.stat - should be boolean")
  out = .check_argument_correct(dots, "p.minkowski", .check_positive_integer, "Check p.minkowski - should be positive integer")
  out = .check_argument_correct(dots, "n_genes_total", .check_positive_integer, "Check n_genes_total - should be positive integer")
  out = .check_argument_correct(dots, "n_genes.step", .check_positive_integer, "Check n_genes.step - should be positive integer")
  out = .check_argument_correct(dots, "n", .check_positive_integer_or_null, "Check n - should be positive integer")
  out = .check_argument_correct(dots, "FDR.thresh", is.numeric, "Check FDR.thresh - should be numeric")
  out = .check_argument_correct(dots, "var.thresh", is.numeric, "Check var.thresh - should be numeric")
  out = .check_argument_correct(dots, "corr_all.thresh", is.numeric, "Check corr_all.thresh - should be numeric")
  out = .check_argument_correct(dots, "library.size_type", function(x) .check_arg_within_options(x, c("single", "series")),
                                "Check library.size_type - should be either 'single' or 'series'")
  out = .check_argument_correct(dots, "method", function(x) .check_arg_within_options(x, c("spearman", "pearson")),
                                "Check method - should be either 'spearman' or 'pearson'")
  out = .check_argument_correct(dots, "test.type", function(x) .check_arg_within_options(x, c("binom", "wilcox", "t")),
                                "Check test.type - should be either 'binom', 'wilcox' or 't'")
  out = .check_argument_correct(dots, "pval.type", function(x) .check_arg_within_options(x, c("all", "some", "any")),
                                "Check pval.type - should be either 'all', 'some' or 'any'")
  out = .check_argument_correct(dots, "which_genes_to_use", function(x) .check_arg_within_options(x, c("all", "DE")),
                                "Check which_genes_to_use - should be either 'all' or 'DE'")
  return(out)
}


.check_arg_within_options = function(x , options){
  out = TRUE
  if (is.null(x)){
    out = FALSE
  }
  else  if (!x %in% options){
    out = FALSE
  }
  return(out)
}


.check_positive_integer = function(x){
  out = TRUE
  if (!is.numeric(x)){
    out = FALSE
  } else if (!x%%1 == 0 | x <= 0){
    out = FALSE
  }
  return(out)
}

.check_n.neigh = function(x){
  out = TRUE
  if (!is.numeric(x)){
    out = FALSE
  } else if (!x%%1 == 0 | x <= 1){
    out = FALSE
  }
  return(out)
}


.check_integer_or_all = function(x){
  out = TRUE
  if (!x == "all"){
    if (!is.numeric(x)){
      out = FALSE
    } else if (!x%%1 == 0){
      out = FALSE
    }
  }
  return(out)
}

.check_integer_or_null = function(x){
  out = TRUE
  if (!is.null(x)){
    if (!is.numeric(x)){
      out = FALSE
    } else if (!x%%1 == 0){
      out = FALSE
    }
  }
  return(out)
}

.check_positive_integer_or_null = function(x){
  out = TRUE
  if (!is.null(x)){
    if (!is.numeric(x)){
      out = FALSE
    } else if (!x%%1 == 0 | x <= 0){
      out = FALSE
    }
  }
  return(out)
}

.check_string_or_null = function(x){
  out = TRUE
  if (!is.null(x)){
    if (!is.character(x)){
      out = FALSE
    }
  }
  return(out)
}


.check_boolean = function(x){
  out = TRUE
  if (!out %in% c(TRUE, FALSE)){
    out = FALSE
  }
  return(out)
}
