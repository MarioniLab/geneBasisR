# context("Testing get_DE_genes")
library(geneBasisR)

### set up inputs

## libs
library(SingleCellExperiment)
set.seed(32)

## toy data
counts = matrix(c(rep(0,10),1,1,1,0,0, 1,1,1,0,0, 1,1,0,0,0,1,0,0,0,0,0,1,1,0,0), nrow = 5, ncol = 7)

sce_correct = SingleCellExperiment(list(logcounts = counts))
colnames(sce_correct) = c(1:7)
rownames(sce_correct) = c(1:5)
sce_correct$celltype = c(1,1,"2","2","1","1","2")

sce_wrong_rownames = SingleCellExperiment(list(logcounts = counts))
colnames(sce_wrong_rownames) = c(1:7)
rownames(sce_wrong_rownames) = c(1,1,2,3,"a")
sce_wrong_rownames$celltype = c(1,1,"2","2","1","1","2")

sce_no_celltype_field = SingleCellExperiment(list(logcounts = counts))
colnames(sce_no_celltype_field) = c(1:7)
rownames(sce_no_celltype_field) = c(1:5)
sce_no_celltype_field$Celltype = c(1,1,"2","2",1,1,"2")


sce_correct_w_batch = SingleCellExperiment(list(logcounts = counts))
colnames(sce_correct_w_batch) = c(1:7)
rownames(sce_correct_w_batch) = c(1:5)
sce_correct_w_batch$celltype = c(1,1,"2","2","1","1","2")
sce_correct_w_batch$batch = c(4,4,4,1,1,1,1)


test_that("Return is the correct class", {

  # correct class, return.stat = T
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, return.stat = T)
  expect_type(out, "list")

  # correct class, return.stat = F
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, return.stat = F)
  expect_type(out, "list")


  # correct length of list, return.stat = T
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, return.stat = T)
  out = length(out)
  expect_equal(out, 2)


  # correct length of list, return.stat = F
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, return.stat = F)
  out = length(out)
  expect_equal(out, 1)


  #correct colnames
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, return.stat = T)
  out = out$mapping
  out = colnames(out)
  out_expect = c("cell" , "celltype" , "mapped_celltype")
  expect_identical(out, out_expect)

  #correct colnames
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, return.stat = T)
  out = out$stat
  out = colnames(out)
  out_expect = c("celltype" , "frac_correctly_mapped")
  expect_identical(out, out_expect)

  #correct colnames
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, return.stat = F)
  out = out$mapping
  out = colnames(out)
  out_expect = c("cell" , "celltype" , "mapped_celltype")
  expect_identical(out, out_expect)
})



test_that("Change in FDR does not affect the option 'all' but does affect the option 'DE'", {
  expect_error( get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, return.stat = F, FDR = .1) , NA)
  expect_error( get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2,
                                     return.stat = F, FDR = .1 , which_genes_to_use = "DE") ,
                "No DE genes discovered with these settings - consider option 'all' or tuning test, type and/or FDR threshold.",
                fixed = TRUE)
  expect_error( get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, test.type = "t" , return.stat = F, FDR = .1) , NA)

})



test_that("Return of the correct output", {

  # generally works
  expect_error(get_celltype_mapping(sce_correct_w_batch, genes.selection = rownames(sce_correct), n.neigh = 2, batch = "batch"),
               NA
  )

  # n.neigh = 2
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), nPC.selection = NULL, n.neigh = 2)
  out = out$mapping$mapped_celltype
  out_expect = as.character(c(1,1,2,2,2,1,2))
  expect_equal(out, out_expect)

  # n.neigh = 5
  out = get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), nPC.selection = NULL, n.neigh = 5)
  out = out$mapping$mapped_celltype
  out_expect = as.character(rep(1,7))
  expect_equal(out, out_expect)


})


test_that("Wrong input, sce", {
  # should be unique rownames
  expect_error(get_celltype_mapping(sce_wrong_rownames, genes.selection = rownames(sce_wrong_rownames)),
               "SCE should have unique rownames.",
               fixed=TRUE
  )
  # sce should be sce
  expect_error(get_celltype_mapping(logcounts(sce_correct), genes.selection = rownames(sce_correct)),
               "SCE should be a SingleCellExperiment object.",
               fixed=TRUE
  )

  # sce should contain field celltype
  expect_error(get_celltype_mapping(sce_no_celltype_field, genes.selection = rownames(sce_no_celltype_field)),
               "'celltype' field should be in colData(sce)",
               fixed=TRUE
  )

})




test_that("Wrong input, genes.selection", {
  # genes should be character
  expect_error(get_celltype_mapping(sce_correct, genes.selection = c(1:5)),
               "Check genes.selection - should be character vector",
               fixed=TRUE
  )
  # genes should be character
  expect_error(get_celltype_mapping(sce_correct, genes.selection = NULL),
               "Check genes.selection - should be character vector",
               fixed=TRUE
  )
  # genes.predict should be a subset of rownames in sce_correct
  expect_error(get_celltype_mapping(sce_correct, genes.selection = as.character(c(1,2,3,6))),
               "Some gene names are missing from SCE.",
               fixed=TRUE
  )

})



test_that("Wrong input, batch", {
  # batch is NULL or character
  expect_error(get_celltype_mapping(sce_correct_w_batch, genes = rownames(sce_correct_w_batch), n.neigh = 2, batch = 1),
               "Check batch - should be NULL or string",
               fixed=TRUE
  )
  # batch should be the field in sce
  expect_error(get_celltype_mapping(sce_correct_w_batch, genes = rownames(sce_correct_w_batch), n.neigh = 2, batch = "sample"),
               "Batch should be one the colnames in colData(sce).",
               fixed=TRUE
  )
})



test_that("Wrong input, n.neigh", {
  # n.neigh - positive scalar > 1
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), n.neigh = -1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh - positive scalar > 1
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh - positive scalar > 1
  expect_error(calc_Minkowski_distances(sce_correct, genes = rownames(sce_correct), n.neigh = "random"),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )

  # internal check - n.neigh can be 'all' but not for users
  expect_error(calc_Minkowski_distances(sce_correct, genes = rownames(sce_correct), n.neigh = "all"),
               NA
  )
  # n.neigh should be < ncol(sce)
  expect_error(get_celltype_mapping(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), nPC.selection = NULL, n.neigh = 7, batch = NULL),
               "n.neigh should be less than number of cells. Decrease n.neigh.",
               fixed=TRUE
  )

  expect_error(get_celltype_mapping(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), nPC.selection = NULL, n.neigh = 6, batch = NULL),
               NA
  )


})




test_that("Wrong input, pval.type", {
  # test.type should be: all, some or any
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), FDR.thresh = 1, pval.type = NULL),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), FDR.thresh = 1, pval.type = 0),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), FDR.thresh = 1, pval.type = "random"),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), FDR.thresh = 1, pval.type = "any"),
               NA
  )


  # test.type should be: all, some or any
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), which_genes_to_use = "DE" , FDR.thresh = 1, pval.type = NULL),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), which_genes_to_use = "DE" , FDR.thresh = 1, pval.type = 0),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), which_genes_to_use = "DE" , FDR.thresh = 1, pval.type = "random"),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_celltype_mapping(sce_correct, genes.selection = rownames(sce_correct), which_genes_to_use = "DE" , FDR.thresh = 1, pval.type = "any"),
               NA
  )

})




# test_that("Wrong input, FDR.thresh", {
#   # FDR.thresh should be numeric
#   expect_error(get_DE_genes(sce_correct, FDR.thresh = NULL),
#                "Check FDR.thresh - should be numeric",
#                fixed=TRUE
#   )
#   # genes should be a subset of rownames in sce_correct
#   expect_error(get_DE_genes(sce_correct, FDR.thresh = "random"),
#                "Check FDR.thresh - should be numeric",
#                fixed=TRUE
#   )
# })
#
#
