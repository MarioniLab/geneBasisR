# context("Testing gene_prediction_scores")
library(geneBasisR)

### set up inputs

## libs
library(SingleCellExperiment)
set.seed(32)

## toy data
counts = matrix(c(rep(0,10),1,1,1,0,0, 1,1,1,0,0, 1,1, 0,0,0,1,0,0,0,0,0,1,1,0,0), nrow = 5, ncol = 7)
sce_correct = SingleCellExperiment(list(logcounts = counts))
colnames(sce_correct) = c(1:7)
rownames(sce_correct) = c(1:5)

sce_correct_w_batch = SingleCellExperiment(list(logcounts = counts))
colnames(sce_correct_w_batch) = c(1:7)
rownames(sce_correct_w_batch) = c(1:5)
sce_correct_w_batch$batch = c(1,1,1,"all","all",1,"all")

sce_wrong_rownames = SingleCellExperiment(list(logcounts = counts))
colnames(sce_wrong_rownames) = c(1:7)
rownames(sce_wrong_rownames) = c(1,1,2,3,"a")

counts = matrix(c(rep(0,35)) , nrow = 5, ncol = 7)
sce_zeros = SingleCellExperiment(list(logcounts = counts))
colnames(sce_zeros) = c(1:7)
rownames(sce_zeros) = c(1:5)

# load mouse ds
data("sce_mouseEmbryo", package = "geneBasisR")



test_that("Return is the correct class", {
  # right class
  out = get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2)
  expect_s3_class(out, "data.frame")

  # right colnames exist
  out = get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2)
  out = colnames(out)
  out_expect = c("cell" , "cell_score")
  out = sum(out_expect %in% out)
  expect_equal(out, 2)
})



test_that("Return of the correct output, simple", {
  # all genes, 0s
  out = get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2)
  out = out[, c("cell" , "cell_score")]
  out_expect = data.frame(cell = as.character(c(1:7)) , cell_score = rep(1,7))
  expect_equal(out$cell_score, out_expect$cell_score)

  # all genes, but matrix is 0L --> NaNs
  out = get_neighborhood_preservation_scores(sce_zeros, genes.selection = rownames(sce_zeros) , n.neigh = 2)
  out = out[, c("cell" , "cell_score")]
  out_expect = data.frame(cell = as.character(c(1:7)) , cell_score = rep(NaN,7))
  expect_equal(out$cell_score, out_expect$cell_score)

})

test_that("neighs.all_stat matches number of batches", {
  #out = get_neighs_all_stat(sce_mouseEmbryo, batch = "sample")
  #names_out = names(out)
  #expect_equal(length(out), 3)
  #expect_equal(names_out, c("counts", "neighs.all" , "mean_dist"))

  out = get_neighs_all_stat(sce_mouseEmbryo, batch = NULL)
  names_out = names(out)
  expect_equal(length(out), 3)
  expect_equal(names_out, c("counts", "neighs.all" , "mean_dist"))

})


test_that("Wrong input, sce", {
  # should be unique rownames
  expect_error(get_neighborhood_preservation_scores(sce_wrong_rownames, genes.selection = rownames(sce_correct)),
               "sce should have unique rownames.",
               fixed=TRUE
  )
  # sce should be sce
  expect_error(get_neighborhood_preservation_scores(logcounts(sce_correct), genes.selection = rownames(sce_correct)),
               "sce should be a SingleCellExperiment object.",
               fixed=TRUE
  )
})


test_that("Wrong input, genes.selection", {
  # genes.selection should be character
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = c(1,2,3)),
               "Check genes.selection - should be character vector",
               fixed=TRUE
  )
  # genes.selection should be a subset of rownames in sce_correct
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(1,2,3,6))),
               "Some gene names are missing from SCE.",
               fixed=TRUE
  )
  # genes.selection is not NULL
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = NULL),
               "Check genes.selection - should be character vector",
               fixed=TRUE
  )
  # is ok if genes.selections are ok
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2),
               NA)
})


test_that("Wrong input, genes.all", {
  # genes.selection should be character
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(1:5)), genes.all = c(1,4,5)),
               "Check genes.all - should be character vector",
               fixed=TRUE
  )
  # genes.selection should be a subset of rownames in sce_correct
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(1:5)), genes.all = as.character(c(1,2,3,6))),
               "Some gene names are missing from SCE.",
               fixed=TRUE
  )
  # genes.selection is not NULL
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(1:5)) , genes.all = NULL),
               "Check genes.all - should be character vector",
               fixed=TRUE
  )
  # is ok if genes.selections are ok
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(1:5)), genes.all = as.character(c(1,4,5))),
               NA)

  # is ok if genes.selections are ok
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(1:5)), corr_all.thresh = 1, genes.all = as.character(c(1,4,5))),
               NA)
  # single gene in predict
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(1:5)), corr_all.thresh = 0,
                                          genes.predict ="1", genes.all = as.character(c(1:5))),
               NA)

})



test_that("Wrong input, batch", {
  # batch is NULL or character
  expect_error(get_neighborhood_preservation_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = 1),
               "Check batch - should be NULL or string",
               fixed=TRUE
  )
  # batch should be the field in sce
  expect_error(get_neighborhood_preservation_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = "sample"),
               "Batch should be one the colnames in colData(sce).",
               fixed=TRUE
  )
  # no errors if batch correct
  expect_error(get_neighborhood_preservation_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = "batch"),
               NA
  )
})




test_that("Wrong input, n.neigh", {
  # n.neigh - positive scalar > 1
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = -1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh - positive scalar > 1
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh - positive scalar > 1
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = "random"),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # internal check - n.neigh can be 'all' but not for users
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = "all"),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh should be < min(size(batch)) - 1
  # expect_error(get_neighborhood_preservation_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 3, batch = "batch"),
  #              "Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.",
  #              fixed=TRUE
  # )
})


test_that("Wrong input, nPC.all", {
  # nPC.all - NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 0),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = -10),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 5.5),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = "all"),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC -  NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = NULL),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 3),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 30000),
               NA)
})


test_that("Wrong input, nPC.selection", {
  # nPC.all - NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = 0),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = -10),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = 5.5),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = "all"),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC -  NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = NULL),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = 3),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 3 , nPC.selection = 30000),
               NA)
})



test_that("neighs.all_stat is of correct format", {

  # neighs.all can be NULL
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2),
               NA)
  # neighs.all can be NULL
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, neighs.all_stat = NULL),
               NA
  )
  # neighs.all_stat should be list w right colnames
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = "2", n.neigh = 2,
                                                    neighs.all_stat = data.frame(counts = 1 , neighs.all = 1 , mean_dist = 1)),
               "neighs.all_stat should be a list. Precompute neighs.all_stat with geneBasisR::get_neighs_all_stat.",
               fixed=TRUE
  )
  # neighs.all_stat should be list w right colnames
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2,
                                                    neighs.all_stat = list(counts = 1 , neigh.all = 1 , mean_dist = 1)),
               "neighs.all_stat should contain fields 'counts', 'neighs.all' and 'mean_dist'. Precompute neighs.all_stat with geneBasisR::get_neighs_all_stat.",
               fixed=TRUE
  )

  # neighs.all_stat can be calculated w no errors
  neighs.all_stat = get_neighs_all_stat(sce_correct , batch = NULL, option = "approx")
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2,
                                          neighs.all_stat = neighs.all_stat),
               NA
  )
  # neighs.all_stat can be calculated w no errors
  neighs.all_stat = get_neighs_all_stat(sce_correct_w_batch , batch = "batch", option = "approx", n.neigh = 2)
  expect_error(get_neighborhood_preservation_scores(sce_correct_w_batch, genes.selection = as.character(c(2)), n.neigh = 2, batch = "batch",
                                          neighs.all_stat = neighs.all_stat),
               NA
  )

  # neighs.all_stat can be calculated w no errors
  neighs.all_stat = get_neighs_all_stat(sce_correct , batch = NULL, option = "exact")
  expect_error(get_neighborhood_preservation_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2,
                                                    neighs.all_stat = neighs.all_stat),
               NA
  )
  # neighs.all_stat can be calculated w no errors
  neighs.all_stat = get_neighs_all_stat(sce_correct_w_batch , batch = "batch", option = "exact", n.neigh = 2)
  expect_error(get_neighborhood_preservation_scores(sce_correct_w_batch, genes.selection = as.character(c(2)), n.neigh = 2, batch = "batch",
                                                    neighs.all_stat = neighs.all_stat),
               NA
  )
})




