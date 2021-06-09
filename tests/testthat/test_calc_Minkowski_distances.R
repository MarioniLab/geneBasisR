context("Testing calc_Minkowski_distances")
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

# load mouse ds
data("sce_mouseEmbryo", package = "geneBasisR")



test_that("Return of the correct output", {
  out = calc_Minkowski_distances(sce_correct, genes = rownames(sce_correct) , n.neigh = 2, p = 1)
  out_expect = data.frame(gene = as.character(c(1:5)) , dist = c(3,0,2,0,0))
  expect_identical(out, out_expect)
})

test_that("Return of the correct output, run 2", {
  out = calc_Minkowski_distances(sce_correct, genes = as.character(c(1,2)) , n.neigh = 2, p = 1)
  out_expect = data.frame(gene = as.character(c(1:5)) , dist = c(0,3,5,0,0))
  expect_identical(out, out_expect)
})

test_that("Return of the correct output, run 3", {
  out = calc_Minkowski_distances(sce_correct, genes = as.character(c(1,2)) , n.neigh = 2, p = 1, genes.predict = as.character(c(3,5)))
  out_expect = data.frame(gene = as.character(c(3,5)) , dist = c(5,0))
  expect_identical(out, out_expect)
})


test_that("Wrong input gives errors", {
  # should be unique rownames
  expect_error(calc_Minkowski_distances(sce_wrong_rownames, genes = rownames(sce_wrong_rownames)),
               "SCE should have unique rownames.",
               fixed=TRUE
  )

  # sce should be sce
  expect_error(calc_Minkowski_distances(logcounts(sce_correct), genes = rownames(sce_correct)),
               "SCE should be a SingleCellExperiment object.",
               fixed=TRUE
  )

  # genes should be character
  expect_error(calc_Minkowski_distances(sce_correct, genes = c(1,2,3)),
               "Check genes - should be NULL or character vector",
               fixed=TRUE
  )

  # n.neigh - positive scalar > 1
  expect_error(calc_Minkowski_distances(sce_correct, genes = rownames(sce_correct), n.neigh = -1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )

  # n.neigh - positive scalar > 1
  expect_error(calc_Minkowski_distances(sce_correct, genes = rownames(sce_correct), n.neigh = 1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )

  # n.neigh - positive scalar > 1
  expect_error(calc_Minkowski_distances(sce_correct, genes = rownames(sce_correct), n.neigh = "all"),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )


  # batch should be the field in sce
  expect_error(calc_Minkowski_distances(sce_correct_w_batch, genes = rownames(sce_correct), n.neigh = 2, batch = "sample"),
               "Batch should be one the colnames in colData(sce).",
               fixed=TRUE
  )

  # n.neigh should be < min(size(batch)) - 1
  expect_error(calc_Minkowski_distances(sce_correct_w_batch, genes = rownames(sce_correct), n.neigh = 3, batch = "batch"),
               "Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.",
               fixed=TRUE
  )



})

