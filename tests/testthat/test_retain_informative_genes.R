context("Testing retain_informative_genes")
library(geneBasisR)

### set up inputs

## libs
library(SingleCellExperiment)
set.seed(32)

## toy data
counts = t(matrix(c(0,0,0,0,1,0,1,1,1), nrow = 3, ncol = 3))
sce_1 = SingleCellExperiment(list(logcounts = counts))
rownames(sce_1) = c(1,2,3)

sce_2 = SingleCellExperiment(list(logcounts = counts))
rownames(sce_2) = c(1,1,2)

# load mouse ds
data("sce_mouseEmbryo", package = "geneBasisR")


test_that("Wrong input gives errors", {

  # should be unique rownames
  expect_error(retain_informative_genes(sce_2),
               "SCE should have unique rownames.",
               fixed=TRUE
  )
  # not good for scran's modelGeneVar
  expect_error(retain_informative_genes(sce_1),
               "Can not perform modelGeneVar on this counts matrix - check your input data.",
               fixed=TRUE
  )

  # n should be positive integer
  expect_error(retain_informative_genes(sce_mouseEmbryo, n = -3),
               "Check n - should be positive integer",
               fixed=TRUE
  )

  # n should be positive integer
  expect_error(retain_informative_genes(sce_mouseEmbryo, n = 3.5),
               "Check n - should be positive integer",
               fixed=TRUE
  )

  # n should be positive integer
  expect_error(retain_informative_genes(sce_mouseEmbryo, n = "test"),
               "Check n - should be positive integer",
               fixed=TRUE
  )

  # var.thresh should be a numeric
  expect_error(retain_informative_genes(sce_mouseEmbryo, var.thresh = "test"),
               "Check var.thresh - should be numeric",
               fixed=TRUE
  )

  # if n is too small, no inputes left to work w
  expect_error(retain_informative_genes(sce_mouseEmbryo, n = 1),
               "Less than 2 genes are selected. Consider checking your counts data, increasing n and/or decreasing var.thresh",
               fixed=TRUE
  )

  # if var.thresh is too big, no inputes left to work w
  expect_error(retain_informative_genes(sce_mouseEmbryo, n = 10, var.thresh = 1000000),
               "Less than 2 genes are selected. Consider checking your counts data, increasing n and/or decreasing var.thresh",
               fixed=TRUE
  )

  # if var.thresh is too big, no inputes left to work w
  expect_error(retain_informative_genes(sce_mouseEmbryo, var.thresh = 1000000),
               "Less than 2 genes are selected. Consider checking your counts data, increasing n and/or decreasing var.thresh",
               fixed=TRUE
  )
})


test_that("Large n will give same results", {
  sce_1 = retain_informative_genes(sce_mouseEmbryo, n = 20000)
  sce_2 = retain_informative_genes(sce_mouseEmbryo, n = 10000)
  expect_identical(sce_1, sce_2)
})


