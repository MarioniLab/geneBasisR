context("Testing get_DE_genes")
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


# load mouse ds
data("sce_mouseEmbryo", package = "geneBasisR")



test_that("Return is the correct class", {
  # correct class
  out = get_DE_genes(sce_correct, genes = rownames(sce_correct) , test.type = "t")
  expect_is(out, "data.frame")

  # correct colnames
  out = get_DE_genes(sce_correct, genes = rownames(sce_correct) , test.type = "t", FDR.thresh = 1)
  out = colnames(out)
  out_expect = c("celltype" , "gene" , "summary.logFC")
  expect_identical(out, out_expect)
})




test_that("Return of the correct output", {

  # if no genes are left - NULL
  out = get_DE_genes(sce_correct, genes = rownames(sce_correct) , test.type = "t", FDR.thresh = -.001)
  out_expect = NULL
  expect_identical(out, out_expect)

  # if FDR = 1, everything (at least where there are no NAs if there is no statistics on the test)
  out = get_DE_genes(sce_correct, genes = rownames(sce_correct) , test.type = "t", FDR.thresh = 1)
  out = nrow(out)
  out_expect = 10
  expect_equal(out, out_expect)

  # for toy data, expect correct gene to be selected
  out = get_DE_genes(sce_correct, genes = rownames(sce_correct) , test.type = "t", FDR.thresh = .01)
  out = out[, c("celltype" , "gene")]
  out_expect = data.frame(celltype = "2" , gene = "3")
  expect_identical(out, out_expect)

})


test_that("Wrong input, sce", {
  # should be unique rownames
  expect_error(get_DE_genes(sce_wrong_rownames),
               "SCE should have unique rownames.",
               fixed=TRUE
  )
  # sce should be sce
  expect_error(get_DE_genes(logcounts(sce_correct)),
               "SCE should be a SingleCellExperiment object.",
               fixed=TRUE
  )

  # sce shold contain field celltype
  expect_error(get_DE_genes(sce_no_celltype_field),
               "'celltype' field should be in colData(sce)",
               fixed=TRUE
  )

})



test_that("Wrong input, test.type", {
  # test.type should be: binom, wilcox or t
  expect_error(get_DE_genes(sce_correct, test.type = NULL),
               "Check test.type - should be either 'binom', 'wilcox' or 't'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_DE_genes(sce_correct, test.type = 0),
               "Check test.type - should be either 'binom', 'wilcox' or 't'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_DE_genes(sce_correct, test.type = "random"),
               "Check test.type - should be either 'binom', 'wilcox' or 't'",
               fixed=TRUE
  )
})




test_that("Wrong input, pval.type", {
  # test.type should be: all, some or any
  expect_error(get_DE_genes(sce_correct, pval.type = NULL),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_DE_genes(sce_correct, pval.type = 0),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_DE_genes(sce_correct, pval.type = "random"),
               "Check pval.type - should be either 'all', 'some' or 'any'",
               fixed=TRUE
  )
})




test_that("Wrong input, FDR.thresh", {
  # FDR.thresh should be numeric
  expect_error(get_DE_genes(sce_correct, FDR.thresh = NULL),
               "Check FDR.thresh - should be numeric",
               fixed=TRUE
  )
  # genes should be a subset of rownames in sce_correct
  expect_error(get_DE_genes(sce_correct, FDR.thresh = "random"),
               "Check FDR.thresh - should be numeric",
               fixed=TRUE
  )
})


