context("Testing gene_search")
library(geneBasisR)

### set up inputs

## libs
library(SingleCellExperiment)
set.seed(32)

## toy data
counts = t(matrix(c(0,0,0,0,1,0,1,0,1,1,1,1), nrow = 3, ncol = 4))
sce_1 = SingleCellExperiment(list(logcounts = counts))
rownames(sce_1) = c(1,2,3,4)

sce_2 = SingleCellExperiment(list(logcounts = counts))
rownames(sce_2) = c(1,1,2,3)

stat_all_1 = data.frame(gene = c(1:4) , dist_all = c(0,1,1,0))
stat_all_2 = data.frame(gene = c(1:4) , dist = c(0,1,1,0))

# load mouse ds
data("sce_mouseEmbryo", package = "geneBasisR")

test_that("Wrong input gives errors", {
  # should be unique rownames
  expect_error(gene_search(sce_2),
               "SCE should have unique rownames.",
               fixed=TRUE
  )
  # stat_all should be the right format
  expect_error(gene_search(sce_1, n_genes_total = 1, stat_all = stat_all_2, n.neigh = 2),
               "stat_all is of the wrong format - should contain fields gene and dist_all.",
               fixed=TRUE
  )

  # n_genes_total should be positive integer
  expect_error(gene_search(sce_mouseEmbryo, n_genes_total = -3),
               "Check n_genes_total - should be positive integer",
               fixed=TRUE
  )

  # n_genes_total should be positive integer
  expect_error(gene_search(sce_mouseEmbryo, n_genes_total = 3.5),
               "Check n_genes_total - should be positive integer",
               fixed=TRUE
  )

  # n_genes_total should be positive integer
  expect_error(gene_search(sce_mouseEmbryo, n_genes_total = "test"),
               "Check n_genes_total - should be positive integer",
               fixed=TRUE
  )

  # n_genes_total should be bigger than nrow(sce)
  expect_error(gene_search(sce_1, n_genes_total = 4, n.neigh = 2),
               "Selected library size is bigger than number of genes in the counts matrix",
               fixed=TRUE
  )

})


test_that("Return of the correct output", {
  out = gene_search(sce_1, stat_all = stat_all_1 , n_genes_total = 1, n.neigh = 2, nPC.all = NULL)
  expect_identical(out$gene, "2")
})


