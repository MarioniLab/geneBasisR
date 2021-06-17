# context("Testing evaluate_library")
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

sce_correct_w_celltype = sce_correct
sce_correct_w_celltype$celltype = c(1,3,2,1,2,2,3)

sce_correct_w_batch = SingleCellExperiment(list(logcounts = counts))
colnames(sce_correct_w_batch) = c(1:7)
rownames(sce_correct_w_batch) = c(1:5)
sce_correct_w_batch$batch = c(1,1,1,"all","all",1,"all")
sce_correct_w_batch$celltype = c(1,3,2,1,2,2,3)

sce_wrong_rownames = SingleCellExperiment(list(logcounts = counts))
colnames(sce_wrong_rownames) = c(1:7)
rownames(sce_wrong_rownames) = c(1,1,2,3,"a")

counts = matrix(c(rep(0,35)) , nrow = 5, ncol = 7)
sce_zeros = SingleCellExperiment(list(logcounts = counts))
colnames(sce_zeros) = c(1:7)
rownames(sce_zeros) = c(1:5)
sce_zeros$celltype = rep(1,7)

# load mouse ds
data("sce_mouseEmbryo", package = "geneBasisR")



test_that("Return is the correct class", {
  # right class
  out = evaluate_library(sce_correct, genes.selection = rownames(sce_correct) , return.celltype_stat = F)
  expect_type(out, "list")

  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) , return.celltype_stat = F)
  expect_type(out, "list")

  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) , return.celltype_stat = T)
  expect_type(out, "list")

  # right colnames exist
  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) , n.neigh = 2)
  out = names(out)
  out_expect = c("cell_score_stat" , "gene_score_stat" , "celltype_stat")
  expect_equal(out, out_expect)
})


test_that("Correct names for correct assays to return", {
  #
  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) ,
                         return.celltype_stat = F , return.gene_score_stat = F , return.cell_score_stat = T)
  out = names(out)
  out_expect = c("cell_score_stat")
  expect_equal(out, out_expect)
  #
  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) ,
                         return.celltype_stat = F , return.gene_score_stat = T , return.cell_score_stat = F, corr_all.thresh = 0)
  out = names(out)
  out_expect = c("gene_score_stat")
  expect_equal(out, out_expect)
  #
  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) ,
                         return.celltype_stat = T , return.gene_score_stat = F , return.cell_score_stat = F)
  out = names(out)
  out_expect = c("celltype_stat")
  expect_equal(out, out_expect)
  #
  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) ,
                         return.celltype_stat = F , return.gene_score_stat = T , return.cell_score_stat = T, corr_all.thresh = 0)
  out = names(out)
  out_expect = c("cell_score_stat", "gene_score_stat")
  expect_equal(out, out_expect)
  #
  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) ,
                         return.celltype_stat = T , return.gene_score_stat = F , return.cell_score_stat = T, corr_all.thresh = 0)
  out = names(out)
  out_expect = c("cell_score_stat", "celltype_stat")
  expect_equal(out, out_expect)
  #
  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) ,
                         return.celltype_stat = T , return.gene_score_stat = T , return.cell_score_stat = F, corr_all.thresh = 0)
  out = names(out)
  out_expect = c("gene_score_stat", "celltype_stat")
  expect_equal(out, out_expect)
  #
  out = evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype) ,
                         return.celltype_stat = T , return.gene_score_stat = T , return.cell_score_stat = T, corr_all.thresh = 0)
  out = names(out)
  out_expect = c("cell_score_stat", "gene_score_stat", "celltype_stat")
  expect_equal(out, out_expect)
  # error if none are to be returned
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype),
                                return.celltype_stat = F, return.gene_score_stat = F, return.cell_score_stat = F),
               "Select at least one stat to return.",
               fixed=TRUE
  )
})



test_that("Return of the correct output, simple", {
  # all genes, 0s
  out = evaluate_library(sce_zeros, genes.selection = rownames(sce_zeros) , n.neigh = 2)
  out_expect = list(cell_score_stat = data.frame(cell = as.character(c(1:7)) , cell_score = rep(NaN, 7) , n_genes = as.factor(5)) ,
                    celltype_stat = data.frame(celltype = "1", frac_correctly_mapped = 1, n_genes = as.factor(5)))
  expect_equal(out, out_expect)
})



test_that("Wrong input, sce", {
  # should be unique rownames
  expect_error(evaluate_library(sce_wrong_rownames, genes.selection = rownames(sce_wrong_rownames), return.celltype_stat = F),
               "sce should have unique rownames.",
               fixed=TRUE
  )
  # sce should be sce
  expect_error(evaluate_library(logcounts(sce_correct_w_celltype), genes.selection = rownames(sce_correct_w_celltype)),
               "sce should be a SingleCellExperiment object.",
               fixed=TRUE
  )
})


test_that("Wrong input, genes.selection", {
  # genes.selection should be character
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = c(1,2,3)),
               "Check genes.selection - should be character vector",
               fixed=TRUE
  )
  # genes.selection should be a subset of rownames in sce_correct
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = as.character(c(1,2,3,6))),
               "Some gene names are missing from SCE.",
               fixed=TRUE
  )
  # genes.selection is not NULL
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = NULL),
               "Check genes.selection - should be character vector",
               fixed=TRUE
  )
  # is ok if genes.selections are ok
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 2),
               NA)
})


test_that("Wrong input, genes.all", {
  # genes.selection should be character
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = as.character(c(1:5)), genes.all = c(1,4,5)),
               "Check genes.all - should be character vector",
               fixed=TRUE
  )
  # genes.selection should be a subset of rownames in sce_correct
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = as.character(c(1:5)), genes.all = as.character(c(1,2,3,6))),
               "Some gene names are missing from SCE.",
               fixed=TRUE
  )
  # genes.selection is not NULL
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = as.character(c(1:5)) , genes.all = NULL),
               "Check genes.all - should be character vector",
               fixed=TRUE
  )
  # is ok if genes.selections are ok
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = as.character(c(1:5)), genes.all = as.character(c(1,4,5))),
               NA)

  # is ok if genes.selections are ok
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = as.character(c(1:5)), corr_all.thresh = 1, genes.all = as.character(c(1,4,5))),
               NA)
  # single gene in predict
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = as.character(c(1:5)), corr_all.thresh = 0,
                                          genes.predict ="1", genes.all = as.character(c(1:5))),
               NA)
})



test_that("Wrong input, batch", {
  # batch is NULL or character
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = 1),
               "Check batch - should be NULL or string",
               fixed=TRUE
  )
  # batch should be the field in sce
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = "sample"),
               "Batch should be one the colnames in colData(sce).",
               fixed=TRUE
  )
  # no errors if batch correct
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = "batch"),
               NA
  )
})




test_that("Wrong input, n.neigh", {
  # n.neigh - positive scalar > 1
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = -1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh - positive scalar > 1
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh - positive scalar > 1
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = "random"),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # internal check - n.neigh can be 'all' but not for users
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = "all"),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh should be < min(size(batch)) - 1
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 7, batch = NULL),
               "Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.",
               fixed=TRUE
  )
  # n.neigh should be < min(size(batch)) - 1
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 3, batch = "batch"),
               "Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.",
               fixed=TRUE
  )
})


test_that("Wrong input, nPC.all", {
  # nPC.all - NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 2 , nPC.all = 0),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 2 , nPC.all = -10),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 2 , nPC.all = 5.5),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 2 , nPC.all = "all"),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC -  NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 2 , nPC.all = NULL),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 2 , nPC.all = 3),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_celltype, genes.selection = rownames(sce_correct_w_celltype), n.neigh = 2 , nPC.all = 30000),
               NA)
})


test_that("Wrong input, nPC.selection", {
  # nPC.all - NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2 , nPC.selection = 0),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2 , nPC.selection = -10),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2 , nPC.selection = 5.5),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2 , nPC.selection = "all"),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC -  NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2 , nPC.selection = NULL),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2 , nPC.selection = 3),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2 , nPC.all = 3 , nPC.selection = 30000),
               NA)
})



test_that("neighs.all is of correct format", {
  # neighs.all can be NULL
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = as.character(c(2)), n.neigh = 2),
               NA)
  # neighs.all can be NULL
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = as.character(c(2)), n.neigh = 2, neighs.all = NULL),
               NA
  )
  # neighs.all should be list w right colnames
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = "2", n.neigh = 2,
                                                    neighs.all = data.frame(cells_mapped = 1 , distances = 1)),
               "Something is wrong with neighs.all argument. For each batch, neighs.all should be a list containing 'cells_mapped' and 'distances' entries; nrow for each entry == n-cells in the batch. Consider recalculating using get_z_scaled_distances function.",
               fixed=TRUE
  )
  # neighs.all should be list w right colnames
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = as.character(c(2)), n.neigh = 2,
                                                    neighs.all = list(cell_mapped = 1 , distances = 1)),
               "Something is wrong with neighs.all argument. For each batch, neighs.all should be a list containing 'cells_mapped' and 'distances' entries; nrow for each entry == n-cells in the batch. Consider recalculating using get_z_scaled_distances function.",
               fixed=TRUE
  )
  # neighs.all should be list w right colnames
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = as.character(c(2)), n.neigh = 2,
                                                    neighs.all = list(cells_mapped = 1 , distances = 1)),
               "Something is wrong with neighs.all argument. For each batch, neighs.all should be a list containing 'cells_mapped' and 'distances' entries; nrow for each entry == n-cells in the batch. Consider recalculating using get_z_scaled_distances function.",
               fixed=TRUE
  )
  # neighs.all should be list w right colnames
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = as.character(c(2)), batch = "batch", n.neigh = 2,
                                                    neighs.all = list("0" = 1 , "all" = 1)),
               "Something is wrong with neighs.all argument. When batch is specified, it should be a list, named after batches. Consider recalculating using get_z_scaled_distances function.",
               fixed=TRUE
  )
  # neighs.all can be calculated w no errors
  neighs.all = get_z_scaled_distances(sce_correct_w_batch , batch = NULL)
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = as.character(c(2)), n.neigh = 2,
                                          neighs.all = neighs.all),
               NA
  )
  # neighs.all can be calculated w no errors
  neighs.all = get_z_scaled_distances(sce_correct_w_batch , batch = "batch")
  expect_error(evaluate_library(sce_correct_w_batch, genes.selection = as.character(c(2)), n.neigh = 2, batch = "batch",
                                          neighs.all = neighs.all),
               NA
  )
})

