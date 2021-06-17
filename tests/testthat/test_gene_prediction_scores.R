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

# load mouse ds
data("sce_mouseEmbryo", package = "geneBasisR")



test_that("Return is the correct class", {
  # right class
  out = get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2)
  expect_s3_class(out, "data.frame")

  # right colnames exist
  out = get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2)
  out = colnames(out)
  out_expect = c("gene" , "gene_score")
  out = sum(out_expect %in% out)
  expect_equal(out, 2)

  # you will get a message if corr-thresh is too high
  expect_message(get_gene_prediction_scores(sce_mouseEmbryo, genes.selection = rownames(sce_mouseEmbryo) , n.neigh = 2 , corr_all.thresh = 1),
                 "No genes are retained for this corr_all.thresh. Consider decreasing the threshold."
                 )
})



test_that("Return of the correct output", {
  out = get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, genes.predict = rownames(sce_correct))
  out = out[, c("gene" , "gene_score")]
  out$gene_score = round(out$gene_score)
  out_expect = data.frame(gene = as.character(c(1:3)) , gene_score = c(1,1,1))
  expect_identical(out, out_expect)
})

test_that("Return of the correct output, run 2", {
  out = get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , n.neigh = 2, genes.predict = as.character(c(1,2)))
  out = out[, c("gene" , "gene_score")]
  out$gene_score = round(out$gene_score)
  out_expect = data.frame(gene = as.character(c(1:2)) , gene_score = c(1,1))
  expect_identical(out, out_expect)
})




test_that("Wrong input, sce", {
  # should be unique rownames
  expect_error(get_gene_prediction_scores(sce_wrong_rownames, genes.selection = rownames(sce_correct)),
               "sce should have unique rownames.",
               fixed=TRUE
  )
  # sce should be sce
  expect_error(get_gene_prediction_scores(logcounts(sce_correct), genes.selection = rownames(sce_correct)),
               "sce should be a SingleCellExperiment object.",
               fixed=TRUE
  )
})


test_that("Wrong input, genes.selection", {
  # genes.selection should be character
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = c(1,2,3)),
               "Check genes.selection - should be character vector",
               fixed=TRUE
  )
  # genes.selection should be a subset of rownames in sce_correct
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,2,3,6))),
               "Some gene names are missing from SCE.",
               fixed=TRUE
  )
  # genes.selection is not NULL
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = NULL),
               "Check genes.selection - should be character vector",
               fixed=TRUE
  )
  # is ok if genes.selections are ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2),
               NA)
})


test_that("Wrong input, genes.all", {
  # genes.selection should be character
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1:5)), genes.all = c(1,4,5)),
               "Check genes.all - should be character vector",
               fixed=TRUE
  )
  # genes.selection should be a subset of rownames in sce_correct
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1:5)), genes.all = as.character(c(1,2,3,6))),
               "Some gene names are missing from SCE.",
               fixed=TRUE
  )
  # genes.selection is not NULL
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1:5)) , genes.all = NULL),
               "Check genes.all - should be character vector",
               fixed=TRUE
  )
  # is ok if genes.selections are ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1:5)), genes.all = as.character(c(1,4,5))),
               NA)

  # is ok if genes.selections are ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1:5)), corr_all.thresh = 1, genes.all = as.character(c(1,4,5))),
               NA)
  # single gene in predict
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1:5)), corr_all.thresh = 0,
                                          genes.predict ="1", genes.all = as.character(c(1:5))),
               NA)

})



test_that("Wrong input, batch", {
  # batch is NULL or character
  expect_error(get_gene_prediction_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = 1),
               "Check batch - should be NULL or string",
               fixed=TRUE
  )
  # batch should be the field in sce
  expect_error(get_gene_prediction_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = "sample"),
               "Batch should be one the colnames in colData(sce).",
               fixed=TRUE
  )
  # no errors if batch correct
  out = get_gene_prediction_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = "batch")
  expect_error(get_gene_prediction_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 2, batch = "batch"),
               NA
  )
})




test_that("Wrong input, n.neigh", {
  # n.neigh - positive scalar > 1
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = -1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh - positive scalar > 1
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 1),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )
  # n.neigh - positive scalar > 1
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = "random"),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )

  # internal check - n.neigh can be 'all' but not for users
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = "all"),
               "Check n.neigh - should be positive integer > 1",
               fixed=TRUE
  )

  # n.neigh should be < min(size(batch)) - 1
  expect_error(get_gene_prediction_scores(sce_correct_w_batch, genes.selection = rownames(sce_correct_w_batch), n.neigh = 3, batch = "batch"),
               "Each batch should contain at least > n.neigh cells. Check your dataset or decrease n.neigh.",
               fixed=TRUE
  )

})


test_that("Wrong input, nPC.all", {
  # nPC.all - NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 0),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = -10),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 5.5),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = "all"),
               "Check nPC.all - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC -  NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = NULL),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 3),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 30000),
               NA)
})


test_that("Wrong input, nPC.selection", {
  # nPC.all - NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = 0),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = -10),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = 5.5),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC.all - NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = "all"),
               "Check nPC.selection - should be NULL or positive integer",
               fixed=TRUE
  )
  # nPC -  NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = NULL),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.selection = 3),
               NA)

  # nPC -  NULL or positive scalar
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct), n.neigh = 2 , nPC.all = 3 , nPC.selection = 30000),
               NA)
})



test_that("Wrong input, genes.predict", {
  # genes should be character
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , genes.predict = c(1:5)),
               "Check genes.predict - should be character vector",
               fixed=TRUE
  )
  # genes should be character
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , genes.predict = NULL),
               "Check genes.predict - should be character vector",
               fixed=TRUE
  )
  # genes.predict should be a subset of rownames in sce_correct
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , genes.predict = as.character(c(1,2,3,6))),
               "Some gene names are missing from SCE.",
               fixed=TRUE
  )
  # genes.predict should be a subset of rownames in sce_correct
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , genes.all = as.character(c(1:5)) , genes.predict = as.character(c(4,5))),
               NA
  )
  # genes.predict should be a subset of rownames in sce_correct
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = rownames(sce_correct) , genes.all = as.character(c(1:5)) , genes.predict = as.character(c(4))),
               NA
  )
})



test_that("Wrong input, method", {
  # method is of the right entity
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,4,5)) , method = NULL),
               "Check method - should be either 'spearman', 'pearson' or 'kendall'",
               fixed=TRUE
  )
  # method is of the right entity
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,4,5)) , method = 0),
               "Check method - should be either 'spearman', 'pearson' or 'kendall'",
               fixed=TRUE
  )
  # method is of the right entity
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,4,5)) , method = "random"),
               "Check method - should be either 'spearman', 'pearson' or 'kendall'",
               fixed=TRUE
  )
  #
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,4,5)) , method = "pearson"),
               NA
  )
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,4,5)) , method = "spearman"),
               NA
  )
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,4,5)) , method = "kendall"),
               NA
  )
})



test_that("Wrong input, corr_all.thresh", {
  # corr_all.thresh should be numeric
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,4,5)) , corr_all.thresh = NULL),
               "Check corr_all.thresh - should be numeric",
               fixed=TRUE
  )
  # corr_all.thresh should be numeric
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = "random"),
               "Check corr_all.thresh - should be numeric",
               fixed=TRUE
  )

  # corr_all.thresh numeric is ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = 0),
               NA
  )
  # corr_all.thresh numeric is ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = 1),
               NA
  )
  # corr_all.thresh numeric is ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = 1000),
               NA
  )
  # corr_all.thresh numeric is ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = -1000),
               NA
  )

})


test_that("Wrong input, gene_stat_all is of correct format", {
  # corr_all.thresh should be numeric
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(1,4,5)) , corr_all.thresh = NULL),
               "Check corr_all.thresh - should be numeric",
               fixed=TRUE
  )
  # corr_all.thresh should be numeric
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = "random"),
               "Check corr_all.thresh - should be numeric",
               fixed=TRUE
  )

  # corr_all.thresh numeric is ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = 0),
               NA
  )
  # corr_all.thresh numeric is ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = 1),
               NA
  )
  # corr_all.thresh numeric is ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = 1000),
               NA
  )
  # corr_all.thresh numeric is ok
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2,3,5)) , corr_all.thresh = -1000),
               NA
  )

})


test_that("gene_stat_all is of correct format", {

  # gene_stat_all can be NULL
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2),
               NA
  )
  # gene_stat_all can be NULL
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, gene_stat_all = NULL),
               NA
  )
  # gene_stat_all should be data.frame with right columns
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, gene_stat_all = c(1,3)),
               "gene_stat_all is of the wrong format - should contain fields gene and corr_all.",
               fixed=TRUE
  )

  # gene_stat_all should be data.frame with right columns
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, gene_stat_all = data.frame(gene = 1 , corr = 2)),
               "gene_stat_all is of the wrong format - should contain fields gene and corr_all.",
               fixed=TRUE
  )

  # gene_stat_all should be data.frame with right columns
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, gene_stat_all = data.frame(gene = 1 , corr_all = 2)),
               NA
  )

  # gene_stat_all should be data.frame with right columns
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, gene_stat_all = data.frame(gene = "1" , corr_all = 2)),
               NA
  )

  # gene_stat_all should be data.frame with right columns
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, gene_stat_all = data.frame(gene = "1" , corr_all = 2, random = "ok")),
               NA
  )

  # gene_stat_all should be data.frame with right columns
  expect_error(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, gene_stat_all = data.frame(gene = "0" , corr_all = 2, random = "ok")),
               NA
  )

  # gonna give message if stat_all is sort of wrong
  expect_message(get_gene_prediction_scores(sce_correct, genes.selection = as.character(c(2)), n.neigh = 2, gene_stat_all = data.frame(gene = "0" , corr_all = 2, random = "ok")),
               "No genes are retained for this corr_all.thresh. Consider decreasing the threshold.",
               fixed = TRUE
  )
})



