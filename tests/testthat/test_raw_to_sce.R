# context("Testing raw_to_sce")
library(geneBasisR)

### set up inputs

## libs
library(SingleCellExperiment)
set.seed(32)

counts_dir = system.file("extdata", "raw_spleen.txt", package = "geneBasisR")
meta_dir = system.file("extdata", "raw_spleen_meta.txt", package = "geneBasisR")

counts_dir_wrong = "/Users/alsu/Develop/geneBasisR/inst/extdata/raw_spleen_wrong.txt"
meta_dir_wrong = "/Users/alsu/Develop/geneBasisR/inst/extdata/raw_spleen_meta_wrong.txt"


test_that("Return is the correct class", {
  out = raw_to_sce(counts_dir , counts_type = "logcounts" )
  expect_s4_class(out, "SingleCellExperiment")

  out = raw_to_sce(counts_dir , counts_type = "counts" )
  expect_s4_class(out, "SingleCellExperiment")

  out = raw_to_sce(counts_dir , counts_type = "counts" , transform_counts_to_logcounts = F)
  expect_s4_class(out, "SingleCellExperiment")

  out = raw_to_sce(counts_dir , counts_type = "counts" , transform_counts_to_logcounts = T)
  expect_s4_class(out, "SingleCellExperiment")

  out = raw_to_sce(counts_dir , counts_type = "logcounts" , meta_dir = meta_dir)
  expect_s4_class(out, "SingleCellExperiment")

  out = raw_to_sce(counts_dir , counts_type = "counts" , transform_counts_to_logcounts = T, meta_dir = meta_dir)
  expect_s4_class(out, "SingleCellExperiment")

  out = raw_to_sce(counts_dir , counts_type = "counts" , transform_counts_to_logcounts = F, meta_dir = meta_dir)
  expect_s4_class(out, "SingleCellExperiment")

  out = raw_to_sce(counts_dir , counts_type = "logcounts" , meta_dir = meta_dir, batch = "celltype")
  expect_s4_class(out, "SingleCellExperiment")

  out = raw_to_sce(counts_dir , counts_type = "counts" , meta_dir = meta_dir, batch = "celltype", d = 10)
  expect_s4_class(out, "SingleCellExperiment")
})


test_that("if counts_dir does not exist, the correct error is returned", {

})



