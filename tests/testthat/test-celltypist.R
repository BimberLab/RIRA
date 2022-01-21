library(SeuratData)
library(testthat)

testthat::context("celltypist")

test_that("celltypist runs", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- Seurat::NormalizeData(seuratObj)
  
  seuratObj <- RIRA::RunCellTypist(seuratObj)
  
  expect_equal(9, length(unique(seuratObj$majority_voting)))
  expect_equal(30, length(unique(seuratObj$predicted_labels)))
  expect_equal(323, unname(table(seuratObj$predicted_labels)['B cells']))
})