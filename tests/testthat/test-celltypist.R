library(SeuratData)
library(testthat)

testthat::context("celltypist")

test_that("celltypist runs", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  
  seuratObj <- RIRA::RunCellTypist(seuratObj)
  
  expect_equal(9, length(unique(seuratObj$majority_voting)), info = 'using default model')
  expect_equal(82, length(unique(seuratObj$predicted_labels)))
  expect_equal(254, unname(table(seuratObj$predicted_labels)['B cells']))

  # NOTE: this is very slow, so skip in automated testing for now
  modelFile <- 'myModel.pkl'
  seuratObjForTraining <- subset(seuratObj, cells = sort(rownames(seuratObj@meta.data))[1:500])
  TrainCellTypist(seuratObjForTraining, 'majority_voting', paste0(getwd(), '/', modelFile))

  # A bit circular, but this is just a test case:
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  seuratObj <- RIRA::RunCellTypist(seuratObj, modelName = modelFile)

  print(table(seuratObj$majority_voting))
  print(table(seuratObj$predicted_labels))
  expect_equal(8, length(unique(seuratObj$majority_voting)), info = 'using custom model')
  expect_equal(54, length(unique(seuratObj$predicted_labels)), tolerance = 3)
  expect_equal(356, unname(table(seuratObj$predicted_labels)['B cells']), tolerance = 2)
})