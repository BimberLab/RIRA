library(Seurat)
library(SeuratData)
library(testthat)

testthat::context("celltypist")

test_that("celltypist is installed", {
    # this provides more useful information if there is an error loading celltypist
    reticulate::import("celltypist")
})

test_that("celltypist runs", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  print(seuratObj)

  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  seuratObj <- RunCellTypist(seuratObj)

  # This appears to differ in bioconductor 3.15 vs devel
  # This should be identical to the test below
  expect_equal(11, length(unique(seuratObj$majority_voting)), info = 'using default model', tolerance = 1)
  expect_equal(110, length(unique(seuratObj$predicted_labels)))
  expect_equal(289, unname(table(seuratObj$predicted_labels)['B cells']))

  # ensure the RIRA model works:
  seuratObj <- RunCellTypist(seuratObj, modelName = 'RIRA_Immune_v1', columnPrefix = 'RIRA.')
  print(sort(table(seuratObj$RIRA.majority_voting)))
  expect_equal(4, length(unique(seuratObj$RIRA.majority_voting)), info = 'using RIRA model', tolerance = 1)
  expect_equal(346, unname(table(seuratObj$RIRA.majority_voting)['Bcell']), tolerance = 1)
  expect_equal(1653, unname(table(seuratObj$RIRA.majority_voting)['T_NK']), tolerance = 1)
  expect_equal(686, unname(table(seuratObj$RIRA.majority_voting)['MoMacDC']))

  # NOTE: this is very slow, so skip in automated testing for now
  modelFile <- 'myModel.pkl'
  cells <- rownames(seuratObj@meta.data)[!is.na(seuratObj$majority_voting)]
  seuratObjForTraining <- subset(seuratObj, cells = cells)
  seuratObjForTraining <- subset(seuratObjForTraining, cells = sort(rownames(seuratObjForTraining@meta.data))[1:500])
  print(paste0('cells for training: ', ncol(seuratObjForTraining)))
  TrainCellTypist(seuratObjForTraining, 'majority_voting', paste0(getwd(), '/', modelFile))

  # A bit circular, but this is just a test case:
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  seuratObj <- RunCellTypist(seuratObj, modelName = modelFile)

  print(table(seuratObj$majority_voting))
  print(table(seuratObj$predicted_labels))
  expect_equal(10, length(unique(seuratObj$majority_voting)), info = 'using custom model', tolerance = 1)
  expect_equal(54, length(unique(seuratObj$predicted_labels)), tolerance = 3)
  expect_equal(356, unname(table(seuratObj$predicted_labels)['B cells']), tolerance = 2)
})

test_that("celltypist runs with batchSize", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  print(seuratObj)

  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  seuratObj <- RunCellTypist(seuratObj, maxBatchSize = 250)

  print('cell class:')
  print(table(seuratObj$cellclass))

  print('majority voting:')
  print(table(seuratObj$majority_voting))

  # This should be identical to the test above
  expect_equal(12, length(unique(seuratObj$cellclass)), info = 'using default model', tolerance = 0)
  expect_equal(28, length(unique(seuratObj$majority_voting)), info = 'using default model', tolerance = 1)  # NOTE: getting different outcomes on devel vs. 3.16, perhaps due to some package difference? the difference is ambugious calls
  expect_equal(110, length(unique(seuratObj$predicted_labels)))
  expect_equal(289, unname(table(seuratObj$predicted_labels)['B cells']))
})

test_that("celltypist runs for RIRA models", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  seuratObj <- Classify_TNK(seuratObj, retainProbabilityMatrix = TRUE)

  print(table(seuratObj$RIRA_TNK_v2.predicted_labels))
  
  expect_equal(4, length(unique(seuratObj$RIRA_TNK_v2.predicted_labels)), info = 'using RIRA T_NK', tolerance = 1)
  expect_equal(215, unname(table(seuratObj$RIRA_TNK_v2.predicted_labels)['CD4+ T Cells']), tolerance = 1)
  expect_equal(710, unname(table(seuratObj$RIRA_TNK_v2.predicted_labels)['CD8+ T Cells']), tolerance = 1)
  expect_equal(118, unname(table(seuratObj$RIRA_TNK_v2.predicted_labels)['NK Cells']), tolerance = 1)
  expect_equal(1657, unname(table(seuratObj$RIRA_TNK_v2.predicted_labels)['Other']), tolerance = 1)

  expect_equal(6.64e-08, min(seuratObj$RIRA_TNK_v2.prob.NK.Cells), tolerance = 0.00001)
})
