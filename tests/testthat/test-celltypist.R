library(Seurat)
library(SeuratData)
library(testthat)

testthat::context("celltypist")

test_that("celltypist is installed", {
  # this provides more useful information if there is an error loading celltypist
  testthat::expect_no_condition(print(reticulate::py_exe()))
  testthat::expect_no_condition(print(reticulate::py_version()))
  testthat::expect_no_condition(print(reticulate::import('celltypist')))
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
  seuratObj <- RunCellTypist(seuratObj, modelName = 'RIRA_Immune_v2', columnPrefix = 'RIRA.')
  print(sort(table(seuratObj$RIRA.majority_voting)))
  expect_equal(4, length(unique(seuratObj$RIRA.majority_voting)), info = 'using RIRA model', tolerance = 1)
  expect_equal(346, unname(table(seuratObj$RIRA.majority_voting)['Bcell']), tolerance = 1)
  expect_equal(1653, unname(table(seuratObj$RIRA.majority_voting)['T_NK']), tolerance = 1)
  expect_equal(570, unname(table(seuratObj$RIRA.majority_voting)['Myeloid']))

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
  expect_equal(13, length(unique(seuratObj$cellclass)), info = 'using default model', tolerance = 0)
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

test_that("FilterDisallowedClasses works as expected", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)

  seuratObj <- RunScGateWithRhesusModels(seuratObj)
  seuratObj <- Classify_ImmuneCells(seuratObj)

  print('RIRA_Immune_v2.cellclass:')
  print(table(seuratObj$RIRA_Immune_v2.cellclass))

  expect_equal(255, sum(seuratObj$RIRA_Immune_v2.cellclass == 'Bcell', na.rm = T))
  expect_equal(570, sum(seuratObj$RIRA_Immune_v2.cellclass == 'Myeloid', na.rm = T))
  expect_equal(1289, sum(seuratObj$RIRA_Immune_v2.cellclass == 'T_NK', na.rm = T))

  print('DisallowedUCellCombinations:')
  print(table(seuratObj$DisallowedUCellCombinations))

  # NOTE: these are producing different results on 3.16 vs devel. This is possibly scGate versions?
  expect_equal(347, sum(seuratObj$DisallowedUCellCombinations == 'NeutrophilLineage.RM_UCell', na.rm = T), tolerance = 3)
  expect_equal(21, sum(seuratObj$DisallowedUCellCombinations == 'Erythrocyte.RM_UCell', na.rm = T))
  expect_equal(57, sum(seuratObj$DisallowedUCellCombinations == 'NK.RM_UCell', na.rm = T))
  expect_equal(55, sum(seuratObj$DisallowedUCellCombinations == 'Platelet.RM_UCell', na.rm = T), tolerance = 1)
})
