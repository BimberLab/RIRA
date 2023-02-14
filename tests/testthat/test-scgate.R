library(Seurat)
library(SeuratData)
library(testthat)
library(dplyr)

# Note: the default of 2 workers is erroring on github actions
future::plan(future::sequential)

testthat::context("scGate")

test_that("scGates load", {
  gates <- GetAvailableScGates()
  expect_equal(length(gates), 23)

  gate <- GetScGateModel('demo_gate')
  expect_equal(length(gate$levels), 59)
})

test_that("scGate runs with custom models", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)

  seuratObj <- RunScGateForModels(seuratObj, modelNames = c('Bcell', 'Tcell', 'NK', 'Myeloid', 'Stromal', 'pDC', 'Erythrocyte', 'Epithelial', 'Platelet_MK'), labelRename = list(Tcell = 'T_NK', NK = 'T_NK'))
  print('scGate runs with custom models')
  print(table(seuratObj$scGateConsensus))
  dat <- table(seuratObj$scGateConsensus)
  expect_equal(unname(dat[['Myeloid']]), 682)
})

test_that("scGate Runs", {
  gate <- GetScGateModel('demo_gate')

  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  
  # Try without reductions present:
  seuratObj <- RunScGate(seuratObj, gate)
  expect_equal(sum(seuratObj$is.pure == 'Pure'), 1489, info = 'Before DimRedux')

  # Try with aliasing of models:
  seuratObj <- RunScGateForModels(seuratObj, modelNames = c('Bcell', 'Tcell', 'NK', 'Myeloid'), labelRename = list(Tcell = 'T_NK', NK = 'T_NK'))
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)
  expect_equal(unname(dat[['Bcell']]), 244, info = 'With aliasing', tolerance = 2)
  expect_equal(unname(dat[['T_NK']]), 1663, info = 'With aliasing')

  expect_false('Tcell' %in% names(dat), info = 'With aliasing')
  expect_false('NK' %in% names(dat), info = 'With aliasing')

  # Add reductions so plotting will work:
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- suppressWarnings(Seurat::NormalizeData(seuratObj))
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, nfeatures = 2000)
  seuratObj <- Seurat::ScaleData(seuratObj)
  seuratObj <- Seurat::RunPCA(seuratObj, features = Seurat::VariableFeatures(object = seuratObj))
  seuratObj <- Seurat::FindNeighbors(seuratObj, dims = 1:10)
  seuratObj <- Seurat::FindClusters(seuratObj, resolution = 0.5, random.seed = GetSeed())

  seuratObj <- RunScGate(seuratObj, gate)
  expect_equal(sum(seuratObj$is.pure == 'Pure'), 1504, info = 'After DimRedux')

  #At least execute this code once, so overt errors are caught
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = 1:10)
  PlotImmuneMarkers(seuratObj, reductions = 'umap')
})


test_that("scGate works with built-in gates", {
  # Use with built-in gate:
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- RunScGate(seuratObj, model = 'Bcell')
  expect_equal(sum(seuratObj$is.pure == 'Pure'), 293)

})


test_that("scGates runs on all", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)

  seuratObj <- RunScGateWithDefaultModels(seuratObj)
  expect_false('Bcell.is.pure.level4' %in% names(seuratObj@meta.data))
  dat <- table(seuratObj$scGateConsensus)
  print('RunScGateWithDefaultModels, using dropAmbiguousConsensusValues = FALSE')
  print(dat)

  expect_equal(unname(dat[['Bcell,PanBcell']]), 270)

  # Now with ambiguous cleanup:
  seuratObj <- RunScGateWithDefaultModels(seuratObj, dropAmbiguousConsensusValues = TRUE)
  dat <- table(seuratObj$scGateConsensus)
  print('RunScGateWithDefaultModels, using dropAmbiguousConsensusValues = TRUE')
  print(dat)
  expect_false('MoMacDC,Myeloid' %in% names(dat))
  expect_equal(unname(dat[['Myeloid']]), 212)
})

test_that("scGate Runs", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)

  # Try with aliasing of models:
  seuratObj <- RunScGateForModels(seuratObj, modelNames = c('Bcell.RM', 'Tcell.RM', 'NK.RM', 'Myeloid.RM', 'AvEp.RM', 'Epithelial.RM', 'Erythrocyte.RM', 'pDC.RM', 'Stromal.RM'), labelRename = list(Tcell.RM = 'T_NK', NK.RM = 'T_NK'))
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)

  expected <- c(
    Bcell.RM = 341,
    Myeloid.RM = 679,
    T_NK = 1649,
    'Bcell.RM,T_NK' = 9
  )

  for (pop in names(expected)) {
    expect_equal(unname(dat[[pop]]), expected[[pop]], info = paste0('RM models: ', pop))
  }

  # Now use wrapper
  print('Running RunScGateWithRhesusModels')
  seuratObj <- RunScGateWithRhesusModels(seuratObj)
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)

  expected <- c(
    Bcell = 341,
    Myeloid = 682,
    T_NK = 1649,
    'Bcell,T_NK' = 9,
    Platelet = 13
  )

  for (pop in names(expected)) {
    expect_equal(unname(dat[[pop]]), expected[[pop]], info = 'RM models using wrapper')
  }

  print(sort(table(seuratObj$scGateRaw)))
  dat <- table(seuratObj$scGateRaw)

  expected <- c(
    NK.RM = 66,
    Myeloid.RM = 679,
    Tcell.RM = 1311,
    'NK.RM,Tcell.RM' = 272
  )

  for (pop in names(expected)) {
    expect_equal(unname(dat[[pop]]), expected[[pop]], info = 'RM models, raw calls')
  }
})