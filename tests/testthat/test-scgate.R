library(Seurat)
library(SeuratData)
library(testthat)
library(dplyr)

testthat::context("scGate")

test_that("scGates load", {
  gates <- GetAvailableScGates()
  expect_equal(length(gates), 17)

  gate <- GetScGateModel('demo_gate')
  expect_equal(length(gate$levels), 59)
})

test_that("scGate runs with custom models", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)

  # Try with aliasing of models:
  seuratObj <- RIRA::RunScGateForModels(seuratObj, modelNames = c('Bcell', 'Tcell', 'NK', 'Myeloid', 'Stromal', 'pDC', 'Erythrocyte', 'Epithelial', 'Platelet_MK'), labelRename = list(Tcell = 'T_NK', NK = 'T_NK'))
  print(table(seuratObj$scGateConsensus))
  dat <- table(seuratObj$scGateConsensus)
  expect_equal(unname(dat[['pDC']]), 1)
})

test_that("scGate Runs", {
  gate <- GetScGateModel('demo_gate')

  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  
  # Try without reductions present:
  seuratObj <- RunScGate(seuratObj, gate)
  expect_equal(sum(seuratObj$is.pure == 'Pure'), 1478, info = 'Before DimRedux')

  # Try with aliasing of models:
  seuratObj <- RIRA::RunScGateForModels(seuratObj, modelNames = c('Bcell', 'Tcell', 'NK', 'Myeloid'), labelRename = list(Tcell = 'T_NK', NK = 'T_NK'))
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)
  expect_equal(unname(dat[['Bcell']]), 240, info = 'With aliasing')
  expect_equal(unname(dat[['T_NK']]), 1660, info = 'With aliasing')

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
  expect_equal(sum(seuratObj$is.pure == 'Pure'), 1478, info = 'After DimRedux')

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
  expect_equal(sum(seuratObj$is.pure == 'Pure'), 242)
  
})


test_that("scGates runs on all", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  
  seuratObj <- RunScGateWithDefaultModels(seuratObj)
  expect_false('Bcell.is.pure.level4' %in% names(seuratObj@meta.data))
  dat <- table(seuratObj$scGateConsensus)
  dat
  
  expect_equal(unname(dat[['Bcell,PanBcell']]), 239)
  expect_equal(unname(dat[['NK']]), 69)
  
  # Now with ambiguous cleanup:
  seuratObj <- RunScGateWithDefaultModels(seuratObj, dropAmbiguousConsensusValues = TRUE)
  dat <- table(seuratObj$scGateConsensus)
  dat
  expect_equal(unname(dat[['NK']]), 67)
  expect_equal(unname(dat[['Myeloid']]), 18)
})

test_that("scGate Runs", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)

  # Try with aliasing of models:
  seuratObj <- RIRA::RunScGateForModels(seuratObj, modelNames = c('Bcell.RM', 'Tcell.RM', 'NK.RM', 'Myeloid.RM', 'AvEp.RM', 'Epithelial.RM', 'Erythrocyte.RM', 'pDC.RM', 'Stromal.RM'), labelRename = list(Tcell.RM = 'T_NK', NK.RM = 'T_NK'))
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)

  expected <- c(
    Bcell.RM = 339,
    Myeloid.RM = 668,
    T_NK = 1647,
    'Bcell.RM,T_NK' = 14
  )

  for (pop in names(expected)) {
    expect_equal(unname(dat[[pop]]), expected[[pop]], info = 'RM models')
  }

  # Now use wrapper
  seuratObj <- RIRA::RunScGateWithRhesusModels(seuratObj)
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)

  expected <- c(
    Bcell = 338,
    MoMacDC = 653,
    T_NK = 1647,
    'Bcell,T_NK' = 14
  )

  for (pop in names(expected)) {
    expect_equal(unname(dat[[pop]]), expected[[pop]], info = 'RM models using wrapper')
  }
})