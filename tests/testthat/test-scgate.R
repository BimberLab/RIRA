library(Seurat)
library(SeuratData)
library(testthat)
library(dplyr)

# Note: the default of 2 workers is erroring on github actions
future::plan(future::sequential)

testthat::context("scGate")

source('installSeuratData.R')

test_that("scGates load", {
  gates <- GetAvailableScGates()
  expect_equal(length(gates), 26)

  gate <- GetScGateModel('demo_gate')
  expect_equal(length(gate$levels), 59)
})

test_that("scGate runs with custom models", {
  seuratObj <- getBaseSeuratData()

  seuratObj <- RunScGateForModels(seuratObj, modelNames = c('Bcell', 'Tcell', 'NK', 'Myeloid', 'Stromal', 'pDC', 'Erythrocyte', 'Epithelial', 'Platelet_MK'), labelRename = list(Tcell = 'T_NK', NK = 'T_NK'))
  print('scGate runs with custom models')
  print(table(seuratObj$scGateConsensus))
  dat <- table(seuratObj$scGateConsensus)

  # Intended to capture 2.4.x:
  if (packageVersion('UCell') < '2.5.0') {
    print('UCell version is less than 2.5.0')
    expect_equal(unname(dat[['Myeloid']]), 680)
  } else {
    expect_equal(unname(dat[['Myeloid']]), 578)
  }

})

test_that("scGate Runs", {
  gate <- GetScGateModel('demo_gate')

  seuratObj <- getBaseSeuratData()
  
  # Try without reductions present:
  seuratObj <- RunScGate(seuratObj, gate)
  expect_equal(sum(seuratObj$is.pure == 'Pure'), 1486, info = 'Before DimRedux', tolerance = 1)

  # Try with aliasing of models:
  seuratObj <- RunScGateForModels(seuratObj, modelNames = c('Bcell', 'Tcell', 'NK', 'Myeloid'), labelRename = list(Tcell = 'T_NK', NK = 'T_NK'))
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)
  expect_equal(unname(dat[['Bcell']]), 244, info = 'With aliasing', tolerance = 2)
  expect_equal(unname(dat[['T_NK']]), 1657, info = 'With aliasing', tolerance = 1)

  expect_false('Tcell' %in% names(dat), info = 'With aliasing')
  expect_false('NK' %in% names(dat), info = 'With aliasing')

  # Add reductions so plotting will work:
  seuratObj <- getBaseSeuratData()
  seuratObj <- suppressWarnings(Seurat::NormalizeData(seuratObj))
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, nfeatures = 2000)
  seuratObj <- Seurat::ScaleData(seuratObj)
  seuratObj <- Seurat::RunPCA(seuratObj, features = Seurat::VariableFeatures(object = seuratObj))
  seuratObj <- Seurat::FindNeighbors(seuratObj, dims = 1:10)
  seuratObj <- Seurat::FindClusters(seuratObj, resolution = 0.5, random.seed = GetSeed())

  seuratObj <- RunScGate(seuratObj, gate)
  if (packageVersion('UCell') < '2.5.0') {
    expect_equal(sum(seuratObj$is.pure == 'Pure'), 1505, info = 'After DimRedux')
  } else {
    expect_equal(sum(seuratObj$is.pure == 'Pure'), 1493, info = 'After DimRedux', tolerance = 1)
  }

  #At least execute this code once, so overt errors are caught
  seuratObj <- suppressWarnings(Seurat::RunUMAP(seuratObj, dims = 1:10))
  PlotImmuneMarkers(seuratObj, reductions = 'umap')
})


test_that("scGate works with built-in gates", {
  # Use with built-in gate:
  seuratObj <- getBaseSeuratData()
  seuratObj <- RunScGate(seuratObj, model = 'Bcell')
  expect_equal(sum(seuratObj$is.pure == 'Pure'), 340, tolerance = 1)

})


test_that("scGates runs on all", {
  seuratObj <- getBaseSeuratData()

  seuratObj <- RunScGateWithDefaultModels(seuratObj)
  expect_false('Bcell.is.pure.level4' %in% names(seuratObj@meta.data))
  dat <- table(seuratObj$scGateConsensus)
  print('RunScGateWithDefaultModels, using dropAmbiguousConsensusValues = FALSE')
  print(dat)

  expect_equal(unname(dat[['Bcell,Bcell.NonGerminalCenter,Immune,PanBcell']]), 331, tolerance = 1)

  # Now with ambiguous cleanup:
  seuratObj <- RunScGateWithDefaultModels(seuratObj, dropAmbiguousConsensusValues = TRUE)
  dat <- table(seuratObj$scGateConsensus)
  print('RunScGateWithDefaultModels, using dropAmbiguousConsensusValues = TRUE')
  print(dat)
  expect_false('MoMacDC,Myeloid' %in% names(dat))
  if (packageVersion('UCell') < '2.5.0') {
    expect_equal(unname(dat[['Immune']]), 7, tolerance = 1)
  } else {
    expect_equal(unname(dat[['Immune']]), 25, tolerance = 1)
  }
})

test_that("scGate Runs", {
  seuratObj <- getBaseSeuratData()

  # Try with aliasing of models:
  seuratObj <- RunScGateForModels(seuratObj, modelNames = c('Bcell.RM', 'Tcell.RM', 'NK.RM', 'Myeloid.RM', 'AvEp.RM', 'Epithelial.RM', 'Erythrocyte.RM', 'pDC.RM', 'Stromal.RM'), labelRename = list(Tcell.RM = 'T_NK', NK.RM = 'T_NK'))
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)

  if (packageVersion('UCell') < '2.5.0') {
    expected <- c(
      Bcell.RM = 337,
      Myeloid.RM = 678,
      T_NK = 1648,
      'Bcell.RM,T_NK' = 15
    )
  } else {
    expected <- c(
      Bcell.RM = 337,
      Myeloid.RM = 676,
      T_NK = 1647,
      'Bcell.RM,T_NK' = 15
    )
  }

  for (pop in names(expected)) {
    expect_equal(unname(dat[[pop]]), expected[[pop]], info = paste0('RM models: ', pop), tolerance = 3)
  }

  # Now use wrapper
  print('Running RunScGateWithRhesusModels')
  seuratObj <- RunScGateWithRhesusModels(seuratObj)
  print(sort(table(seuratObj$scGateConsensus)))
  dat <- table(seuratObj$scGateConsensus)

  expected <- c(
    Bcell = 337,
    Myeloid = 680,
    T_NK = 1647,
    'Bcell,T_NK' = 15,
    Platelet = 14
  )

  for (pop in names(expected)) {
    expect_equal(unname(dat[[pop]]), expected[[pop]], info = 'RM models using wrapper', tolerance = 4)
  }

  print(sort(table(seuratObj$scGateRaw)))
  dat <- table(seuratObj$scGateRaw)

  expected <- c(
    NK.RM = 72,
    Myeloid.RM = 679,
    Tcell.RM = 1301,
    'NK.RM,Tcell.RM' = 274
  )

  for (pop in names(expected)) {
    expect_equal(unname(dat[[pop]]), expected[[pop]], info = 'RM models, raw calls', tolerance = 4)
  }
})