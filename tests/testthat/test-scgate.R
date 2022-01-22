library(SeuratData)
library(testthat)

testthat::context("scGate")

test_that("scGates load", {
  gates <- GetAvailableScGates()
  expect_equal(length(gates), 1)
  
  gate <- GetScGateModel('demo_gate')
  expect_equal(length(gate$levels), 61)
})


test_that("scGate Runs", {
  gate <- GetScGateModel('demo_gate')

  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  
  # Try without reductions present:
  seuratObj <- RunScGate(seuratObj, gate)
  expect_equal(sum(seuratObj$is.pure.level1 == 'Pure'), 2336, info = 'Before DimRedux')
  
  # Add reductions so plotting will work:
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- Seurat::NormalizeData(seuratObj)
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, nfeatures = 2000, seed.use = GetSeed())
  seuratObj <- Seurat::ScaleData(seuratObj)
  seuratObj <- Seurat::RunPCA(seuratObj, features = Seurat::VariableFeatures(object = seuratObj))
  seuratObj <- Seurat::FindNeighbors(seuratObj, dims = 1:10)
  seuratObj <- Seurat::FindClusters(seuratObj, resolution = 0.5, random.seed = GetSeed())
  
  seuratObj <- RunScGate(seuratObj, gate)
  expect_equal(sum(seuratObj$is.pure.level1 == 'Pure'), 2343, info = 'After DimRedux')
})


test_that("scGate works with build-in gates", {
  # Use with built-in gate:
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  seuratObj <- RunScGate(seuratObj, model = 'Bcell')
  expect_equal(sum(seuratObj$is.pure.level4 == 'Pure'), 245)
  
})


test_that("scGates runs on all", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  
  seuratObj <- RunScGateWithDefaultModels(seuratObj)
  expect_false('Bcell.is.pure.level4' %in% names(seuratObj@meta.data))
  dat <- table(seuratObj$scGateConsensus)
  dat
  
  expect_equal(unname(dat[['Bcell,PanBcell']]), 244)
  expect_equal(unname(dat[['NK']]), 161)
})

