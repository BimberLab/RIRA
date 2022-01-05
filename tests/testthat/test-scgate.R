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

  seuratObj <- RunScGate(seuratObj, gate)
  expect_equal(sum(seuratObj$is.pure.level1 == 'Pure'), 2343)
})