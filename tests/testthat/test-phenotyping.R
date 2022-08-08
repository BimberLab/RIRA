library(Seurat)
library(SeuratData)
library(testthat)

context("Phenotyping")

test_that("Gene sets work", {
  expect_equal(length(GetGeneSet('MMul10TcrGenes')), 107)
  expect_equal(length(GetGeneSet('TandNK_Activation.1')), 8)
})

test_that("ExpandGeneList works", {
  expect_equal(length(ExpandGeneList(c('MMul10TcrGenes', 'TBX21', 'CD8A'))), 109)
  expect_equal(length(ExpandGeneList(c('MMul10TcrGenes', 'TBX21', 'CD8A', 'MMul10_MHC'))), 144)
})

test_that("ScoreUsingSavedComponent works", {
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)

  df <- ExtractGeneWeights(seuratObj, componentNum = 1)

  seuratObj <- ScoreUsingSavedComponent(seuratObj, componentOrName = df, fieldName = 'Test1')
  seuratObj <- ScoreUsingSavedComponent(seuratObj, componentOrName = 'Tcell_NaiveToEffector', fieldName = 'Test2')
})