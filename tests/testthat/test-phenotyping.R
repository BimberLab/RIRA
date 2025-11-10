library(Seurat)
library(SeuratData)
library(testthat)

context("Phenotyping")

source('installSeuratData.R')

test_that("Gene sets work", {
  expect_equal(length(GetGeneSet('MMul10TcrGenes')), 102)
  expect_equal(length(GetGeneSet('TandNK_Activation.1')), 8)
})

test_that("ExpandGeneList works", {
  expect_equal(length(ExpandGeneList(c('MMul10TcrGenes', 'TBX21', 'CD8A'))), 104)
  expect_equal(length(ExpandGeneList(c('MMul10TcrGenes', 'TBX21', 'CD8A', 'MMul10_MHC'))), 139)
})

test_that("ScoreUsingSavedComponent works", {
  set.seed(GetSeed())
  seuratObj <- getBaseSeuratData()
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, nfeatures = 2000)
  seuratObj <- Seurat::ScaleData(seuratObj)
  seuratObj <- Seurat::RunPCA(seuratObj, features = Seurat::VariableFeatures(object = seuratObj))

  df <- ExtractGeneWeights(seuratObj, componentNum = 1)
  expect_equal(nrow(df), 200)

  seuratObj <- ScoreUsingSavedComponent(seuratObj, componentOrName = df, fieldName = 'Test1')
  # NOTE: the sign of this component seems to flip-flop
  expect_equal(round(min(abs(seuratObj$Test1)), 2), 2.52)
  expect_equal(round(max(abs(seuratObj$Test1)), 2), 199.21)
  
  seuratObj <- ScoreUsingSavedComponent(seuratObj, componentOrName = 'Tcell_NaiveToEffector', fieldName = 'Test2')
  expect_equal(round(max(seuratObj$Test2), 2), 41.61)
  expect_equal(round(min(seuratObj$Test2), 2), -6.69)

  seuratObj <- CalculateUCellScores(seuratObj)
  expect_true('TandNK_Activation_UCell' %in% names(seuratObj@meta.data))
  expect_true('MAIT_Markers_UCell' %in% names(seuratObj@meta.data))
})