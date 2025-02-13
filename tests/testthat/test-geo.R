library(Seurat)
library(SeuratData)
library(testthat)

context("GEO")

test_that("GEO Files Exist", {
  expected <- c(
    "GSE277821_feature_reference.xlsx",
    "GSE277821_RIRA.All.ADT.counts.rds",
    "GSE277821_RIRA.All.MergedClonotypes.txt.gz",
    "GSE277821_RIRA.All.Metadata.rds",
    "GSE277821_RIRA.All.RNA.counts.rds",
    "GSE277821_RIRA.All.seurat.rds",
    "GSE277821_RIRA.Bcell.seurat.rds",
    "GSE277821_RIRA.Myeloid.seurat.rds",
    "GSE277821_RIRA.T_NK.seurat.rds"
  )

  fn <- RCurl::getURL('ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277821/suppl/', .opts = RCurl::curlOptions(ftplistonly=TRUE))
  fn <- sort(gsub(unlist(strsplit(fn, split = '\n')), pattern = '\r', replacement = ''))
  expect_equal(fn, expected)
})