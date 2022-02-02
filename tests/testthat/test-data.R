context("RawData")

test_that("Atlas loading works", {
  mat <- matrix(sample(x = 0:10, size = 100, replace = T), ncol = 10)
  colnames(mat) <- LETTERS[1:10]
  rownames(mat) <- paste0('Gene', LETTERS[1:10])
  mat <- Seurat::as.sparse(mat)
  
  fn <- './AtlasTemp'
  dir.create(fn)
  fn <- normalizePath(fn)

  SetAtlasDir(fn)

  version <- '0.0.1'
  fn <- paste0(fn, '/', version)
  dir.create(fn)

  DropletUtils::write10xCounts(x = mat, path = paste0(fn, '/counts'), overwrite = T)

  meta <- data.frame(Col1 = 1:ncol(mat))
  rownames(meta) <- colnames(mat)
  utils::write.table(meta, paste0(fn, '/meta.csv'), sep = ',')

  seuratObj <- GetRiraCountMatrix(version)
  expect_equal(mat, Seurat::GetAssayData(seuratObj, slot = 'counts'))

  unlink(fn, recursive = T)
})