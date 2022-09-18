library(Seurat)
library(SeuratData)

testthat::context("Classification")

getBaseSeuratData <- function(){
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)

  return(seuratObj)
}

prepareTrainingData <- function(){
  seuratObj <- getBaseSeuratData()
  set.seed(GetSeed())
  toKeep <- sample(1:ncol(seuratObj), size = 2000)[1:1000]
  seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[toKeep])
  
  seuratObj <- Seurat::NormalizeData(seuratObj)
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, nfeatures = 2000)
  seuratObj <- Seurat::ScaleData(seuratObj)
  seuratObj <- Seurat::RunPCA(seuratObj, features = Seurat::VariableFeatures(object = seuratObj))
  seuratObj <- Seurat::FindNeighbors(seuratObj, dims = 1:10)
  seuratObj <- Seurat::FindClusters(seuratObj, resolution = 0.5)
  
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = 1:10)
  print(Seurat::DimPlot(seuratObj, reduction = "umap", label = T))
  print(Seurat::FeaturePlot(seuratObj, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"), label = T))
  
  seuratObj$CellType <- NA
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 0] <- 'TorNK' #Naive CD4+ T'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 1] <- 'TorNK' #Memory CD4+'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 2] <- 'Myeloid' #'CD14 Mono'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 3] <- 'B'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 4] <- 'TorNK' #CD8+ T'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 5] <- 'Myeloid' #'FCGR3A Mono'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 6] <- 'TorNK'
  
  # These no longer form clusters after downsample:
  #seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 6] <- 'NK'
  #seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 7] <- 'DC'
  #seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 8] <- 'Platelet'
  
  return(seuratObj)
}

prepareTestData <- function(){
  seuratObj <- getBaseSeuratData()
  set.seed(GetSeed())
  toKeep <- sample(1:ncol(seuratObj), size = 2000)[1001:2000]
  seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[toKeep])
  
  return(seuratObj)
}

test_that("Cell type classification works", {
  fn <- 'seurat3k.rds'
  if (file.exists(fn)) {
    seuratObjTrain <- readRDS(fn)
  } else {
    seuratObjTrain <- prepareTrainingData()  
    saveRDS(seuratObjTrain, file = fn)
  }

  TrainModelsFromSeurat(seuratObj = seuratObjTrain, celltype_column = 'CellType', n_cores = 2, output_dir = "./testClassifiers")
  modelFiles <- list.files("./testClassifiers/models/")
  expect_equal(length(modelFiles), 3)

  names(modelFiles) <- sapply(modelFiles, function(x){
    return(tools::file_path_sans_ext(basename(x)))
  })

  print(modelFiles)

  # Use new data:
  seuratObj <- prepareTestData()
  seuratObj <- PredictCellTypeProbability(seuratObj = seuratObj, models = modelFiles)
  
  print(sort(table(seuratObj$RIRA_Consensus)))
  expected <- list(
    'B' = 104,
    'TorNK' = 609,
    'Myeloid' = 231,
    'Unknown' = 56
  )
  
  for (cellType in names(expected)) {
    testthat::expect_equal(sum(seuratObj$RIRA_Consensus == cellType), expected[[cellType]])
  }

  # Test batchSize:
  seuratObj2 <- prepareTestData()
  seuratObj2 <- PredictCellTypeProbability(seuratObj = seuratObj2, models = modelFiles, batchSize = 38)
  testthat::expect_equal(0, sum(seuratObj$RIRA_Consensus != seuratObj2$RIRA_Consensus))
  
  feats <- c("IFI30", "CD7", "CD3E",   "MS4A1", "CD79A",   "VCAN", "MNDA",   "C1QB", "C1QA")
  TrainModelsFromSeurat(seuratObj = seuratObjTrain, celltype_column = 'CellType', n_cores = 2, output_dir = './classifiers2', gene_list = feats)
  
  # Note: this is incredibly slow, so use the feature-limited version:
  InterpretModels(output_dir = './classifiers2')
})

test_that("Built-in models work", {
  seuratObj <- prepareTestData()
  seuratObj <- PredictCellTypeProbability(seuratObj = seuratObj, models = list(
    'T_LR' = 'RIRA_CD4vCD8_LR'
  ), fieldName = 'RIRA_TCell')

  print(sort(table(seuratObj$RIRA_TCell)))
})