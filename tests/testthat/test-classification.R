library(RIRA)
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
  set.seed(RIRA::GetSeed())
  toKeep <- sample(1:ncol(seuratObj), size = 2000)[1:1000]
  seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[toKeep])
  
  seuratObj <- Seurat::NormalizeData(seuratObj)
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, nfeatures = 2000)
  seuratObj <- Seurat::ScaleData(seuratObj)
  seuratObj <- Seurat::RunPCA(seuratObj, features = Seurat::VariableFeatures(object = seuratObj))
  seuratObj <- Seurat::FindNeighbors(seuratObj, dims = 1:10)
  seuratObj <- Seurat::FindClusters(seuratObj, resolution = 0.5)
  
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = 1:10)
  Seurat::DimPlot(seuratObj, reduction = "umap", label = T)
  Seurat::FeaturePlot(seuratObj, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"), label = T)
  
  seuratObj$CellType <- NA
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 0] <- 'TorNK' #Naive CD4+ T'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 1] <- 'Myeloid' #'CD14 Mono'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 2] <- 'TorNK' #Memory CD4+'
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
  set.seed(RIRA::GetSeed())
  toKeep <- sample(1:ncol(seuratObj), size = 2000)[1001:2000]
  seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[toKeep])
  
  return(seuratObj)
}


test_that("Cell type classification works", {
  fn <- 'seurat3k.rds'
  if (file.exists(fn)) {
    seuratObj <- readRDS(fn)
  } else {
    seuratObj <- prepareTrainingData()  
    saveRDS(seuratObj, file = fn)
  }
  
  RIRA::TrainAllModels(seuratObj = seuratObj, celltype_column = 'CellType', n_cores = 2)
  
  # Use new data:
  seuratObj <- prepareTestData()
  seuratObj <- PredictCellTypeProbability(seuratObj = seuratObj)
  seuratObj <- AssignCellType(seuratObj = seuratObj)
  
  table(seuratObj$Classifier_Consensus_Celltype)
  expected <- list(
    'B' = 743,
    'TorNK' = 61,
    'Unknown' = 196
  )
  
  for (cellType in names(expected)) {
    testthat::expect_equal(sum(seuratObj$Classifier_Consensus_Celltype == cellType), expected[[cellType]])
  }

  # Note: this is incredibly slow, so omit from the test:
  #RIRA::InterpretModels()
})

