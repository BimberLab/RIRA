library(RIRA)
library(SeuratData)

context("Classification")

prepareExampleData <- function(){
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)
  
  set.seed(RIRA::GetSeed())
  toKeep <- sample(1:ncol(seuratObj), size = 1000)
  seuratObj <- subset(seuratObj, cells = colnames(seuratObj)[toKeep])
  
  seuratObj <- Seurat::NormalizeData(seuratObj)
  seuratObj <- Seurat::FindVariableFeatures(seuratObj, nfeatures = 2000)
  seuratObj <- Seurat::ScaleData(seuratObj)
  seuratObj <- Seurat::RunPCA(seuratObj, features = Seurat::VariableFeatures(object = seuratObj))
  seuratObj <- Seurat::FindNeighbors(seuratObj, dims = 1:10)
  seuratObj <- Seurat::FindClusters(seuratObj, resolution = 0.5)
  
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = 1:10)
  Seurat::DimPlot(seuratObj, reduction = "umap")
  Seurat::FeaturePlot(seuratObj, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"), label = T)
  
  seuratObj$CellType <- NA
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 0] <- 'Naive CD4+ T'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 1] <- 'CD14+ Mono'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 2] <- 'Memory CD4+'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 3] <- 'B'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 4] <- 'CD8+ T'
  seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 5] <- 'FCGR3A+ Mono'
  
  # These no longer form clusters after downsample:
  #seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 6] <- 'NK'
  #seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 7] <- 'DC'
  #seuratObj$CellType[seuratObj$RNA_snn_res.0.5 == 8] <- 'Platelet'
  
  return(seuratObj)
}

test_that("Cell type classification works", {
  fn <- 'seurat3k.rds'
  if (file.exists(fn)) {
    seuratObj <- readRDS(fn)
  } else {
    seuratObj <- prepareExampleData()  
    saveRDS(seuratObj, file = fn)
  }
  
  RIRA::TrainAllModels(seuratObj = seuratObj, celltype_column = 'CellType', n_cores = 2)
  seuratObj <- PredictCellTypeProbability(seuratObj = seuratObj)
  seuratObj <- AssignCellType(seuratObj = seuratObj)
  
  table(table(seuratObj$Classifier_Consensus_Celltype))
  
})

