library(Seurat)
library(SeuratData)

testthat::context("Classification")

source('installSeuratData.R')

SimulateSeuratData <- function() {
  set.seed(GetSeed())

  #union features from all activation component files
  features <- c()
  
  #grab CD8 component files
  for (i in 1:6) {
    savedFile <- system.file(paste0("components/CD8_Activation_sPLSDA_component", i, ".tsv"), package = "RIRA")
    component <- utils::read.table(savedFile, header = T, sep = '\t')
    features <- c(features, component$feature)
  }
  
  #grab CD4 component files
  for (i in 1:6) {
    savedFile <- system.file(paste0("components/CD4_Activation_sPLSDA_component", i, ".tsv"), package = "RIRA")
    component <- utils::read.table(savedFile, header = T, sep = '\t')
    features <- c(features, component$feature)
  }
  
  #grab General component files (9 components)
  for (i in 1:9) {
    savedFile <- system.file(paste0("components/GeneralizedTCR_Activation_sPLSDA_component", i, ".tsv"), package = "RIRA")
    component <- utils::read.table(savedFile, header = T, sep = '\t')
    features <- c(features, component$feature)
  }
  
  features <- unique(features)

  #simulate a seurat object with all features
  mtx <- matrix(
    rpois(200 * length(features), lambda = 2),
    nrow = 200,
    ncol = length(features)
  )

  colnames(mtx) <- features
  rownames(mtx) <- paste0("Cell_", seq_len(nrow(mtx)))
  
  #suppress warnings about dgCMatrix coercion
  seuratObj <- suppressWarnings(CreateSeuratObject(
    counts = t(mtx),
    min.cells = 3,
    min.features = 200
  ))

  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  seuratObj <- Seurat::ScaleData(seuratObj, features = features, verbose = FALSE)
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

  seuratObj <- suppressWarnings(Seurat::RunUMAP(seuratObj, dims = 1:10))
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
  seuratObjTrain <- prepareTrainingData()

  TrainModelsFromSeurat(seuratObj = seuratObjTrain, celltype_column = 'CellType', n_cores = 1, output_dir = "./testClassifiers")
  modelFiles <- paste0("./testClassifiers/models/", list.files("./testClassifiers/models/"))
  expect_equal(length(modelFiles), 3)

  names(modelFiles) <- sapply(modelFiles, function(x){
    return(tools::file_path_sans_ext(basename(x)))
  })

  # Use new data:
  seuratObj <- prepareTestData()
  seuratObj <- PredictCellTypeProbability(seuratObj = seuratObj, models = modelFiles)

  print(sort(table(seuratObj$RIRA_Consensus)))
  expected <- list(
    'B' = 109,
    'TorNK' = 607,
    'Myeloid' = 230,
    'Unknown' = 54
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


test_that("Predict TCell Activation works ", {

  SimulateSeuratObj <- SimulateSeuratData()
  SimulateSeuratObj <- PredictTcellActivation(SimulateSeuratObj)
  SimulateSeuratObj <- PredictTcellActivation(SimulateSeuratObj, modelName = 'CD4')
  SimulateSeuratObj <- PredictTcellActivation(SimulateSeuratObj, modelName = 'CD8')

  #exhaustive check on component scores, since all are necessary for the default parameterization of the function
  for (i in 1:6) {
    colName <- paste0("CD8_sPLSDA_Score_", i)
    testthat::expect_true(colName %in% names(SimulateSeuratObj@meta.data))
  }
  
  for (i in 1:6) {
    colName <- paste0("CD4_sPLSDA_Score_", i)
    testthat::expect_true(colName %in% names(SimulateSeuratObj@meta.data), label = paste0('Missing: colName'))
  }
  
  for (i in 1:9) {
    colName <- paste0("GeneralizedTCR_sPLSDA_Score_", i)
    testthat::expect_true(colName %in% names(SimulateSeuratObj@meta.data))
  }
  
  #check that the scores are numeric
  testthat::expect_type(SimulateSeuratObj@meta.data[5, 'CD8_sPLSDA_Score_2'], "double")
  testthat::expect_type(SimulateSeuratObj@meta.data[5, 'CD4_sPLSDA_Score_3'], "double")
  testthat::expect_type(SimulateSeuratObj@meta.data[5, 'GeneralizedTCR_sPLSDA_Score_7'], "double")

  #check that prediction portion of the function worked

  #probabilities
  testthat::expect_equal(4, sum(grepl("^CD8_sPLSDA_prob_", colnames(SimulateSeuratObj@meta.data))))
  testthat::expect_equal(5, sum(grepl("^CD4_sPLSDA_prob_", colnames(SimulateSeuratObj@meta.data))))
  testthat::expect_equal(14, sum(grepl("^GeneralizedTCR_sPLSDA_prob_", colnames(SimulateSeuratObj@meta.data))))

  #classes
  testthat::expect_equal(1, sum(grepl("CD8_sPLSDA_class", colnames(SimulateSeuratObj@meta.data))))
  testthat::expect_equal(1, sum(grepl("CD4_sPLSDA_class", colnames(SimulateSeuratObj@meta.data))))
  testthat::expect_equal(1, sum(grepl("GeneralizedTCR_sPLSDA_class", colnames(SimulateSeuratObj@meta.data))))

  #verify specific probability value
  testthat::expect_equal(SimulateSeuratObj@meta.data[5, 'CD4_sPLSDA_prob_Resting1'], expected = 0.0086, tolerance = 0.001)
  
  table(SimulateSeuratObj$GeneralizedTCR_sPLSDA_ConsensusClass, SimulateSeuratObj$GeneralizedTCR_sPLSDA_class)
  testthat::expect_equal(198, sum(SimulateSeuratObj$GeneralizedTCR_sPLSDA_ConsensusClass == 'Cultured T cell'))
  testthat::expect_equal(2, sum(SimulateSeuratObj$GeneralizedTCR_sPLSDA_ConsensusClass == 'Bulk Tissue T cell'))
  testthat::expect_equal(73, sum(SimulateSeuratObj$GeneralizedTCR_sPLSDA_class == 'Cultured_Bystander_BFA'))
})


test_that("CombineTcellActivationClasses combines classes and writes versioned columns", {
  seuratObj <- SimulateSeuratData()
  seuratObj <- PredictTcellActivation(seuratObj, combineClasses = FALSE)

  classMapping <- list(
    "Th1" = c("Th1_MIP1B.neg_CD137.neg", "Th1_MIP1B.neg_CD137.pos", "Th1_MIP1B.pos"),
    "Th17" = c("Th17"), 
    "Bystander Activated" = c("NonSpecificActivated_L1", "NonSpecificActivated_L2"), 
    "Bulk Tissue T cell" = c("Uncultured"), 
    "Naive T cell" = c("Cultured_Bystander_NoBFA", "Cultured_Bystander_BFA")
  )

  seuratObj <- .CombineTcellActivationClasses(seuratObj, classMapping = classMapping, sourceFieldPrefix = 'GeneralizedTCR_sPLSDA', outputFieldPrefix = 'CustomizedClasses')

  #confirm combined class column exists with model/version naming
  testthat::expect_true("CustomizedClasses_ConsensusClass" %in% colnames(seuratObj@meta.data))

  #confirm combined probability columns exist for each mapping
  for (newClassName in names(classMapping)) {
    testthat::expect_true(paste0("CustomizedClasses_prob_", newClassName) %in% colnames(seuratObj@meta.data))
  }

  #test numeric values
  testthat::expect_type(seuratObj@meta.data[1, "CustomizedClasses_ConsensusClass"], "character")
  
  #validate error cases
  
  #empty names not allowed
  badMapping1 <- list(NULL = c("Th17"))
  expect_error(.CombineTcellActivationClasses(seuratObj, classMapping = badMapping1, sourceFieldPrefix = 'GeneralizedTCR_sPLSDA', outputFieldPrefix = 'CustomizedClasses'))
  #non-character values not allowed
  badMapping2 <- list("Th17" = 1)
  expect_error(.CombineTcellActivationClasses(seuratObj, classMapping = badMapping2, sourceFieldPrefix = 'GeneralizedTCR_sPLSDA', outputFieldPrefix = 'CustomizedClasses'))
  #empty vector not allowed
  badMapping3 <- list("Th17" = character(0))
  expect_error(.CombineTcellActivationClasses(seuratObj, classMapping = badMapping3, sourceFieldPrefix = 'GeneralizedTCR_sPLSDA', outputFieldPrefix = 'CustomizedClasses'))
  #duplicate classes across mappings not allowed
  badMapping4 <- list("A" = c("Th17"), "B" = c("Th17"))
  expect_error(.CombineTcellActivationClasses(seuratObj, classMapping = badMapping4, sourceFieldPrefix = 'GeneralizedTCR_sPLSDA', outputFieldPrefix = 'CustomizedClasses'))
  #classes not present in predictions should error
  badMapping5 <- list("X" = c("NotAClass"))
  expect_warning(.CombineTcellActivationClasses(seuratObj, classMapping = badMapping5, sourceFieldPrefix = 'GeneralizedTCR_sPLSDA', outputFieldPrefix = 'CustomizedClasses'))
})


test_that("Activation class mapping registry works and validates structure", {
  #known keys should exist and return a named list
  basic <- GetActivationClassMapping('TcellActivation.Basic')
  expect_true(is.list(basic))
  expect_true(length(basic) > 0)
  expect_true(!is.null(names(basic)))
  expect_false(any(is.na(names(basic)) | names(basic) == ""))

  #values should be non-empty character vectors with non-empty entries
  for (nm in names(basic)) {
    v <- basic[[nm]]
    expect_true(is.character(v), label = paste0("Value for '", nm, "' should be character"))
    expect_true(length(v) > 0, label = paste0("Value for '", nm, "' should not be empty"))
    expect_false(any(is.na(v) | v == ""), label = paste0("Value for '", nm, "' should not contain empty/NA"))
  }

  # Unknown key should warn and return NULL
  expect_error(x <- GetActivationClassMapping('This.Key.Does.Not.Exist'))
})

