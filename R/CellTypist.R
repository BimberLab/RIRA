#' @include Utils.R


#' @title Run Celltypist
#'
#' @description Runs celltypist on a seurat object and stores the calls as metadata
#' @param seuratObj The seurat object
#' @param modelName The model name or path to celltypist model
#'
#' @export
RunCellTypist <- function(seuratObj, modelName = "Immune_All_Low.pkl") {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed!')
  }

  print(paste0('Running celltypist using model: ', modelName))
  outFile <- tempfile()
  outDir <- dirname(outFile)
  seuratObj <- Seurat::DietSeurat(seuratObj)
  seuratAnnData <- SeuratToAnnData(seuratObj, paste0(outFile, '.seurat.annData'))

  #exe <- reticulate::py_exe()
  exe <- "celltypist"
  #"-m", "celltypist.command_line",
  pyOut <- system2(exe, c("--update-models", "--quiet"), stdout = TRUE, stderr = TRUE)
  print(pyOut)

  # "-m", "celltypist.command_line",
  # "--plot-results"
  args <- c("-i", seuratAnnData, "-m", modelName, "--outdir", outDir, "--majority-voting", "--prefix", "celltypist.", "--quiet")
  pyOut <- system2(exe, args, stdout = TRUE, stderr = TRUE)
  print(pyOut)

  labels <- paste0(outDir, '/celltypist.predicted_labels.csv')
  if (!file.exists(labels)) {
    stop(paste0('Missing file: ', labels))
  }

  labels <- utils::read.csv(labels, header = T, row.names = 1)
  seuratObj <- Seurat::AddMetaData(seuratObj, labels)

  unlink(seuratAnnData)
  unlink(labels)

  return(seuratObj)
}

#' @title Train Celltypist
#'
#' @description Train celltypist using the data in a seurat object and saves the resulting model
#' @param seuratObj The seurat object
#' @param labelField The field in seuratObj@meta.data holding the labels for training
#' @param modelFile The path to save the model
#'
#' @export
TrainCellTypist <- function(seuratObj, labelField, modelFile) {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed!')
  }

  if (!labelField %in% names(seuratObj@meta.data)) {
    stop('The labelField must be a field in seuratObj@meta.data')
  }

  modelFile <- gsub(modelFile, pattern = '\\\\', replacement = '/')
  outFile <- gsub(tempfile(), pattern = '\\\\', replacement = '/')
  outDir <- dirname(outFile)
  seuratObj <- Seurat::DietSeurat(seuratObj)
  trainData <- SeuratToAnnData(seuratObj, paste0(outFile, '-seurat-annData'))
  labelFile <- paste0(outDir, '.seurat.labels.txt')
  scriptFile <- paste0(outDir, '.seurat.train.py')

  write.table(seuratObj@meta.data[[labelField]], row.names = F, sep = '\t', quote = F, col.names = F, file = labelFile)

  typistCommand <- c(
    "import celltypist;",
    paste0("new_model = celltypist.train('", trainData, "', labels = '", labelFile, "', use_SGD = False, solver = 'saga', feature_selection = True, top_genes = 300);"),
    paste0("new_model.write('", modelFile, "');")
  )

  write(typistCommand, file = scriptFile)

  pyOut <- system2(reticulate::py_exe(), c(scriptFile), stdout = TRUE, stderr = TRUE)
  print(pyOut)

  #unlink(labelFile)
  #unlink(scriptFile)
  #unlink(trainData)
}

