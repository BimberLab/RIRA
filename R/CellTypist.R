#' @include Utils.R


#' @title Run Celltypist
#'
#' @description Runs celltypist on a seurat object and stores the calls as metadata
#' @param seuratObj The seurat object
#' @param modelName The model name or path to celltypist model
#' @param extraArgs An optional list of additional arguments passed directly on the command line to cell typist
#' @param assayName The name of the assay to use. Others will be dropped
#'
#' @export
RunCellTypist <- function(seuratObj, modelName = "Immune_All_Low.pkl", extraArgs = c("--majority-voting", "--mode", "prob_match", "--p-thres", 0.5), assayName = 'RNA') {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed! If you believe it has been installed, run reticulate::import("celltypist") to get more information and debug')
  }

  if (!.HasNormalizationBeenPerformed(seuratObj, assayName)) {
    print('Seurat object does not appear to be normalized, running now:')
    seuratObj <- NormalizeData(seuratObj, verbose = FALSE)
  }

  print(paste0('Running celltypist using model: ', modelName))
  outFile <- tempfile()
  outDir <- dirname(outFile)
  seuratAnnData <- SeuratToAnnData(seuratObj, paste0(outFile, '-seurat-annData'), assayName, doDiet = TRUE)

  #exe <- reticulate::py_exe()
  exe <- "celltypist"
  #"-m", "celltypist.command_line",
  system2(exe, c("--update-models", "--quiet"))

  # "-m", "celltypist.command_line",
  # "--plot-results"
  args <- c("-i", seuratAnnData, "-m", modelName, "--outdir", outDir, "--prefix", "celltypist.", "--quiet")
  if (!all(is.na(extraArgs))) {
    args <- c(args, extraArgs)
  }

  system2(exe, args)

  labels <- paste0(outDir, '/celltypist.predicted_labels.csv')
  if (!file.exists(labels)) {
    stop(paste0('Missing file: ', labels))
  }

  labels <- utils::read.csv(labels, header = T, row.names = 1)
  labels$majority_voting[labels$majority_voting == 'Unassigned'] <- NA
  seuratObj <- Seurat::AddMetaData(seuratObj, labels)

  unlink(seuratAnnData)
  unlink(labels)

  print(ggplot(seuratObj@meta.data, aes(x = majority_voting, fill = majority_voting)) +
          geom_bar(color = 'black') +
          egg::theme_presentation(base_size = 12) +
          ggtitle('Celltypist Call') +
          labs(x = 'Celltypist Call', y = '# Cells') +
          theme(
            legend.position = 'none',
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
          )
  )

  return(seuratObj)
}

#' @title Train Celltypist
#'
#' @description Train celltypist using the data in a seurat object and saves the resulting model
#' @param seuratObj The seurat object
#' @param labelField The field in seuratObj@meta.data holding the labels for training
#' @param modelFile The path to save the model
#' @param minCellsPerClass If provided, any classes (and corresponding cells) with fewer than this many cells will be dropped from the training data
#' @param assayName The name of the assay to use
#' @param tempFileLocation The location where temporary files (like the annData version of the seurat object), will be written.
#'
#' @export
TrainCellTypist <- function(seuratObj, labelField, modelFile, minCellsPerClass = 20, assayName = 'RNA', tempFileLocation = NULL) {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed! If you believe it has been installed, run reticulate::import("celltypist") to get more information and debug')
  }

  if (!labelField %in% names(seuratObj@meta.data)) {
    stop('The labelField must be a field in seuratObj@meta.data')
  }

  modelFile <- gsub(modelFile, pattern = '\\\\', replacement = '/')
  outDir <- ifelse(is.null(tempFileLocation), yes = tempdir(), no = tempFileLocation)
  if (endsWith(outDir, "/")){
    outDir <- gsub(outDir, pattern = "/$", replacement = "")
  }
  outFile <- tempfile(tmpdir = outDir)

  if (!is.null(minCellsPerClass) && minCellsPerClass > 0) {
    seuratObj <- .DropLowCountClasses(seuratObj, labelField, minCellsPerClass)
  }

  if (sum(is.null(seuratObj@meta.data[[labelField]]) | is.na(seuratObj@meta.data[[labelField]])) > 0) {
    initialCells <- ncol(seuratObj@assays[[assayName]])
    toKeep <- colnames(seuratObj@assays[[assayName]])[!is.null(seuratObj@meta.data[[labelField]]) & !is.na(seuratObj@meta.data[[labelField]])]
    print(paste0('Dropping cells with NA or NULL labels, total dropped: ', (initialCells - length(toKeep))))
    seuratObj <- subset(seuratObj, cells = toKeep)
  }

  trainData <- SeuratToAnnData(seuratObj, paste0(outFile, '-seurat-annData'), assayName = assayName, doDiet = TRUE)

  # potentially replace windows slashes with forward slash
  trainData <- gsub(trainData, pattern = '\\\\', replacement = '/')
  outDir <- gsub(outDir, pattern = '\\\\', replacement = '/')

  labelFile <- paste0(outDir, '.seurat.labels.txt')
  scriptFile <- paste0(outDir, '.seurat.train.py')

  utils::write.table(seuratObj@meta.data[[labelField]], row.names = F, sep = '\t', quote = F, col.names = F, file = labelFile)

  typistCommand <- c(
    "import celltypist;",
    paste0("new_model = celltypist.train('", trainData, "', labels = '", labelFile, "', use_SGD = False, solver = 'saga', feature_selection = True, top_genes = 300);"),
    paste0("new_model.write('", modelFile, "');")
  )

  write(typistCommand, file = scriptFile)

  print('Running celltypist.train:')
  print(typistCommand)

  system2(reticulate::py_exe(), c(scriptFile))

  if (!file.exists(modelFile)) {
    stop(paste0('Unable to find file: ', modelFile))
  }

  unlink(labelFile)
  unlink(scriptFile)
  unlink(trainData)
}

