#' @include Utils.R

utils::globalVariables(
  names = c('majority_voting'),
  package = 'RIRA',
  add = TRUE
)

#' @title Run Celltypist
#'
#' @description Runs celltypist on a seurat object and stores the calls as metadata
#' @param seuratObj The seurat object
#' @param modelName The model name or path to celltypist model
#' @param extraArgs An optional list of additional arguments passed directly on the command line to cell typist
#' @param assayName The name of the assay to use. Others will be dropped
#' @param columnPrefix A prefix that will be added to the beginning of the resulting columns, added the seurat@meta.data
#' @param convertAmbiguousToNA If true, any values for majority_voting with commas (indicating they are ambiguous) will be converted to NA
#' @param pThreshold By default, this would be passed to the --p-thres argument. However, if you also provide extraArgs, this is ignored.
#' @param minProp By default, this would be passed to the --min-prop argument. However, if you also provide extraArgs, this is ignored.
#'
#' @export
RunCellTypist <- function(seuratObj, modelName = "Immune_All_Low.pkl", pThreshold = 0.5, minProp = 0, extraArgs = c("--majority-voting", "--mode", "prob_match", "--p-thres", pThreshold, "--min-prop", minProp), assayName = 'RNA', columnPrefix = NULL, convertAmbiguousToNA = FALSE) {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed! If you believe it has been installed, run reticulate::import("celltypist") to get more information and debug')
  }

  if (!.HasNormalizationBeenPerformed(seuratObj, assayName)) {
    print('Seurat object does not appear to be normalized, running now:')
    seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  }

  print(paste0('Running celltypist using model: ', modelName))
  outFile <- tempfile()
  outDir <- dirname(outFile)
  seuratAnnData <- SeuratToAnnData(seuratObj, paste0(outFile, '-seurat-annData'), assayName, doDiet = TRUE)

  # Ensure models present:
  system2(reticulate::py_exe(), c("-m", "celltypist.command_line", "--update-models", "--quiet"))

  # Now run celltypist itself:
  args <- c("-m", "celltypist.command_line", "-i", seuratAnnData, "-m", modelName, "--outdir", outDir, "--prefix", "celltypist.", "--quiet")

  # NOTE: this produces a series of PDFs, one per class. Consider either providing an argument on where to move these, or reading/printing them
  #if (generatePlots) {
  #  args <- c(args, "--plot-results")
  #}

  if (!all(is.na(extraArgs))) {
    args <- c(args, extraArgs)
  }

  system2(reticulate::py_exe(), args)

  labels <- paste0(outDir, '/celltypist.predicted_labels.csv')
  if (!file.exists(labels)) {
    stop(paste0('Missing file: ', labels))
  }

  labels <- utils::read.csv(labels, header = T, row.names = 1)
  labels$majority_voting[labels$majority_voting == 'Unassigned'] <- NA

  if (convertAmbiguousToNA) {
    toDrop <- grepl(labels$majority_voting, pattern = ',')
    if (sum(toDrop) > 0) {
      print(paste0('Converting ', sum(toDrop), ' cells with ambiguous values for majority_voting to NAs'))
      labels$majority_voting[toDrop] <- NA
    }
  }

  if (!is.null(columnPrefix)) {
    names(labels) <- paste0(columnPrefix, names(labels))
  }

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
#' @param excludedClasses A vector of labels to discard.
#' @param tempFileLocation The location where temporary files (like the annData version of the seurat object), will be written.
#' @param dropAmbiguousLabelValues If true, and label value with a comma will be dropped.
#'
#' @export
TrainCellTypist <- function(seuratObj, labelField, modelFile, minCellsPerClass = 20, assayName = 'RNA', tempFileLocation = NULL, dropAmbiguousLabelValues = TRUE, excludedClasses = NULL) {
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
  if (is.null(tempFileLocation)) {
    outDir <- tempdir()
  } else {
    outDir <- tempFileLocation
  }

  if (endsWith(outDir, "/")){
    outDir <- gsub(outDir, pattern = "/$", replacement = "")
  }
  outFile <- tempfile(tmpdir = outDir)

  if (!is.null(excludedClasses)) {
    for (label in excludedClasses) {
      print(paste0('Dropping label: ', label))
        toKeep <- rownames(seuratObj@meta.data)[seuratObj@meta.data[[label]] != label]
        seuratObj <- subset(seuratObj, cells = toKeep)
      print(paste0('Cells remaining: ', ncol(seuratObj)))
    }
  }

  if (sum(is.null(seuratObj@meta.data[[labelField]]) | is.na(seuratObj@meta.data[[labelField]])) > 0) {
    initialCells <- ncol(seuratObj@assays[[assayName]])
    toKeep <- colnames(seuratObj@assays[[assayName]])[!is.null(seuratObj@meta.data[[labelField]]) & !is.na(seuratObj@meta.data[[labelField]])]
    print(paste0('Dropping cells with NA or NULL labels, total dropped: ', (initialCells - length(toKeep))))
    seuratObj <- subset(seuratObj, cells = toKeep)
  }

  if (dropAmbiguousLabelValues) {
    toDrop <- grepl(seuratObj@meta.data[[labelField]], pattern = ',')
    if (sum(toDrop) > 0) {
      print('Dropping the following ambiguous consensus labels:')
      print(sort(table(seuratObj@meta.data[[labelField]][toDrop])))

      initialCells <- ncol(seuratObj@assays[[assayName]])
      toKeep <- colnames(seuratObj@assays[[assayName]])[!toDrop]
      print(paste0('Total cells dropped: ', (initialCells - length(toKeep))))
      seuratObj <- subset(seuratObj, cells = toKeep)
    }
  }

  # Perform this after other filters are applied, so low-count classes are not binned into 'Other'
  if (!is.null(minCellsPerClass) && minCellsPerClass > 0) {
    seuratObj <- .DropLowCountClasses(seuratObj, labelField, minCellsPerClass)
  }

  trainData <- SeuratToAnnData(seuratObj, paste0(outFile, '-seurat-annData'), assayName = assayName, doDiet = TRUE)

  # potentially replace windows slashes with forward slash
  trainData <- gsub(trainData, pattern = '\\\\', replacement = '/')
  outFile <- gsub(outFile, pattern = '\\\\', replacement = '/')

  labelFile <- paste0(outFile, '.seurat.labels.txt')
  scriptFile <- paste0(outFile, '.seurat.train.py')

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

