#' @include Utils.R

utils::globalVariables(
  names = c('majority_voting', 'Fraction', 'PropPerCluster', 'over_clustering', 'predicted_labels', 'totalPerCluster', 'totalPerLabel', 'propPerLabel', 'sortOrder'),
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
#' @param pThreshold By default, this would be passed to the --p-thres argument. However, if you also provide extraArgs, this is ignored.
#' @param minProp By default, this would be passed to the --min-prop argument. However, if you also provide extraArgs, this is ignored.
#' @param useMajorityVoting If true, the celltypist --majority-voting option will be added.
#' @param mode Passed directly to the celltypist --mode argument. Either best_match or prob_match
#' @param maxAllowableClasses Celltypist can assign a cell to many classes, creating extremely long labels. Any cell with more than this number of labels will be set to NA
#' @param minFractionToInclude If non-null, any labels with fewer than this fraction of cells will be set to NA.
#' @param minCellsToRun If the input seurat object has fewer than this many cells, NAs will be added for all expected columns and celltypist will not be run.
#' @param maxBatchSize If more than this many cells are in the object, it will be split into batches of this size and run in serial.
#' @param retainProbabilityMatrix If true, the celltypist probability_matrix with per-class probabilities will be stored in meta.data
#' @param runCelltypistUpdate If true, --update-models will be run for celltypist prior to scoring cells.
#'
#' @importFrom dplyr %>%
#' @export
RunCellTypist <- function(seuratObj, modelName = "Immune_All_Low.pkl", pThreshold = 0.5, minProp = 0, useMajorityVoting = TRUE, mode = "prob_match", extraArgs = c("--mode", mode, "--p-thres", pThreshold, "--min-prop", minProp), assayName = Seurat::DefaultAssay(seuratObj), columnPrefix = NULL, maxAllowableClasses = 6, minFractionToInclude = 0.01, minCellsToRun = 200, maxBatchSize = 100000, retainProbabilityMatrix = FALSE, runCelltypistUpdate = TRUE) {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed! If you believe it has been installed, run reticulate::import("celltypist") to get more information and debug')
  }

  if (ncol(seuratObj) < minCellsToRun) {
    warning('Too few cells, will not run celltypist. NAs will be added instead')
    expectedCols <- c('predicted_labels', 'over_clustering', 'majority_voting', 'cellclass')
    if (!is.null(columnPrefix)) {
      expectedCols <- paste0(columnPrefix, expectedCols)
    }

    for (colName in expectedCols) {
      seuratObj[[colName]] <- NA
    }

    return(seuratObj)
  }

  if (!.HasNormalizationBeenPerformed(seuratObj, assayName)) {
    print('Seurat object does not appear to be normalized, running now:')
    seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE, assay = assayName)
  }

  print(paste0('Running celltypist using model: ', modelName))

  shouldDownloadModels <- runCelltypistUpdate
  mf <- system.file(paste0("models/", modelName, '.pkl'), package = "RIRA")
  if (file.exists(mf)) {
    print(paste0('Using RIRA model: ', modelName))
    modelName <- mf
    shouldDownloadModels <- FALSE
  }

  nBatches <- 1
  if (ncol(seuratObj) > maxBatchSize) {
    nBatches <- ceiling(ncol(seuratObj) / maxBatchSize)
    print(paste0('The object will be split into ', nBatches, ', batches'))
  }

  labels <- NULL
  if (nBatches == 1) {
    labels <- .RunCelltypistOnSubset(seuratObj = seuratObj, assayName = assayName, modelName = modelName, useMajorityVoting = useMajorityVoting, extraArgs = extraArgs, updateModels = shouldDownloadModels, retainProbabilityMatrix = retainProbabilityMatrix)
  }
  else {
    cellsPerBatch <- .SplitCellsIntoBatches(seuratObj, nBatches = nBatches)
    for (i in 1:nBatches) {
      toRun <- cellsPerBatch[[i]]
      print(paste0('Running batch ', i, ' of ', nBatches, ' with ', length(toRun), ' cells'))
      so <- subset(seuratObj, cells = toRun)
      if (ncol(so) != length(toRun)) {
        stop(paste0('Error subsetting seurat object, batch size does not match cells after subset: ', length(toRun), ' / ', ncol(seuratObj)))
      }

      if (i > 1) {
        shouldDownloadModels <- FALSE
      }

      df <- .RunCelltypistOnSubset(seuratObj = so, assayName = assayName, modelName = modelName, useMajorityVoting = useMajorityVoting, extraArgs = extraArgs, updateModels = shouldDownloadModels, retainProbabilityMatrix = retainProbabilityMatrix)
      if (useMajorityVoting) {
        df$over_clustering <- paste0(i, '-', df$over_clustering)
      }
      rm(so)

      if (all(is.null(labels))) {
        labels <- df
      } else {
        labels <- rbind(labels, df)
      }
    }

    if (nrow(labels) != ncol(seuratObj)) {
      stop('There was an error processing celltypist batches')
    }
  }

  labels$predicted_labels[is.na(labels$predicted_labels)] <- 'Unassigned'
  if (useMajorityVoting) {
    labels$majority_voting[is.na(labels$majority_voting)] <- 'Unassigned'
  }

  if (useMajorityVoting) {
      print(ggplot(data.frame(table(labels$over_clustering)), aes(x = Freq)) +
            geom_histogram() +
            labs(x = 'Cluster Size', y = '# Clusters') +
            egg::theme_presentation(base_size = 12) +
            ggtitle('Distribution of Cluster Size')
      )

      dat <- labels
      dat$Category <- ifelse(!is.na(dat$majority_voting) & dat$majority_voting == 'Heterogeneous', yes = 'Heterogenous', no = 'Not Heterogenous')
      dat <- dat %>% dplyr::group_by(Category, over_clustering) %>% dplyr::mutate(totalPerCluster = dplyr::n())
      dat <- dat %>% dplyr::group_by(Category, over_clustering, totalPerCluster, predicted_labels) %>% dplyr::summarize(totalPerLabel = dplyr::n())
      dat$propPerLabel <- dat$totalPerLabel / dat$totalPerCluster
      dat <- dat %>% dplyr::group_by(Category, over_clustering) %>% dplyr::summarize(PropPerCluster = max(propPerLabel))
      P1 <- ggplot(dat, aes(x = PropPerCluster)) +
              geom_density() +
              labs(x = 'Max Prop Per Cluster', y = '# Clusters') +
              egg::theme_presentation(base_size = 12) +
              ggtitle('Proportion of Highest Class Per Cluster') +
              facet_grid(Category ~ .)

      if (!is.null(minProp) && minProp > 0) {
        P1 <- P1 + geom_vline(xintercept = minProp, linetype = 'dashed', color = 'red')
      }

      print(P1)
  }

  if (!is.na(maxAllowableClasses)) {
    for (fieldName in c('majority_voting', 'predicted_labels')) {
      if (!fieldName %in% names(labels)) {
        next
      }

      dat <- unique(labels[[fieldName]])
      toDrop <- dat[lengths(regmatches(x = dat, gregexpr("\\|", dat))) > maxAllowableClasses]
      toDrop <- toDrop[!is.na(toDrop)]
      if (length(toDrop) > 0) {
        print(paste0('Dropping cells with more than ', maxAllowableClasses, ' calls for: ', fieldName, '. These were:'))
        print(paste0(toDrop, collapse = ', '))
      }

      # NOTE: %in% doesnt handle NAs well
      labels[[fieldName]][!is.na(labels[[fieldName]]) & labels[[fieldName]] %in% toDrop] <- 'Ambiguous'
    }
  }

  plotColname <- ifelse('majority_voting' %in% names(labels), yes = 'majority_voting', no = 'predicted_labels')
  if (!is.null(columnPrefix)) {
    names(labels) <- paste0(columnPrefix, names(labels))
    plotColname <- paste0(columnPrefix, plotColname)
  }

  # Create a single output column that will be consistent no matter whether majority_voting used or not.
  # This will simplify some of the original values, but the source columns are also preserved:
  consensusColName <- paste0(columnPrefix, 'cellclass')
  labels[[consensusColName]] <- as.character(labels[[plotColname]])
  labels[[consensusColName]][grepl(labels[[consensusColName]], pattern = '\\|')] <- 'Ambiguous'
  if ('Heterogenous' %in% labels[[consensusColName]]) {
    labels[[consensusColName]][labels[[consensusColName]] == 'Heterogeneous'] <- 'Ambiguous'
  }
  labels[[consensusColName]] <- naturalsort::naturalfactor(labels[[consensusColName]])
  plotColname <- consensusColName

  seuratObj <- Seurat::AddMetaData(seuratObj, labels)

  seuratObj <- .FilterLowCalls(seuratObj, plotColname, minFractionToInclude)
  if (length(seuratObj@reductions) == 0) {
    print('No reductions calculated, cannot plot tSNE/UMAP')
  }
  else if (sum(!is.na(seuratObj[[plotColname]])) == 0) {
    print('No values returned, skipping plot')
  } else {
    print(DimPlot(seuratObj, group.by = plotColname, shuffle = TRUE))
  }

  toPlot <- seuratObj[[plotColname]]
  if (sum(!is.na(toPlot)) != 0) {
    print(ggplot(toPlot, aes(x = !!rlang::sym(plotColname), fill = !!rlang::sym(plotColname))) +
            geom_bar(color = 'black') +
            egg::theme_presentation(base_size = 12) +
            ggtitle('Celltypist Call') +
            labs(x = 'Celltypist Call', y = '# Cells') +
            theme(
              legend.position = 'none',
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            )
    )
  }
  return(seuratObj)
}

.RunCelltypistOnSubset <- function(seuratObj, assayName, modelName, useMajorityVoting, extraArgs, retainProbabilityMatrix, updateModels = TRUE) {
  outDir <- tempfile(fileext = '')
  matrixFile <- SeuratToMatrix(seuratObj, outDir = outDir, assayName = assayName)

  # Ensure models present:
  if (updateModels) {
    system2(reticulate::py_exe(), c("-m", "celltypist.command_line", "--update-models", "--quiet"))
  }

  # Now run celltypist itself:
  args <- c("-m", "celltypist.command_line", "-i", matrixFile, "-m", modelName, "--outdir", outDir, "--prefix", "celltypist.", "--quiet")

  # NOTE: this produces a series of PDFs, one per class. Consider either providing an argument on where to move these, or reading/printing them
  #if (generatePlots) {
  #  args <- c(args, "--plot-results")
  #}

  if (useMajorityVoting) {
    args <- c(args, "--majority-voting")
  }

  if (!all(is.na(extraArgs))) {
    args <- c(args, extraArgs)
  }

  system2(reticulate::py_exe(), args)

  labelFile <- paste0(outDir, '/celltypist.predicted_labels.csv')
  if (!file.exists(labelFile)) {
    stop(paste0('Missing file: ', labelFile, '. files present: ', paste0(list.files(outDir), collapse = ', ')))
  }

  labels <- utils::read.csv(labelFile, header = T, row.names = 1, stringsAsFactors = FALSE)

  if ('majority_voting' %in% names(labels)) {
    labels$majority_voting[labels$majority_voting == 'Unassigned'] <- NA
  }

  probabilityMatrixFile <- paste0(outDir, '/celltypist.probability_matrix.csv')
  if (retainProbabilityMatrix) {
    if (!file.exists(probabilityMatrixFile)) {
      stop(paste0('Missing file: ', probabilityMatrixFile, '. files present: ', paste0(list.files(outDir), collapse = ', ')))
    }

    probabilityMatrix <- utils::read.csv(probabilityMatrixFile, header = T, row.names = 1, check.names = T, stringsAsFactors = FALSE)
    names(probabilityMatrix) <- paste0('prob.', names(probabilityMatrix))
    labels <- cbind(labels, probabilityMatrix)
  }

  unlink(outDir, recursive = TRUE)
  unlink(labelFile)
  unlink(probabilityMatrixFile)
  unlink(paste0(outDir, '/celltypist.decision_matrix.csv'))

  return(labels)
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
#' @param featureInclusionList If provided, the input count matrix will be subset to just these features. If used, Seurat::NormalizeData will be re-run.
#' @param featureExclusionList If provided, the input count matrix will be subset to remove these features. If used, Seurat::NormalizeData will be re-run.
#'
#' @export
TrainCellTypist <- function(seuratObj, labelField, modelFile, minCellsPerClass = 20, assayName = Seurat::DefaultAssay(seuratObj), tempFileLocation = NULL, dropAmbiguousLabelValues = TRUE, excludedClasses = NULL, featureInclusionList = NULL, featureExclusionList = NULL) {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed! If you believe it has been installed, run reticulate::import("celltypist") to get more information and debug')
  }

  if (!labelField %in% names(seuratObj@meta.data)) {
    stop(paste0('The labelField must be a field in seuratObj@meta.data. not found: ', labelField))
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

  print(paste0('Initial cells: ', ncol(seuratObj)))
  if (!is.null(excludedClasses)) {
    for (label in excludedClasses) {
      print(paste0('Dropping label: ', label))
      toKeep <- rownames(seuratObj@meta.data)[seuratObj@meta.data[[labelField]] != label]
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

  print('Summary of cells per class:')
  print(sort(table(seuratObj@meta.data[[labelField]])))

  shouldNormalize <- FALSE
  if (!all(is.null(featureInclusionList))) {
    ad <- seuratObj@assays[[assayName]]
    featureInclusionList <- RIRA::ExpandGeneList(featureInclusionList)
    preExisting <- intersect(rownames(ad), featureInclusionList)
    print(paste0('Limiting to ', length(featureInclusionList), ' features, of which ', length(preExisting), ' exist in this assay'))
    ad <- subset(ad, features = preExisting)
    print(paste0('Total features after: ', nrow(ad)))
    seuratObj@assays[[assayName]] <- ad
    shouldNormalize <- TRUE
  }

  if (!all(is.null(featureExclusionList))){
    ad <- seuratObj@assays[[assayName]]
    featureExclusionList <- RIRA::ExpandGeneList(featureExclusionList)
    toDrop <- intersect(rownames(ad), featureExclusionList)
    print(paste0('Excluding ', length(featureExclusionList), ' features(s) from the input assay, of which ', length(toDrop), ' exist in this assay'))
    if (length(toDrop) == 0) {
      stop(paste0('None of the featureExclusionList features were found in this object: ', paste0(featureExclusionList, collapse = ',')))
    }

    featuresToKeep <- rownames(ad)[!rownames(ad) %in% toDrop]
    ad <- subset(ad, features = featuresToKeep)
    print(paste0('Total features after: ', nrow(ad)))
    seuratObj@assays[[assayName]] <- ad
    shouldNormalize <- TRUE
  }

  if (shouldNormalize) {
    seuratObj <- Seurat::NormalizeData(seuratObj, verbose = FALSE)
  }

  trainData <- SeuratToAnnData(seuratObj, paste0(outFile, '-seurat-annData'), assayName = assayName, exportMinimalObject = TRUE, allowableMetaCols = labelField)

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

#' @title Classify T/NK
#'
#' @description Runs celltypist using the RIRA T/NK 4-class model to score cells using CellTypist with optimized parameters.
#' @param seuratObj The seurat object
#' @param assayName The name of the assay to use. Others will be dropped
#' @param columnPrefix A prefix that will be added to the beginning of the resulting columns, added the seurat@meta.data
#' @param maxAllowableClasses Celltypist can assign a cell to many classes, creating extremely long labels. Any cell with more than this number of labels will be set to NA
#' @param minFractionToInclude If non-null, any labels with fewer than this fraction of cells will be set to NA.
#' @param minCellsToRun If the input seurat object has fewer than this many cells, NAs will be added for all expected columns and celltypist will not be run.
#' @param maxBatchSize If more than this many cells are in the object, it will be split into batches of this size and run in serial.
#' @param retainProbabilityMatrix If true, the celltypist probability_matrix with per-class probabilities will be stored in meta.data
#'
#' @export
Classify_TNK <- function(seuratObj, assayName = Seurat::DefaultAssay(seuratObj), columnPrefix = 'RIRA_TNK_v2.', maxAllowableClasses = 6, minFractionToInclude = 0.01, minCellsToRun = 200, maxBatchSize = 600000, retainProbabilityMatrix = FALSE) {
  return(RunCellTypist(seuratObj = seuratObj,
         modelName = "RIRA_TNK_v2",
         # These are optimized for this model:
         pThreshold = 0.5, minProp = 0, useMajorityVoting = FALSE, mode = "prob_match",

         assayName = assayName,
         columnPrefix = columnPrefix,
         maxAllowableClasses = maxAllowableClasses,
         minFractionToInclude = minFractionToInclude,
         minCellsToRun = minCellsToRun,
         maxBatchSize = maxBatchSize,
         retainProbabilityMatrix = retainProbabilityMatrix
  ))
}


#' @title Classify Bulk Immune cells
#'
#' @description Runs celltypist using the RIRA bulk immune model to score cells using CellTypist with optimized parameters.
#' @param seuratObj The seurat object
#' @param assayName The name of the assay to use. Others will be dropped
#' @param columnPrefix A prefix that will be added to the beginning of the resulting columns, added the seurat@meta.data
#' @param maxAllowableClasses Celltypist can assign a cell to many classes, creating extremely long labels. Any cell with more than this number of labels will be set to NA
#' @param minFractionToInclude If non-null, any labels with fewer than this fraction of cells will be set to NA.
#' @param minCellsToRun If the input seurat object has fewer than this many cells, NAs will be added for all expected columns and celltypist will not be run.
#' @param maxBatchSize If more than this many cells are in the object, it will be split into batches of this size and run in serial.
#' @param retainProbabilityMatrix If true, the celltypist probability_matrix with per-class probabilities will be stored in meta.data
#' @param filterDisallowedClasses If true, this will run FilterDisallowedClasses() on the output.
#'
#' @export
Classify_ImmuneCells <- function(seuratObj, assayName = Seurat::DefaultAssay(seuratObj), columnPrefix = 'RIRA_Immune_v2.', maxAllowableClasses = 6, minFractionToInclude = 0.01, minCellsToRun = 200, maxBatchSize = 600000, retainProbabilityMatrix = FALSE, filterDisallowedClasses = TRUE) {
  seuratObj <- RunCellTypist(seuratObj = seuratObj,
                       modelName = 'RIRA_Immune_v2',

                       # These are optimized for this model:
                       minProp = 0.5, useMajorityVoting = TRUE, mode = "prob_match",

                       assayName = assayName,
                       columnPrefix = columnPrefix,
                       maxAllowableClasses = maxAllowableClasses,
                       minFractionToInclude = minFractionToInclude,
                       minCellsToRun = minCellsToRun,
                       maxBatchSize = maxBatchSize,
                       retainProbabilityMatrix = retainProbabilityMatrix
  )

  # Create a simplified final column:
  targetField <- paste0(columnPrefix, 'cellclass')
  if (!targetField %in% names(seuratObj@meta.data)) {
    foundFields <- sort(unique(c(names(seuratObj@meta.data)[grepl(names(seuratObj@meta.data), pattern = 'cellclass')], names(seuratObj@meta.data)[grepl(names(seuratObj@meta.data), pattern = columnPrefix)])))
    stop(paste0('Missing expected field, something went wrong with celltypist: ', targetField, '. Fields found: ', paste0(foundFields, collapse = ',')))
  }

  seuratObj@meta.data[[targetField]] <- as.character(seuratObj@meta.data[[targetField]])
  seuratObj@meta.data[[targetField]][!is.na(seuratObj@meta.data[[targetField]]) & grepl(seuratObj@meta.data[[targetField]], pattern = '\\|')] <- 'Ambiguous'
  if (filterDisallowedClasses) {
    seuratObj <- FilterDisallowedClasses(seuratObj, sourceField = targetField)
    seuratObj@meta.data[[targetField]][!is.na(seuratObj@meta.data$DisallowedUCellCombinations)] <- 'Unknown'
  }

  seuratObj@meta.data[[targetField]][!is.na(seuratObj@meta.data[[targetField]]) & seuratObj@meta.data[[targetField]] %in% c('AvEp', 'Epithelial', 'Stromal', 'Mesothelial', 'ActiveAvEp', 'Myofibroblast', 'Fibroblast', 'Hepatocyte')] <- 'Non-Immune'
  seuratObj@meta.data[[targetField]][!is.na(seuratObj@meta.data[[targetField]]) & seuratObj@meta.data[[targetField]] %in% c('Unassigned', 'Contamination', 'Ambiguous', 'Heterogeneous', 'Unknown')] <- 'Unknown'
  seuratObj@meta.data[[targetField]] <- naturalsort::naturalfactor(seuratObj@meta.data[[targetField]])

  print(ggplot(seuratObj@meta.data, aes(x = !!rlang::sym(targetField), fill = !!rlang::sym(targetField))) +
      geom_bar(color = 'black') +
      egg::theme_presentation(base_size = 12) +
      ggtitle('Celltypist Call After Rename') +
      labs(x = 'Celltypist Call', y = '# Cells') +
      theme(
        legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  )

  return(seuratObj)
}

#' @title Filter Disallowed Classes
#'
#' @description This is used to flag cells with high UCell score combinations generally indiciative of contamination or doublets
#' @param seuratObj The seurat object
#' @param sourceField The name of the field on which to compare
#' @param outputFieldName The name of the field to store the results
#' @param ucellCutoff Any cells expressing the disallowed UCell above this value will be flagged
#' @param disallowedClasses This is a list where the names are the cell classes (which should match levels in sourceField), and values are a vector of UCell field names.
#'
#' @export
FilterDisallowedClasses <- function(seuratObj, sourceField = 'RIRA_Immune_v2.majority_voting', outputFieldName = 'DisallowedUCellCombinations', ucellCutoff = 0.2, disallowedClasses = list(
  T_NK = c('Bcell.RM_UCell', 'Myeloid.RM_UCell', 'Erythrocyte.RM_UCell', 'Platelet.RM_UCell', 'NeutrophilLineage.RM_UCell'),
  Myeloid = c('Bcell.RM_UCell', 'Tcell.RM_UCell', 'NK.RM_UCell', 'Erythrocyte.RM_UCell', 'Platelet.RM_UCell'),
  Bcell = c('Tcell.RM_UCell', 'NK.RM_UCell', 'Myeloid.RM_UCell', 'Erythrocyte.RM_UCell', 'Platelet.RM_UCell', 'NeutrophilLineage.RM_UCell', 'Complement.RM_UCell')
)) {
  if (!sourceField %in% names(seuratObj@meta.data)) {
    stop(paste0('Missing source field: ', sourceField))
  }

  toDrop <- data.frame(cellbarcode = character(), reason = character())
  for (cls in names(disallowedClasses)) {
    for (ucell in disallowedClasses[[cls]]) {
      if (!ucell %in% names(seuratObj@meta.data) || any(is.na(seuratObj@meta.data[[ucell]]))) {
        modelName <- gsub(ucell, pattern = '_UCell', replacement = '')
        model <- GetScGateModel(modelName)
        if (!is.null(model)) {
          seuratObj <- RunScGate(seuratObj, modelName)
        } else {
          stop(paste0('Missing UCell field: ', ucell))
        }
      }

      x <- colnames(seuratObj)[seuratObj@meta.data[[sourceField]] == cls & seuratObj@meta.data[[ucell]] > ucellCutoff]
      if (length(x) > 0) {
        toDrop <- rbind(toDrop, data.frame(cellbarcode = x, reason = ucell))
      }
    }
  }

  allCells <- data.frame(cellbarcode = colnames(seuratObj), sortOrder = 1:ncol(seuratObj))
  if (nrow(toDrop) > 0) {
    allCells <- merge(allCells, toDrop, by = 'cellbarcode', all.x = T)
    allCells <- dplyr::arrange(allCells, sortOrder)

    toAdd <- allCells$reason
    names(toAdd) <- allCells$cellbarcode
    seuratObj <- Seurat::AddMetaData(seuratObj, toAdd, col.name = outputFieldName)
  } else {
    seuratObj[[outputFieldName]] <- NA
  }

  print(table(seuratObj@meta.data[[sourceField]], seuratObj@meta.data[[outputFieldName]]))

  return(seuratObj)
}