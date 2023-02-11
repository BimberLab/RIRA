
pkg.env <- new.env(parent=emptyenv())

pkg.env$RANDOM_SEED <- 1234
set.seed(pkg.env$RANDOM_SEED)

#' @title Set random seed
#'
#' @description Sets the seed used for R‘s random number generator, which should be used in all internal functions
#' @param seed The random seed
#' @export
SetSeed <- function(seed) {
  pkg.env$RANDOM_SEED <- seed
  set.seed(pkg.env$RANDOM_SEED)
}

#' @title Get random seed
#'
#' @description Sets a random seed, which should be used in all internal functions
#' @export
GetSeed <- function() {
  return(pkg.env$RANDOM_SEED)
}

#' @title Set the directory where atlas data reside
#'
#' @description Sets the directory where atlas data reside.
#' @param folderPath The atlas directory. Within the folder, there should be sub-folders by version.
#' @export
SetAtlasDir <- function(folderPath) {
  if (!dir.exists(folderPath)) {
    stop(paste0('Unable to find folder: ', folderPath))
  }

  folderPath <- gsub(folderPath, pattern = '[/\\]+$', replacement = '')
  pkg.env$ATLAS_DIR <- folderPath
}

.GetAtlasDir <- function() {
  return(pkg.env$ATLAS_DIR)
}

.GetLatestVersion <- function(){
  return(utils::packageVersion('RIRA'))
}

.GetAtlasBaseDir <- function(version) {
  if (is.null(version)) {
    version <- .GetLatestVersion()
  }

  if (is.null(.GetAtlasDir())) {
    parentFolder <- system.file('inst/data', package = 'RIRA')
  } else {
    parentFolder <- .GetAtlasDir()
  }

  parentFolder <- paste0(parentFolder, '/', version)
  if (!dir.exists(parentFolder)) {
    print(paste0('Unknown version: ', version))
  }

  return(parentFolder)
}

SeuratToAnnData <- function(seuratObj, outFileBaseName, assayName = NULL, doDiet = FALSE, allowableMetaCols = NULL) {
  tmpFile <- outFileBaseName
  if (!is.null(assayName)) {
    for (an in names(seuratObj@assays)) {
      if (an != assayName) {
        seuratObj[[an]] <- NULL
      }
    }
  }

  if (doDiet) {
    seuratObj <- Seurat::DietSeurat(seuratObj)
  }

  seuratObj@misc <- list()
  seuratObj@commands <- list()

  if (all(is.na(allowableMetaCols))) {
    print('Removing all metadata prior to save')
    seuratObj@meta.data <- seuratObj@meta.data[,NULL, drop = FALSE]
    seuratObj@meta.data$CellBarcode <- rownames(seuratObj@meta.data)
  } else if (!all(is.null(allowableMetaCols))) {
    if (!all(allowableMetaCols %in% names(seuratObj@meta.data))) {
      stop('Not all columns requested in allowableMetaCols exist in the seurat object')
    }

    seuratObj@meta.data <- seuratObj@meta.data[,allowableMetaCols, drop = FALSE]
  }

  print('Saving to H5Seurat')
  SeuratDisk::SaveH5Seurat(seuratObj, filename = tmpFile)
  h5seurat <- paste0(tmpFile, ".h5seurat")
  print('Converting from H5Seurat to AnnData')
  SeuratDisk::Convert(source = h5seurat, dest = "h5ad", overwrite = T)
  unlink(h5seurat)

  return(paste0(outFileBaseName, ".h5ad"))
}

.DropLowCountClasses <- function(seuratObj, targetColumn, minCells) {
  dat <- table(seuratObj@meta.data[[targetColumn]])
  dat <- dat[dat > 0]
  toDrop <- names(dat)[dat < minCells]
  if (length(toDrop) > 0) {
    print(paste0('The following classes will be dropped:'))
    for (c in toDrop) {
      print(paste0(c, ': ', dat[[c]], ' cells'))
    }

    toKeep <- rownames(seuratObj@meta.data)[!seuratObj@meta.data[[targetColumn]] %in% toDrop]
    seuratObj <- subset(seuratObj, cells = toKeep)
    print(paste0('Cells remaining: ', ncol(seuratObj)))
  }

  return(seuratObj)
}

.HasNormalizationBeenPerformed <- function(seuratObj, assayName = 'RNA'){
    if (!assayName %in% names(seuratObj@assays)) {
        stop(paste0('Assay not found: ', assayName))
    }

    return(!identical(seuratObj@assays[[assayName]]@counts, seuratObj@assays[[assayName]]@data))
}

.FilterLowCalls <- function(seuratObj, label, minFraction, returnAsFactor = TRUE, labelToAssign = 'Unknown') {
  if (is.null(minFraction)){
    return(seuratObj)
  }

  print(paste0('Filtering ', label, ' below: ', minFraction))
  d <- data.frame(table(Label = unlist(seuratObj[[label]])))
  if (nrow(d) == 0) {
    print(paste0('All values for ', label, ', were NA, skipping filter'))
    return(seuratObj)
  }

  names(d) <- c('Label', 'Count')
  d$Fraction <- d$Count / sum(d$Count)

  d <- d %>% dplyr::arrange(dplyr::desc(Fraction))
  print(d)
  toRemove <- d$Label[d$Fraction < minFraction]
  if (length(toRemove) > 0) {
    print(paste0('Will remove: ', paste0(toRemove, collapse = ', ')))
  }

  if (length(toRemove) > 0) {
    l <- as.character(unlist(seuratObj[[label]]))
    names(l) <- colnames(seuratObj)
    l[l %in% toRemove] <- labelToAssign
    if (returnAsFactor) {
      l <- as.factor(l)
    }

    seuratObj[[label]] <- l
  }

  print('After filter:')
  d <- data.frame(table(Label = unlist(seuratObj[[label]])))
  names(d) <- c('Label', 'Count')
  print(d)

  return(seuratObj)
}

.InferBpParam <- function(nThreads, defaultValue = NULL) {
  if (!is.null(nThreads) && nThreads > 1) {
    if (.Platform$OS.type == 'windows') {
      return(BiocParallel::SnowParam(nThreads))
    } else {
      return(BiocParallel::MulticoreParam(nThreads))
    }
  } else {
    return(defaultValue)
  }
}

.SplitCellsIntoBatches <- function(seuratObj, nBatches, seed = GetSeed()) {
  if (nBatches == 1) {
    stop('It does not make sense to call this function with a single batch')
  }

  cellsPerBatch <- floor(ncol(seuratObj) / nBatches)
  remainder <- ncol(seuratObj) - (cellsPerBatch * nBatches)

  ret <- list()
  set.seed(seed)
  allCells <- colnames(seuratObj)
  ret[[1]] <- sample(allCells, (cellsPerBatch + remainder), replace=FALSE)
  allCells <- setdiff(allCells, ret[[1]])
  for (i in 2:nBatches) {
    ret[[i]] <- sample(allCells, cellsPerBatch, replace=FALSE)
    allCells <- setdiff(allCells, ret[[i]])
  }

  return(ret)
}