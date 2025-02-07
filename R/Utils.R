
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

SeuratToMatrix <- function(seuratObj, outDir, assayName, slot = 'counts'){
  if (endsWith(outDir, "/")){
    outDir <- gsub(outDir, pattern = "/$", replacement = "")
  }

  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, slot = slot), path = outDir, overwrite = TRUE, type = 'sparse')

  return(paste0(outDir, '/matrix.mtx'))
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

.HasNormalizationBeenPerformed <- function(seuratObj, assayName = Seurat::DefaultAssay(seuratObj)){
  if (!assayName %in% names(seuratObj@assays)) {
    stop(paste0('Assay not found: ', assayName))
  }

  assayObj <- Seurat::GetAssay(seuratObj, assay = assayName)
  if (class(assayObj)[1] == 'Assay') {
    return(!identical(Seurat::GetAssayData(seuratObj, assay = assayName, slot = 'counts'), Seurat::GetAssayData(seuratObj, assay = assayName, slot = 'data')))
  } else if (class(assayObj)[1] == 'Assay5') {
    return('data' %in% names(assayObj@layers))
  } else {
    stop(paste0('Unknown assay class: ', class(assayObj)[1]))
  }
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

#' @title Download RIRA data from GEO
#'
#' @description This will download the RIRA raw data from GEO, as a Seurat object
#' @param outfile This location to initially download the raw data. If omitted or NULL, it will download the RDS to the temporary file folder and delete it on completion. However, because the download is large it may be useful to download to a location you plan to keep.
#' @param subset The cell subset to download, either All, T_NK, Myeloid, BCell or Other
#' @export
DownloadFromGEO <- function(outfile = NULL, subset = c('All', 'T_NK', 'Myeloid', 'Bcell')) {
  subset <- match.arg(subset)

  deleteDownloadedFile <- FALSE
  if (all(is.null(outfile))) {
    print('Downloading data to a temporary file, which will be deleted after loading. If you expect to load these data multiple times, we recommend specifying a permanent download location using the outfile argument')
    outfile <- tempfile()
    deleteDownloadedFile <- TRUE
  }

  if (subset == 'All'){
    fn <- 'GSE277821_RIRA.All.seurat.rds.gz'
  } else if (subset == 'T_NK') {
    fn <- 'GSE277821_RIRA.T_NK.seurat.rds.gz'
  } else if (subset == 'Myeloid') {
    fn <- 'GSE277821_RIRA.Myeloid.seurat.rds.gz'
  } else if (subset == 'Bcell') {
    fn <- 'GSE277821_RIRA.Bcell.seurat.rds.gz'
  } else {
    stop(paste0('Unknown value for subset: ', subset))
  }

  url <- paste0('ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277821/suppl/', fn)
  message(paste0('Downloading file: ', url))
  httr::GET(url = url, httr::progress(), httr::write_disk(path = outfile, overwrite = TRUE))

  seuratObj <- readRDS(outfile)
  if (deleteDownloadedFile) {
    unlink(outfile)
  }

  return(seuratObj)
}
