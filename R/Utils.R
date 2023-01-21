#' @import Seurat ggplot2 rlang
#' @importFrom lsa cosine

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

#' @title Remove UCell Doublets
#'
#' @description Attempts to remove dump-gate "doublets" from the main three cell type subset objects (Tcell, Bcell, MoMacDC)
#' @param parentUCell The subset's target cell type's UCell score. (Computed by scGate early in processing).
#' @param otherUCellScores A vector of cell type UCell scores to attempt to discard from downstream processing. By default, these are cell types that comprise the dump gate in RIRA's scGate model. 
#' @param clusterResolution The metadata column name that defines cluster membership (e.g. ClusterNames_0.2)
#' @param cosineSimilarityCutoff A number between -1 and 1 that is used to threshold/gate the prospective doublets as defined by the computed score.
#' @return Returns a seurat object with the UCell doublets removed. 
#' @export
#' 
.RemoveUCellDoublets <- function(seuratObj, parentUCell = NULL, otherUCellScores = c("Platelet_MK_UCell", "Erythrocyte_UCell", "Neutrophil_UCell", "Epithelial_UCell", "Stromal_UCell"), clusterResolution = "ClusterNames_1.2", cosineSimilarityCutoff = 0.5){
  #test for valid parentUCell
  if(is.null(parentUCell)){
    stop("Please define a parent population's UCell score column within the seurat object's metadata. For RIRA, this should be one of: 'Tcell_UCell', 'Bcell_UCell', or 'MoMacDC_UCell'.")
  }
  #check for cosineSimilarityCutoff
  if(!is.numeric(cosineSimilarityCutoff)){
    stop("Please ensure cosine similarity cutoff is a double between -1 and 1. ")
  } else if(abs(cosineSimilarityCutoff) >=1){
    stop("Please ensure cosine similarity cutoff is a double between -1 and 1. ")
  }
  #check for valid clusterResolution
  if(!clusterResolution %in% colnames(seuratObj@meta.data)){
    stop("Please specify a valid clusterResolution. clusterResolution should be a column of the seurat object's metadata specifying cluster membership (e.g. 'ClusterNames_0.2').")
  }
  #check for valid otherUCellScore columns
  if(!any(otherUCellScores %in% colnames(seuratObj@meta.data))){
    stop("No valid UCell scores were supplied to otherUCellScores. Please ensure that otherUCellScores is a vector of column names in the seurat object's metadata.")
  }
  #sanitize otherUCellScores for missing scores
  for(celltypeUCell in otherUCellScores){
    if(!celltypeUCell %in% colnames(seuratObj@meta.data)){
      warning(paste0("'",celltypeUCell, "'", " missing from seurat metadata columns. Omitting from doublet removal."))
      otherUCellScores <- otherUCellScores[otherUCellScores != celltypeUCell]
    }
  }
  
  #iterate through the celltypes in otherUCellScores
  for(celltypeUCell in otherUCellScores){
    #ensure cosine_similiarity column is blank
    seuratObj@meta.data[["cosine_similarity"]] <- NULL
    #locally copy metadata and compute cosine similiarity
    meta <- seuratObj@meta.data %>% 
      group_by(!!rlang::sym(clusterResolution)) %>%
      mutate(cosine_similarity = lsa::cosine(!!rlang::sym(parentUCell), !!rlang::sym(celltypeUCell))) 
    #store store cosine similarity as a uniquely named DoubletScore
    doubletScoreName <- paste0(celltypeUCell,"_DoubletScore")
    seuratObj@meta.data[[doubletScoreName]] <- meta[["cosine_similarity"]]
    
    #in the event of cosine(vector_of_zeroes, vector_of_zeroes), the metric is undefined
    #we are using this to threshold downstream, so I think setting this to zero to ensure it can be plotted and doesn't get filtered is fine.
    seuratObj@meta.data[is.nan(seuratObj@meta.data[,doubletScoreName]), doubletScoreName] <- 0
    #calculate centroids for plotting
    centroids <- aggregate(x = cbind(seuratObj@meta.data[[parentUCell]],seuratObj@meta.data[[doubletScoreName]]) ~ seuratObj@meta.data[[clusterResolution]], data=seuratObj@meta.data, FUN=mean)
    #force consistent names
    colnames(centroids) <- c(clusterResolution, parentUCell, doubletScoreName)
    
    #plot cosine similarity per cluster against parent population
    P1 <- ggplot(seuratObj@meta.data, aes(x = !!rlang::sym(parentUCell), y = !!rlang::sym(doubletScoreName), color = !!rlang::sym(clusterResolution))) + 
      geom_point() +
      geom_text(data = centroids, mapping = aes(x = !!rlang::sym(parentUCell), y = !!rlang::sym(doubletScoreName), label = !!rlang::sym(clusterResolution)), color = "black") +
      geom_hline(yintercept = cosineSimilarityCutoff, linetype='dashed', color = 'red') +
      theme_light()
    print(P1)
    #drop cells above the cosineSimilarityCutoff for the current celltypeUCell
    seuratObj <- subset(seuratObj, subset = !!rlang::sym(doubletScoreName) <= cosineSimilarityCutoff)
  }
  return(seuratObj)
}