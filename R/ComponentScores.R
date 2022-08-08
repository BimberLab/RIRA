#' @include Utils.R

.GetSavedComponent <- function(componentName) {
  savedFile <- system.file(paste0("components/", componentName, ".tsv"), package = "RIRA")
  if (file.exists(savedFile)) {
    return(read.table(savedFile, header = T, sep = '\t'))
  }

  stop(paste0('Unknown component: ', componentName))
}

#' @title ScoreUsingSavedComponent
#'
#' @description Scores a seurat object using the gene weights from a saved component
#' @param seuratObj The seurat object
#' @param componentOrName Either a data frame with the columns 'feature' and 'weight', or the name of a previously saved component available in RIRA
#' @param fieldName The name of the field in which to save the resulting scores
#'
#' @export
ScoreUsingSavedComponent <- function(seuratObj, componentOrName, fieldName) {
  if (is.data.frame(componentOrName)) {
    savedComponent <- componentOrName
  } else {
    savedComponent <- .GetSavedComponent(componentOrName)
  }

  geneWeights <- savedComponent$weight
  names(geneWeights) <- savedComponent$feature

  cellScores <- Matrix::t(Matrix::as.matrix(Seurat::GetAssayData(seuratObj)[names(geneWeights), ]))  %*%  geneWeights
  seuratObj <- Seurat::AddMetaData(seuratObj, cellScores, col.name = fieldName)

  print(FeaturePlot(seuratObj, features=fieldName, order = T, cols = c("navy", "dodgerblue", "gold", "red")))

  return(seuratObj)
}

#' @title ExtractGeneWeights
#'
#' @description Extracts the loadings from a given component in a seurat object and returns a data.frame suitable for ScoreUsingSavedComponent
#' @param seuratObj The seurat object
#' @param componentNum The component number to use
#' @param topNGenes The number of top genes to use (ranked by loading)
#' @param bottomNGenes The number of bottom genes to use (ranked by loading)
#' @param reduction The reduction to use
#'
#' @export
ExtractGeneWeights <- function(seuratObj, componentNum, topNGenes = 100, bottomNGenes = 100, reduction = "pca") {
  if (missing(componentNum)) {
    stop('Must provide componentNum argument')
  }

  x <- Seurat::Loadings(seuratObj, reduction = reduction)
  if (nrow(x) < topNGenes) {
    topNGenes <- nrow(x)
    print(paste0('Too many genes requested, reducing to: ', topNGenes))
  }
  topGenes <- sort(x[, componentNum], decreasing = T)[1:topNGenes]

  if (nrow(x) < bottomNGenes) {
    bottomNGenes <- nrow(x)
    print(paste0('Too many genes requested, reducing to: ', bottomNGenes))
  }
  bottomGenes <- sort(x[, componentNum], decreasing = F)[1:bottomNGenes]
  bottomGenes <- bottomGenes[!(names(bottomGenes) %in% names(topGenes))]

  allGenes <- c(topGenes, bottomGenes)
  allGenes <- data.frame(feature = names(allGenes), weight = unname(allGenes))

  return(allGenes)
}