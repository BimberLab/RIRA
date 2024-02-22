#' @include Utils.R

.GetSavedComponent <- function(componentName) {
  savedFile <- system.file(paste0("components/", componentName, ".tsv"), package = "RIRA")
  if (file.exists(savedFile)) {
    return(utils::read.table(savedFile, header = T, sep = '\t'))
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

  ad <- Seurat::GetAssayData(seuratObj)
  toDrop <- names(geneWeights)[!names(geneWeights) %in% rownames(ad)]
  if (length(toDrop) > 0) {
    print(paste0('The following ', length(toDrop), ' genes were in the component but not the assay, skipping: ', paste0(toDrop, collapse = ',')))
    geneWeights <- geneWeights[!names(geneWeights) %in% toDrop]
    print(paste0('Remaining: ', length(geneWeights)))
  }

  cellScores <- Matrix::t(Matrix::as.matrix(ad[names(geneWeights), ]))  %*%  geneWeights
  seuratObj <- Seurat::AddMetaData(seuratObj, cellScores, col.name = fieldName)

  if (length(names(seuratObj@reductions)) > 0) {
    suppressMessages(print(FeaturePlot(seuratObj, features = fieldName, order = T) & ggplot2::scale_colour_gradientn(colours = c("navy", "dodgerblue", "gold", "red"))))
  } else {
    print('No reductions present, cannot plot')
  }

  graphics::hist(seuratObj@meta.data[[fieldName]], breaks = 300, main = fieldName)

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
