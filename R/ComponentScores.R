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
#' @param assayName The assay to use
#' @param layer The layer to use
#' @param referenceForScaling An optional dataframe to use for scaling. This must contain the columns row.means and row.sd, and the rownames must match any required genes
#'
#' @export
ScoreUsingSavedComponent <- function(seuratObj, componentOrName, fieldName, assayName = 'RNA', layer = 'data', referenceForScaling = NULL) {
  if (is.data.frame(componentOrName)) {
    savedComponent <- componentOrName
  } else {
    savedComponent <- .GetSavedComponent(componentOrName)
  }

  geneWeights <- savedComponent$weight
  names(geneWeights) <- savedComponent$feature

  ad <- Seurat::GetAssayData(seuratObj, assay = assayName, layer = layer)
  if (!is.null(referenceForScaling)) {
    if (!all(names(geneWeights) %in% rownames(referenceForScaling))) {
      missingGenes <- names(geneWeights)[!names(geneWeights) %in% rownames(referenceForScaling)]
      stop(paste0('The following genes are missing from referenceForScaling: ', paste0(missingGenes, collapse = ',')))
    }

    ref_means <- referenceForScaling[names(geneWeights), 'row.means']
    ref_sds <- referenceForScaling[names(geneWeights), 'row.sd']

    # 4. Manually scale the QUERY data matrix using the reference metrics
    # Formula: (Query - Ref_Mean) / Ref_SD
    print('Scaling data using custom values')
    scaled_query_matrix <- t(apply(as.matrix(ad[names(geneWeights), ]), 2, function(x) (x - ref_means) / ref_sds))
    scaled_query_matrix <- t(scaled_query_matrix) # Keep genes as rows, cells as columns
    scaled_query_matrix[is.nan(scaled_query_matrix) | is.infinite(scaled_query_matrix)] <- 0
    ad <- Seurat::as.sparse(scaled_query_matrix)
  }

  # Perform the check first for scale.data, since we can add missing genes:
  missingGenes <- names(geneWeights)[!names(geneWeights) %in% rownames(ad)]
  if (length(missingGenes) > 0 && layer == 'scale.data') {
    dataLayer <- Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'data')
    toScale <- missingGenes[missingGenes %in% rownames(dataLayer)]
    if (length(toScale) > 0) {
      print(paste0('The following ', length(toScale), ' genes are missing from scale.data and will be added now: ', paste0(toScale, collapse = ',')))
      scaled <- Seurat::ScaleData(
        dataLayer,
        features = toScale,
        do.scale = TRUE,
        do.center = TRUE
      )

      if (nrow(scaled) != length(toScale)) {
        stop('The number of rows after scaling did not match the expected number')
      }

      ad <- rbind(ad, scaled)
      seuratObj <- suppressWarnings(Seurat::SetAssayData(seuratObj, assay = assayName, layer = layer, new.data = ad))
    }
  }

  toDrop <- names(geneWeights)[!names(geneWeights) %in% rownames(ad)]
  if (length(toDrop) > 0) {
    print(paste0('The following ', length(toDrop), ' genes were in the component but not the assay, skipping: ', paste0(toDrop, collapse = ',')))
    geneWeights <- geneWeights[!names(geneWeights) %in% toDrop]
    print(paste0('Remaining: ', length(geneWeights)))
  }

  if (length(geneWeights) == 0) {
    stop('No genes overlap between the score and this object')
  }

  cellScores <- Matrix::t(Matrix::as.matrix(ad[names(geneWeights), ]))  %*%  geneWeights
  rownames(cellScores) <- colnames(ad)
  seuratObj <- Seurat::AddMetaData(seuratObj, cellScores, col.name = fieldName)

  if (length(names(seuratObj@reductions)) > 0) {
    suppressMessages(print(FeaturePlot(seuratObj, features = fieldName, order = T) & ggplot2::scale_colour_gradientn(colours = c("navy", "dodgerblue", "gold", "red"))))
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
