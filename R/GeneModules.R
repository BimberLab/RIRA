#' @include Utils.R
#' @include Phenotyping.R
#'
#' @title CalculateUCellScores
#'
#' @description This will calculate a handful of standardized UCell scores for a seurat object
#' @param seuratObj The seurat object
#' @param forceRecalculate If true, UCell will always be re-run, even if the field is already present.
#' @export
CalculateUCellScores <- function(seuratObj, forceRecalculate = FALSE) {
  if (forceRecalculate || !'Cytotoxicity_UCell' %in% names(seuratObj@meta.data) || !'TandNK_Activation_UCell' %in% names(seuratObj@meta.data)) {
    seuratObj <- UCell::AddModuleScore_UCell(seuratObj, features = list(
      TandNK_Activation = GetGeneSet('TandNK_Activation.1'),
      Cytotoxicity = GetGeneSet('Cytotoxicity')
    )) } else {
    print('UCell score already present, will not recalculate')
  }

  print(Seurat::FeaturePlot(seuratObj, features = 'TandNK_Activation_UCell', min.cutoff = 'q02', max.cutoff = 'q98') + ggtitle('T Cell Activation Score'))
  print(Seurat::FeaturePlot(seuratObj, features = 'Cytotoxicity_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Cytotoxicity Score'))

  return(seuratObj)
}