#' @include Utils.R
#' @include Phenotyping.R
#'
#' @title CalculateUCellScores
#'
#' @description This will calculate a handful of standardized UCell scores for a seurat object
#' @param seuratObj The seurat object
#' @param forceRecalculate If true, UCell will always be re-run, even if the field is already present.
#' @param seed If non-null, set.seed() will be called prior to running UCell
#' @export
CalculateUCellScores <- function(seuratObj, forceRecalculate = FALSE, seed = GetSeed()) {
  toCalculate <- list(
    TandNK_Activation = GetGeneSet('TandNK_Activation.1'),
    Cytotoxicity = GetGeneSet('Cytotoxicity'),
    EffectorT = GetGeneSet('EffectorT'),
    CentralMemT = GetGeneSet('CentralMemT'),
    NaiveT = GetGeneSet('NaiveT'),
    Glycolysis = GetGeneSet('Glycolysis'),
    Interferon_Response = GetGeneSet('Interferon_Response')
  )

  if (forceRecalculate || any(!paste0(names(toCalculate), '_UCell') %in% names(seuratObj@meta.data))) {
    if (!is.null(seed)) {
      print(paste0('Setting random seed: ', seed))
      set.seed(seed)
    }

    seuratObj <- UCell::AddModuleScore_UCell(seuratObj, features = toCalculate)
  } else {
    print('UCell score already present, will not recalculate')
  }

  print(Seurat::FeaturePlot(seuratObj, features = 'TandNK_Activation_UCell', min.cutoff = 'q02', max.cutoff = 'q98') + ggtitle('T Cell Activation Score'))
  print(Seurat::FeaturePlot(seuratObj, features = 'Cytotoxicity_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Cytotoxicity Score'))
  print(Seurat::FeaturePlot(seuratObj, features = 'EffectorT_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Effector T Score'))
  print(Seurat::FeaturePlot(seuratObj, features = 'NaiveT_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Naive T Score'))
  print(Seurat::FeaturePlot(seuratObj, features = 'CentralMemT_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Central Mem T Score'))
  print(Seurat::FeaturePlot(seuratObj, features = 'Glycolysis_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Glycolysis'))
  print(Seurat::FeaturePlot(seuratObj, features = 'Interferon_Response_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Interferon_Response'))

  return(seuratObj)
}