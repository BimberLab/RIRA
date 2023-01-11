#' @include Utils.R
#' @include Phenotyping.R
#'
#' @title CalculateUCellScores
#'
#' @description This will calculate a handful of standardized UCell scores for a seurat object
#' @param seuratObj The seurat object
#' @param forceRecalculate If true, UCell will always be re-run, even if the field is already present.
#' @param seed If non-null, set.seed() will be called prior to running UCell
#' @param ncores The number of core to use with UCell::AddModuleScore_UCell
#' @param assayName The assay from which to calculate UCell scores.
#' @param storeRanks Passed directly to UCell::AddModuleScore_UCell. Increases object size but makes future calculations quicker.
#' @param plotCor If true, a plot of correlations between the UCell score and each component gene will be shown
#' @export
CalculateUCellScores <- function(seuratObj, forceRecalculate = FALSE, seed = GetSeed(), ncores = 1, assayName = 'RNA', storeRanks = FALSE, plotCor = TRUE) {
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

    BPPARAM <- .InferBpParam(ncores, defaultValue = NULL)
    seuratObj <- UCell::AddModuleScore_UCell(seuratObj, features = toCalculate, BPPARAM = BPPARAM, assay = assayName, storeRanks = storeRanks)
  } else {
    print('UCell score already present, will not recalculate')
  }

  if (plotCor) {
    for (moduleName in names(toCalculate)) {
      geneList <- toCalculate[[moduleName]]
      geneData <- as.data.frame(t(as.matrix(seuratObj@assays[[assayName]]@data[geneList,])))
      geneData$UCell <- unlist(seuratObj[[paste0(moduleName, '_UCell')]])

      ret <- stats::cor(geneData, method = "spearman")
      p.mat <- ggcorrplot::cor_pmat(ret, method = 'spearman')
      print(ggcorrplot::ggcorrplot(ret, hc.order = TRUE, type = "lower", p.mat = p.mat) + ggtitle(paste0(moduleName, " UCell-Gene Correlations")))
    }
  }

  if (any(is.na(seuratObj[['TandNK_Activation_UCell']]))) {
    print('Data has NAs, cannot make feature plot: TandNK_Activation_UCell')
  } else {
    print(Seurat::FeaturePlot(seuratObj, features = 'TandNK_Activation_UCell', min.cutoff = 'q02', max.cutoff = 'q98') + ggtitle('T Cell Activation Score'))
  }

  if (any(is.na(seuratObj[['Cytotoxicity_UCell']]))) {
    print('Data has NAs, cannot make feature plot: Cytotoxicity_UCell')
  } else {
    print(Seurat::FeaturePlot(seuratObj, features = 'Cytotoxicity_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Cytotoxicity Score'))
  }

  if (any(is.na(seuratObj[['EffectorT_UCell']]))) {
    print('Data has NAs, cannot make feature plot: EffectorT_UCell')
  } else {
    print(Seurat::FeaturePlot(seuratObj, features = 'EffectorT_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Effector T Score'))
  }

  if (any(is.na(seuratObj[['NaiveT_UCell']]))) {
    print('Data has NAs, cannot make feature plot: NaiveT_UCell')
  } else {
    print(Seurat::FeaturePlot(seuratObj, features = 'NaiveT_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Naive T Score'))
  }

  if (any(is.na(seuratObj[['CentralMemT_UCell']]))) {
    print('Data has NAs, cannot make feature plot: CentralMemT_UCell')
  } else {
    print(Seurat::FeaturePlot(seuratObj, features = 'CentralMemT_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Central Mem T Score'))
  }

  if (any(is.na(seuratObj[['Glycolysis_UCell']]))) {
    print('Data has NAs, cannot make feature plot: Glycolysis_UCell')
  } else {
    print(Seurat::FeaturePlot(seuratObj, features = 'Glycolysis_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Glycolysis'))
  }

  if (any(is.na(seuratObj[['Interferon_Response_UCell']]))) {
    print('Data has NAs, cannot make feature plot: Interferon_Response_UCell')
  } else {
    print(Seurat::FeaturePlot(seuratObj, features = 'Interferon_Response_UCell', min.cutoff = 'q05', max.cutoff = 'q95') + ggtitle('Interferon_Response'))
  }

  return(seuratObj)
}