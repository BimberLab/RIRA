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
CalculateUCellScores <- function(seuratObj, forceRecalculate = FALSE, seed = GetSeed(), ncores = 1, assayName = 'RNA', storeRanks = TRUE, plotCor = TRUE) {
  toCalculate <- list(
    TandNK_Activation = GetGeneSet('TandNK_Activation.1'),
    TandNK_ActivationCore = GetGeneSet('TandNK_Activation.Core'),
    Cytotoxicity = GetGeneSet('Cytotoxicity'),
    EffectorT = GetGeneSet('EffectorT'),
    CentralMemT = GetGeneSet('CentralMemT'),
    NaiveT = GetGeneSet('NaiveT'),
    Glycolysis = GetGeneSet('Glycolysis'),
    Interferon_Response = GetGeneSet('Interferon_Response'),
    Interferon_Response_IFI6 = GetGeneSet('Interferon_Response_IFI6_correlated'),
    Ribosomal = GetGeneSet('MMul10_Ribosomal'),
    Mitochondrial = GetGeneSet('MMul10_Mitochondrial'),
    EffectorCytokines = GetGeneSet('EffectorCytokines'),
    ExhaustionOrInhibitory = GetGeneSet('ExhaustionOrInhibitory'),
    MAIT_Markers = GetGeneSet('MAIT_Markers')
  )

  needsRecalc <- forceRecalculate || any(!paste0(names(toCalculate), '_UCell') %in% names(seuratObj@meta.data))

  # NOTE: situations like merging two seurat objects could result in these columns existing, but having NAs
  if (!needsRecalc) {
    for (colName in paste0(names(toCalculate), '_UCell')) {
      if (colName %in% names(seuratObj@meta.data) && any(is.na(seuratObj@meta.data[[colName]]))) {
        needsRecalc <- TRUE
        break
      }
    }
  }

  if (needsRecalc) {
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
      geneList <- gsub(toCalculate[[moduleName]], pattern = '-$', replacement = '')
      missingGenes <- dplyr::setdiff(geneList, rownames(seuratObj@assays[[assayName]]))
      if (length(missingGenes) > 0) {
        print(paste0('The following genes were not present in the object: ', paste0(missingGenes, collapse = ',')))
      }

      geneList <- intersect(geneList, rownames(seuratObj@assays[[assayName]]))

      # Drop any genes with all zeros
      genesToSkip <- c()
      for (gene in geneList) {
        if (sum(Seurat::GetAssayData(seuratObj, assay = assayName, slot = 'data')[gene,] > 0) == 0) {
          genesToSkip <- c(genesToSkip, gene)
        }
      }

      if (length(genesToSkip) > 0) {
        print(paste0('Skipping genes with zero counts: ', paste0(genesToSkip, collapse = ',')))
        geneList <- geneList[!geneList %in% genesToSkip]
      }

      if (length(geneList) == 0) {
        warning(paste0('No shared genes, skipping: ', moduleName))
        next
      }

      geneData <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seuratObj, assay = assayName, slot = 'data')[geneList,])))

      geneData$UCell <- unlist(seuratObj[[paste0(moduleName, '_UCell')]])

      ret <- stats::cor(geneData, method = "spearman")
      p.mat <- ggcorrplot::cor_pmat(ret, method = 'spearman')
      print(ggcorrplot::ggcorrplot(ret, hc.order = TRUE, type = "lower", p.mat = p.mat) + ggtitle(paste0(moduleName, " UCell-Gene Correlations")))
    }
  }

  hasReductions <- length(seuratObj@reductions) > 0
  if (!hasReductions) {
    print('No reductions calculated, cannot plot tSNE/UMAP')
  }

  for (geneModule in names(toCalculate)) {
    ucell <- paste0(geneModule, '_UCell')
    if (any(is.na(seuratObj[[ucell]]))) {
      print(paste0('Data has NAs, cannot make feature plot: ', ucell))
    } else {
      if (hasReductions) {
        print(Seurat::FeaturePlot(seuratObj, features = ucell, min.cutoff = 'q02', max.cutoff = 'q98') + ggtitle(geneModule))
      }
    }
  }

  return(seuratObj)
}