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
#' @param dropAllExistingUcells If true, any existing column ending in UCell will be removed
#' @export
CalculateUCellScores <- function(seuratObj, forceRecalculate = FALSE, seed = GetSeed(), ncores = 1, assayName = 'RNA', storeRanks = TRUE, plotCor = TRUE, dropAllExistingUcells = FALSE) {
  toCalculate <- list(
    TandNK_Activation = GetGeneSet('TandNK_Activation.1'),
    TandNK_Activation2 = GetGeneSet('TandNK_Activation.2'),
    TandNK_Activation3 = GetGeneSet('TandNK_Activation.3'),
    TandNK_ActivationCore = GetGeneSet('TandNK_Activation.Core'),
    Cytotoxicity = GetGeneSet('Cytotoxicity'),
    Cytotoxicity.GzmABH = GetGeneSet('Cytotoxicity.GzmABH'),
    Cytotoxicity.GzmKM = GetGeneSet('Cytotoxicity.GzmKM'),
    Glycolysis = GetGeneSet('Glycolysis'),
    Interferon_Response = GetGeneSet('Interferon_Response'),
    Mitochondrial = GetGeneSet('MMul10_Mitochondrial'),
    EffectorCytokines = GetGeneSet('EffectorCytokines'),
    MAIT_Markers = GetGeneSet('MAIT_Markers'),
    Metallothionein = GetGeneSet('Metallothionein'),
    Metallothionein.Core = GetGeneSet('Metallothionein.Core'),
    IEGs = GetGeneSet('IEGs'),
    MHCII = GetGeneSet('MHC-II'),
    TCellMemory = GetGeneSet('TCellMemoryS100'),
    InflMyeloid = GetGeneSet('InflMyeloid')
  )

  if (dropAllExistingUcells) {
    print('Dropping all existing UCell columns:')
    for (cn in grep(names(seuratObj@meta.data), pattern = '_UCell$', value = TRUE)) {
      seuratObj[[cn]] <- NULL
    }
  }

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
    PlotUcellCorrelation(seuratObj, toCalculate)
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

#' @title PlotUcellCorrelation
#'
#' @description This will plot the correlation between a UCell score and a gene set of interest
#' @param seuratObj The seurat object
#' @param toCalculate A named list where each item is a character vector of genes
#' @param assayName The assay to use
#' @return A list with moduleName and the spearman correlation matrix
#' @export
PlotUcellCorrelation <- function(seuratObj, toCalculate, assayName = 'RNA') {
  corData <- list()
  for (moduleName in names(toCalculate)) {
    geneList <- gsub(toCalculate[[moduleName]], pattern = '-$', replacement = '')
    missingGenes <- dplyr::setdiff(geneList, rownames(seuratObj@assays[[assayName]]))
    if (length(missingGenes) > 0) {
      print(paste0('The following genes were not present in the object: ', paste0(missingGenes, collapse = ',')))
    }

    geneList <- intersect(geneList, rownames(seuratObj@assays[[assayName]]))

    # Drop any genes with all zeros
    genesToSkip <- NULL
    for (gene in geneList) {
      if (sum(Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'data')[gene,] > 0) == 0) {
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

    geneData <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seuratObj, assay = assayName, layer = 'data')[geneList,, drop = FALSE])))
    if (! paste0(moduleName, '_UCell') %in% names(seuratObj@meta.data)) {
      stop(paste0('Missing column: ', paste0(moduleName, '_UCell')))
    }

    geneData$UCell <- unlist(seuratObj[[paste0(moduleName, '_UCell')]])

    ret <- stats::cor(geneData, method = "spearman")
    corData[[moduleName]] <- ret

    tryCatch({
      p.mat <- ggcorrplot::cor_pmat(ret, method = 'spearman')
      print(ggcorrplot::ggcorrplot(ret, hc.order = TRUE, type = "lower", p.mat = p.mat) + ggtitle(paste0(moduleName, " UCell-Gene Correlations")))
    }, error = function(e){
      print(paste0('Error generating ggcorrplot for: ', moduleName))
      print(utils::str(geneData))
      print(utils::str(ret))
      print(conditionMessage(e))
    })
  }

  return(corData)
}