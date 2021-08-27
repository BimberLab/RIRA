#' @include Utils.R

#' @title Get the integrated cFIT reference object for an atlas version
#' @description Return the integrated reference for the provided RIRA version
#' @param version The version string
#' @export
GetCFitReference <- function(version = NULL) {
  fn <- .GetExpectedCfitRef(version)
  if (!file.exists(fn)) {
    stop(paste0('Unable to find file: ', fn))
  }

  return(readRDS(fn))
}

.GetExpectedCfitRef <- function(version) {
  parentFolder <- .GetAtlasBaseDir(version)
  fn <- paste0(parentFolder, '/cFIT/ref.rds')

  return(fn)
}


#' @title GenerateIntegratedReference
#' @description This function uses cFIT to create an integrated reference for the provided seurat object
#' @param seuratObj The seurat object
#' @param batchIdCol If provided, data will be subset using this metadata field
#' @param targetBatchSize This is an alternative to batchIdCol. If batchIdCol is not provided, data will be randomly subset
#' @param assay The seurat assay holding count data
#' @param ngenes The number of variable genes to select
#' @param saveFile If provided, the integrated reference object will be saved to this file as RDS
#' @param seed The random seed
#' @param verbose Determined output verbosity
#' @param minBatchSize If any resulting batch has fewer than this many cells, an error is thrown
#' @param selectVariableFeaturesFromWholeInput If true, the variable genes will be calculated from the intact input object. If false, FindVariableFeatures() will be calculated independently on each subset seurat object, followed by Seurat::SelectIntegrationFeatures
#' @export
GenerateIntegratedReference <- function(seuratObj, batchIdCol = NULL, targetBatchSize = 5000, assay = 'RNA', ngenes = 3000, saveFile = NULL, seed = 0, verbose = FALSE, minBatchSize = 100, selectVariableFeaturesFromWholeInput = TRUE) {
  # TODO: should we be able to determine if this already has been done?
  seuratObj <- Seurat::NormalizeData(seuratObj, verbose = verbose)

  print('Subsetting input')
  if (!is.null(batchIdCol)) {
    print(paste0('Batching by: ', batchIdCol))
    seuratObjects <- Seurat::SplitObject(seuratObj, split.by = batchIdCol)
  } else if (!is.null(targetBatchSize)){
    nBatches <- floor(ncol(seuratObj) / targetBatchSize)  #this will tend to err on slightly higher number of cells/batch
    print(paste0('will divide data into ', nBatches))

    seuratObj$BatchCol <- sample(1:nBatches, ncol(seuratObj), replace = TRUE)
    seuratObjects <- Seurat::SplitObject(seuratObj, split.by = 'BatchCol')
  } else {
    stop('Must provide either batchIdCol or targetBatchSize')
  }

  invisible(lapply(seuratObjects, function(so) {
    if (ncol(so) < minBatchSize) {
      stop(paste0('One or more batches was below minBatchSize: ', minBatchSize))
    }
  }))

  print('Selecting genes')
  if (selectVariableFeaturesFromWholeInput) {
    seuratObj <- Seurat::FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = ngenes, verbose = verbose)
    genes <- Seurat::VariableFeatures(seuratObj)
  } else {
    seuratObjects <- lapply(seuratObjects, function(obj){
      obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = ngenes, verbose = verbose)
      return(obj)
    })

    genes <- Seurat::SelectIntegrationFeatures(object.list = seuratObjects, nfeatures = ngenes, verbose = TRUE)
  }

  print("Preparing inputs")
  exprs.list <- lapply(seuratObjects, function(so){
    x <- Seurat::GetAssayData(so, assay = assay, slot = 'data')
    x <- x[genes,]
    x <- Seurat::ScaleData(x, do.center = F, do.scale = T, verbose = verbose)
    return(Matrix::t(x))
  })

  print('Running cFIT Integration')
  int.out <- cFIT::CFITIntegrate(X.list=exprs.list, r = 15, verbose = verbose, max.niter = 100, seed = seed)
  if (!is.null(saveFile)) {
    saveRDS(int.out, file = saveFile)
  }

  return(int.out)
}

#' @title Perform CFITTransfer
#' @description Perform CFITTransfer
#' @param seuratObj The input seurat object
#' @param integratedReference The integrated reference from cFIT
#' @param verbose Controls the verbosity of the output
#' @param seed The random reed
#' @param assay The name of the assay holding count data
#' @importFrom dplyr %>%
#' @export
Perform_CFITTransfer <- function(seuratObj, integratedReference = GetCFitReference(), verbose = TRUE, seed = 0, assay = 'RNA') {
  #prepare target data for transfer
  target.X <- Seurat::GetAssayData(seuratObj, slot = 'counts', assay = assay)
  genes_target <- intersect(rownames(integratedReference$W), rownames(target.X))
  if (length(genes_target) != length(rownames(integratedReference$W))) {
    stop('The genes on the input seurat object and integrated reference do not match')
  }

  x <- target.X[genes_target,]
  x <- Seurat::NormalizeData(x)
  x <- Matrix::t(Seurat::ScaleData(x, do.center = F, do.scale = T, verbose = F))
  target.exprs.list <- x
  names(target.exprs.list) <- names(target.X)
  W_ref <- integratedReference$W[rownames(integratedReference$W) %in% genes_target,]
  idx <- match(rownames(W_ref), colnames(target.exprs.list))
  ordered_targ.exp.list <- target.exprs.list[,idx]

  #transfer learned gene expression signatures from integration onto target dataset
  tf.out <- cFIT::CFITTransfer(Xtarget = ordered_targ.exp.list, Wref = W_ref, max.niter = 100, seed = seed, verbose = verbose, future.plan='sequential')

  return(tf.out)
}

