#' @include Utils.R

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

.GenerateIntegratedReference <- function(seuratObj, labelCol, batchIdCol = 'DatasetId', assay = 'RNA', dataset.name = 'atlas:', ngenes = 3000, saveFile = NULL) {
  countMatrix <- Seurat::GetAssayData(seuratObj, assay = assay, slot = 'counts')

  print('Subsetting input')
  data.list <- cFIT::split_dataset_by_batch(X = t(as.matrix(countMatrix)),
                                            batch = seuratObj@meta.data[[batchIdCol]],
                                            labels = seuratObj@meta.data[[labelCol]],
                                            metadata = seuratObj@meta.data,
                                            dataset.name = dataset.name)

  print('Selecting genes')
  X.list <- data.list$X.list
  obj.list <- lapply(1:length(X.list), function(j) {
    obj <- Seurat::CreateSeuratObject(counts = Matrix::t(X.list[[j]]))
    obj <- Seurat::NormalizeData(obj, verbose=F)
    obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = ngenes, verbose=F)
    return(obj)
  })
  genes <- Seurat::SelectIntegrationFeatures(object.list = obj.list, nfeatures = ngenes, verbose = TRUE)

  print("Preparing inputs")
  preprocess_for_integration <- function(X.list, genes, scale = T, center = F, verbose = F) {
    datasets <- names(X.list)
    for (i in 1:length(X.list)) {
      genes <- intersect(genes, colnames(X.list[[i]]))
    }
    X.list <- lapply(1:length(X.list), function(j) {
      x <- Matrix::t(X.list[[j]][, genes])
      x <- Seurat::NormalizeData(x)
      x <- Matrix::t(Seurat::ScaleData(x, do.center = center, do.scale = scale, verbose = verbose))

      return(x)
    })

    names(X.list) <- datasets
    return(X.list)
  }
  exprs.list <- preprocess_for_integration(data.list$X.list, genes, scale.factor=10^4, scale=T, center=F)

  print('Running cFIT Integration')
  int.out <- cFIT::CFITIntegrate(X.list=exprs.list, r = 15, verbose = F, max.niter = 100, seed = 0)
  if (!is.null(saveFile)) {
    saveRDS(int.out, file = saveFile)
  }

  return(int.out)
}

#' @title Perform CFITTransfer
#' @description Perform CFITTransfer
#' @param seuratObj The input seurat object
#' @param integratedReference The integrated reference from cFIT
#' @importFrom dplyr %>%
#' @export
Perform_CFITTransfer <- function(seuratObj, integratedReference = GetCFitReference(), verbose = TRUE, seed = 0, assay = 'RNA') {
  #prepare target data for transfer
  target.X <- Matrix::t(Seurat::GetAssayData(seuratObj, slot = 'counts', assay = assay))
  genes_target <- intersect(rownames(integratedReference$W), colnames(target.X))
  if (length(genes_target) != length(rownames(integratedReference$W))) {
    stop('The genes on the input seurat object and integrated reference do not match')
  }

  x <- target.X[, genes_target]
  x <- Seurat::NormalizeData(x)
  x <- Matrix::t(Seurat::ScaleData(Matrix::t(x), do.center = F, do.scale = T, verbose = F))
  target.exprs.list <- x
  names(target.exprs.list) <- names(target.X)
  W_ref <- integratedReference$W[rownames(integratedReference$W) %in% genes_target,]
  idx <- match(rownames(W_ref), colnames(target.exprs.list))
  ordered_targ.exp.list <- target.exprs.list[,idx]

  #transfer learned gene expression signatures from integration onto target dataset
  tf.out <- cFIT::CFITTransfer(Xtarget=ordered_targ.exp.list, Wref=W_ref, max.niter=100, seed=0, verbose=F, future.plan='sequential')

  return(tf.out)
}

