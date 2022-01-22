#' @include Utils.R

# See: https://github.com/carmonalab/scGate/pull/1
#' @importFrom UCell AddModuleScore_UCell
#'
#' @title Run scGate
#'
#' @description Helper function to run scGate
#' @param seuratObj The seurat object
#' @param model Either an scGate model, or a character that will be passed to GetScGateModel()
#' @param min.cells Passed directly to scGate::scGate. Stop iterating if fewer than this number of cells is left
#' @param assay Passed directly to scGate::scGate. Seurat assay to use
#' @param pos.thr Passed directly to scGate::scGate. Minimum UCell score value for positive signatures
#' @param neg.thr Passed directly to scGate::scGate. Maximum UCell score value for negative signatures
#' @param ncores Passed directly to scGate::scGate. Number of processors for parallel processing (requires future.apply)
#' @param output.col.name Passed directly to scGate::scGate. Column name with 'pure/impure' annotation
#' @param genes.blacklist Passed directly to scGate::scGate. Genes blacklisted from variable features. The default loads the list of genes in scGate::genes.blacklist.default; you may deactivate blacklisting by setting genes.blacklist=NULL
#'
#' @export
RunScGate <- function(seuratObj, model, min.cells = 10, assay = 'RNA', pos.thr = 0.13, neg.thr = 0.13, ncores = 1, output.col.name = "is.pure", genes.blacklist = 'default') {
  if (is.character(model)) {
    model <- GetScGateModel(model)
    if (is.null(model)) {
      stop(paste0('Unknown gate model: ', model))
    }
  }

  seuratObj <- scGate::scGate(data = seuratObj,
                        model = model,
                        min.cells = min.cells,
                        assay = assay,
                        pos.thr = pos.thr,
                        neg.thr = neg.thr,
                        seed = GetSeed(),
                        ncores = ncores,
                        output.col.name = output.col.name,
                        genes.blacklist = genes.blacklist
  )

  if (length(names(seuratObj@reductions)) > 0) {
    print(Seurat::DimPlot(seuratObj, group.by = output.col.name))
    colNames <- names(seuratObj@meta.data)[grepl(names(seuratObj@meta.data), pattern = paste0('^', output.col.name, '.'))]
    for (col in colNames) {
      print(Seurat::DimPlot(seuratObj, group.by = col))
    }

    colNames <- names(seuratObj@meta.data)[grepl(names(seuratObj@meta.data), pattern = paste0('UCell$'))]
    for (col in colNames) {
      suppressWarnings(print(Seurat::FeaturePlot(seuratObj, features = col, min.cutoff = 'q05', max.cutoff = 'q95')))
    }
  } else {
    print('There are no reductions in this seurat object, cannot create dimplots')
  }

  return(seuratObj)
}

#' @title GetAvailableScGates
#'
#' @description Return a list of available scGate models
#' @export
GetAvailableScGates <- function() {
  dir <- system.file("gates", package = "RIRA")
  files <- list.files(dir, recursive = FALSE, full.names = FALSE)
  files <- files[files != 'master_table.tsv']
  files <- sapply(files, function(x){
    return(gsub(x, pattern = '.tsv', replacement = ''))
  })

  return(files)
}

#' @title GetScGateModel
#'
#' @description Returns the selected scGate model
#' @param modelName The name of the gate to return. See GetAvailableScGates() for a list of known gates
#' @export
GetScGateModel <- function(modelName, allowSCGateDB = TRUE) {
  gateFile <- system.file(paste0("gates/", modelName, ".tsv"), package = "RIRA")
  if (file.exists(gateFile)) {
    masterFile <- system.file("gates/master_table.tsv", package = "RIRA")
    return(scGate::load_scGate_model(gateFile, master.table = masterFile))
  }

  if (!allowSCGateDB) {
    stop(paste0('Unable to find gate: ', modelName))
  }

  models.DB <- scGate::get_scGateDB()
  if (!modelName %in% names(models.DB$human$generic)){
    stop(paste0('Unable to find model: ', modelName))
  }

  print(paste0('Using build-in model: ', modelName))
  return(models.DB$human$generic[[modelName]])
}

#' @title Run scGate With DefaultModels
#'
#' @description Helper function to run scGate, running all human models in scGate::get_sc()
#' @param seuratObj The seurat object
#' @param min.cells Passed directly to scGate::scGate. Stop iterating if fewer than this number of cells is left
#' @param assay Passed directly to scGate::scGate. Seurat assay to use
#' @param pos.thr Passed directly to scGate::scGate. Minimum UCell score value for positive signatures
#' @param neg.thr Passed directly to scGate::scGate. Maximum UCell score value for negative signatures
#' @param ncores Passed directly to scGate::scGate. Number of processors for parallel processing (requires future.apply)
#' @param output.col.name Passed directly to scGate::scGate. Column name with 'pure/impure' annotation
#' @param genes.blacklist Passed directly to scGate::scGate. Genes blacklisted from variable features. The default loads the list of genes in scGate::genes.blacklist.default; you may deactivate blacklisting by setting genes.blacklist=NULL
#'
#' @export
RunScGateWithDefaultModels <- function(seuratObj, min.cells = 10, assay = 'RNA', pos.thr = 0.13, neg.thr = 0.13, ncores = 1, genes.blacklist = 'default') {
  models.DB <- scGate::get_scGateDB()
  for (modelName in names(models.DB$human$generic)){
    print(paste0('Running model: ', modelName))

    seuratObj <- RunScGate(seuratObj = seuratObj,
              model = modelName,
              min.cells = min.cells,
              assay = assay,
              pos.thr = pos.thr,
              neg.thr = neg.thr,
              ncores = ncores,
              output.col.name = paste0(modelName, '.is.pure'),
              genes.blacklist = genes.blacklist
          )
  }

  return(seuratObj)
}