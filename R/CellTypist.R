#' @include Utils.R


#' @title Run Celltypist
#'
#' @description Runs celltypist on a seurat object and stores the calls as metadata
#' @param seuratObj The seurat object
#'
#' @export
RunCellTypist <- function(seuratObj) {
  if (!reticulate::py_available(initialize = TRUE)) {
    stop(paste0('Python/reticulate not configured. Run "reticulate::py_config()" to initialize python'))
  }

  if (!reticulate::py_module_available('celltypist')) {
    stop('The celltypist python package has not been installed!')
  }

  outFile <- tempfile()
  outDir <- dirname(outFile)
  seuratObj <- Seurat::DietSeurat(seuratObj)
  seuratAnnData <- SeuratToAnnData(seuratObj, paste0(outDir, 'seurat.annData'))
  seuratAnnData <- SeuratToAnnData(seuratObj, seuratAnnData)

  #exe <- reticulate::py_exe()
  exe <- "celltypist"
  #"-m", "celltypist.command_line",
  pyOut <- system2(exe, c("--update-models", "--quiet"), stdout = TRUE, stderr = TRUE)
  print(pyOut)

  # "-m", "celltypist.command_line",
  # "--plot-results"
  args <- c("-i", seuratAnnData, "-m", "Immune_All_Low.pkl", "--outdir", outDir, "--majority-voting", "--prefix", "celltypist.", "--quiet")
  print(args)
  pyOut <- system2(exe, args, stdout = TRUE, stderr = TRUE)
  print(pyOut)

  labels <- paste0(outDir, '/celltypist.predicted_labels.csv')
  if (!file.exists(labels)) {
    stop(paste0('Missing file: ', labels))
  }

  labels <- read.csv(labels, header = T, row.names = 1)
  seuratObj <- Seurat::AddMetaData(seuratObj, labels)

  unlink(seuratAnnData)

  return(seuratObj)
}


