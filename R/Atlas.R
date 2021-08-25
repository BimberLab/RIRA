#' @include Utils.R



GetRiraCountMatrix <- function(version = NULL) {
  if (is.null(version)) {
    version <- .GetLatestVersion()
  }

  countMat <- system.file(paste0('inst/data', 'RIRA-', version, '.counts.csv'), package = 'RIRA')
  if (!file.exists(countMat)) {
    print(paste0('Unknown version: ', version))
  }
  countMat <- as.matrix(read.csv(countMat, row.names = 1))

  metadata <- system.file(paste0('inst/data', 'RIRA-', version, '.metadata.csv'), package = 'RIRA')
  if (!file.exists(metadata)) {
    print(paste0('Unknown version: ', version))
  }
  metadata <- read.csv(metadata, header = T, row.names = 1)

  return(Seurat::CreateSeuratObject(counts = countMat, project = paste0('RIRA-', version), meta.data = metadata))
}

