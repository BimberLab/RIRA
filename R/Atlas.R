#' @include Utils.R


#' @title Get the count data for a given RIRA version, as a Seurat object
#'
#' @description Get the count data for a given RIRA version, as a Seurat object
#' @param version The version string. If not provided, defaults to the current version.
#' @export
GetRiraCountMatrix <- function(version = NULL) {
  parentFolder <- .GetAtlasBaseDir(version)

  countMat <- paste0(parentFolder, '/counts')
  countMat <- DropletUtils::read10xCounts(countMat, col.names = T)
  countMat2 <- countMat@assays@data$counts
  rownames(countMat2) <- rownames(countMat)
  colnames(countMat2) <- colnames(countMat)
  
  metadata <- paste0(parentFolder, '/meta.csv')
  metadata <- utils::read.csv(metadata, header = T, row.names = 1)

  return(Seurat::CreateSeuratObject(counts = countMat2, project = paste0('RIRA-', version), meta.data = metadata))
}

