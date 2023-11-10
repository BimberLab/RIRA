getBaseSeuratData <- function(){
  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(Seurat::UpdateSeuratObject(pbmc3k))

  return(seuratObj)
}