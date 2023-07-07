getBaseSeuratData <- function(){
  # See: https://github.com/satijalab/seurat-data/issues/70#issuecomment-1625806697
  install.packages("https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.0.0.tar.gz", repos = NULL, type = "source")

  suppressWarnings(SeuratData::InstallData("pbmc3k"))
  suppressWarnings(data("pbmc3k"))
  seuratObj <- suppressWarnings(pbmc3k)

  return(seuratObj)
}