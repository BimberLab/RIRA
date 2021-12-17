
pkg.env <- new.env(parent=emptyenv())

pkg.env$RANDOM_SEED <- 1234
set.seed(pkg.env$RANDOM_SEED)

#' @title Set random seed
#'
#' @description Sets the seed used for R‘s random number generator, which should be used in all internal functions
#' @param seed The random seed
#' @export
SetSeed <- function(seed) {
  pkg.env$RANDOM_SEED <- seed
  set.seed(pkg.env$RANDOM_SEED)
}

#' @title Get random seed
#'
#' @description Sets a random seed, which should be used in all internal functions
#' @export
GetSeed <- function() {
  return(pkg.env$RANDOM_SEED)
}

#' @title Set the directory where atlas data reside
#'
#' @description Sets the directory where atlas data reside.
#' @param folderPath The atlas directory. Within the folder, there should be sub-folders by version.
#' @export
SetAtlasDir <- function(folderPath) {
  if (!dir.exists(folderPath)) {
    stop(paste0('Unable to find folder: ', folderPath))
  }

  folderPath <- gsub(folderPath, pattern = '[/\\]+$', replacement = '')
  pkg.env$ATLAS_DIR <- folderPath
}

.GetAtlasDir <- function() {
  return(pkg.env$ATLAS_DIR)
}

.GetLatestVersion <- function(){
  return(utils::packageVersion('RIRA'))
}

.GetAtlasBaseDir <- function(version) {
  if (is.null(version)) {
    version <- .GetLatestVersion()
  }

  if (is.null(.GetAtlasDir())) {
    parentFolder <- system.file('inst/data', package = 'RIRA')
  } else {
    parentFolder <- .GetAtlasDir()
  }

  parentFolder <- paste0(parentFolder, '/', version)
  if (!dir.exists(parentFolder)) {
    print(paste0('Unknown version: ', version))
  }

  return(parentFolder)
}