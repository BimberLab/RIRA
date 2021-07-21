
# classifier_Seurat : calls seurat fx to classify a query dataset 

classifier_Seurat <- function(){
 
  
  
  
  # outLS must contain appropriate output for classify_cells()
  outLS = list()
  return(outLS)
}





# classifier_SingleR : calls seurat fx to classify a query dataset 

classifier_SingleR <- function(){
  
  
  
  
  # outLS must contain appropriate output for classify_cells()
  outLS = list()
  return(outLS)
}



#' @title Modified Copy of CellMembrane's RunSingleR
#'
#' @description This copy allows for a local reference to be run, until we host the atlas remotely. Defaults to parallelization (8 cores) and the rhesusAtlas (to be provided from running ConvertAtlasAndApplyMetadata() ).
#' @param seuratObj The Seurat Object to be labeled by SingleR
#' @param datasets One or more datasets to use as a reference. Allowable values are: hpca, blueprint, dice, monaco, rhesusAtlas, and immgen. See cellDex package for available datasets.
#' @param assay The assay in the seurat object to use
#' @param resultTableFile If provided, a table of results will be saved here
#' @param rawDataFile If provided, the complete SingleR results will be saved to this file
#' @param minFraction If provided, any labels present with fraction of this or fewer across cells will be converted to Unknown
#' @param showHeatmap If true, heatmaps will be generated showing the SingleR calls
#' @param maxCellsForHeatmap The heatmap will only be plotted if the total cells is below this number
#' @param nThreads If provided, this integer value is passed to SingleR's BPPARAM argument. On windows ths is passed to BiocParallel::SnowParam(). On other OS it is passed to BiocParallel::MulticoreParam()
#' @return The modified seurat object
#' @importFrom pheatmap pheatmap
#' @import Seurat
#' @import SingleR
#' @export
#' @importFrom scuttle logNormCounts
#' @examples
#' seuratObj.test <- AtlasRunSingleR(seuratObj = seuratObj.test, datasets = c("rhesusAtlas"), sceAtlas = sceAtlas)
AtlasRunSingleR <- function(seuratObj = NULL, datasets = c("rhesusAtlas"), assay = NULL, resultTableFile = NULL, rawDataFile = NULL, minFraction = 0.01, showHeatmap = TRUE, maxCellsForHeatmap = 20000, nThreads = 8, sceAtlas = NULL){
  if (is.null(seuratObj)){
    stop("Seurat object is required")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seuratObj)
  }
  
  if (length(seuratObj@assays[[assay]]@counts) == 0) {
    print('Selected assay has no count data, trying RNA')
    assay <- 'RNA'
    if (length(seuratObj@assays[[assay]]@counts) == 0) {
      warning('Unable to find counts for the seurat object, aborting SingleR')
      return(seuratObj)
    }
  }
  
  allFields <- c()
  completeRawData <- NULL
  
  for (dataset in datasets) {
    print(paste0('Adding dataset: ', dataset))
    if (dataset == 'hpca'){
      ref <- celldex::HumanPrimaryCellAtlasData()
    } else if (dataset == 'immgen') {
      ref <- celldex::ImmGenData()
    } else if (dataset == 'blueprint') {
      ref <- celldex::BlueprintEncodeData()
    } else if (dataset == 'dice') {
      ref <- celldex::DatabaseImmuneCellExpressionData()
    } else if (dataset == 'monaco') {
      ref <- celldex::MonacoImmuneData()
    } else if (dataset == 'rhesusAtlas'){
      ref <- sceAtlas
    } else {
      stop(paste0('unknown reference dataset: ', dataset))
    }
    
    #Subset genes:
    genesPresent <- intersect(rownames(seuratObj@assays[[assay]]), rownames(ref))
    print(paste0('Total genes shared with reference data: ', length(genesPresent), ' of ', nrow(seuratObj)))
    
    if (length(genesPresent) < 100) {
      print(paste0('Too few shared genes, skipping: ', length(genesPresent)))
      next
    }
    
    ref <- ref[genesPresent,]
    
    seuratObjSubset <- Seurat::DietSeurat(seuratObj, assays = c(assay), counts = TRUE)
    seuratObjSubset <- subset(seuratObj, features = genesPresent)
    Seurat::DefaultAssay(seuratObjSubset) <- assay
    
    #Convert to SingleCellExperiment
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = GetAssayData(object = seuratObjSubset, assay = assay, slot = "counts")))
    sce <- scuttle::logNormCounts(sce)
    rm(seuratObjSubset)
    
    refAssay <- 'logcounts'
    if (!('logcounts' %in% names(SummarizedExperiment::assays(ref)))) {
      print('logcount not present, using normcounts as assay type')
      refAssay <- 'normcounts'
    }
    
    if (!is.null(nThreads)) {
      if (.Platform$OS.type == 'windows') {
        BPPARAM <- BiocParallel::SnowParam(nThreads)
      } else {
        BPPARAM <- BiocParallel::MulticoreParam(nThreads)
      }
    } else {
      BPPARAM <- BiocParallel::SerialParam()
    }
    
    tryCatch({
      pred.results <- suppressWarnings(SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main, assay.type.test = 'logcounts', assay.type.ref = refAssay, fine.tune = TRUE, prune = TRUE, BPPARAM = BPPARAM))
      if (length(colnames(seuratObj)) != nrow(pred.results)) {
        stop('Length of SingleR results did not match seurat object')
      }
      
      if (!is.null(rawDataFile)){
        toBind <- data.frame(cellbarcode = colnames(seuratObj), classification_type = 'Main', dataset = dataset, labels = pred.results$labels, pruned.labels = pred.results$pruned.labels)
        if (is.null(completeRawData)) {
          completeRawData <- toBind
        } else {
          completeRawData <- rbind(completeRawData, toBind)
        }
      }
      
      if (showHeatmap) {
        cells.use <- NULL
        if (ncol(seuratObj) > maxCellsForHeatmap) {
          cells.use <- sample(1:ncol(seuratObj), size = maxCellsForHeatmap)
        }
        
        print(SingleR::plotScoreHeatmap(pred.results, cells.use = cells.use))
      }
      
      toAdd <- pred.results$pruned.labels
      toAdd[is.na(toAdd)] <- 'Unknown'
      names(toAdd) <- colnames(seuratObj)
      fn <- paste0(dataset, '.label')
      allFields <- c(allFields, fn)
      seuratObj[[fn]] <- toAdd
      
      pred.results <- suppressWarnings(SingleR::SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.test = 'logcounts', assay.type.ref = refAssay, fine.tune = TRUE, prune = TRUE))
      if (length(colnames(seuratObj)) != nrow(pred.results)) {
        stop('Length of SingleR results did not match seurat object')
      }
      
      if (!is.null(rawDataFile)){
        toBind <- data.frame(cellbarcode = colnames(seuratObj), classification_type = 'Fine', dataset = dataset, labels = pred.results$labels, pruned.labels = pred.results$pruned.labels)
        if (is.null(completeRawData)) {
          completeRawData <- toBind
        } else {
          completeRawData <- rbind(completeRawData, toBind)
        }
      }
      
      if (showHeatmap) {
        cells.use <- NULL
        if (ncol(seuratObj) > maxCellsForHeatmap) {
          cells.use <- sample(1:ncol(seuratObj), size = maxCellsForHeatmap)
        }
        
        print(SingleR::plotScoreHeatmap(pred.results, cells.use = cells.use))
      }
      
      toAdd <- pred.results$pruned.labels
      toAdd[is.na(toAdd)] <- 'Unknown'
      names(toAdd) <- colnames(seuratObj)
      
      fn2 <- paste0(dataset, '.label.fine')
      allFields <- c(allFields, fn2)
      seuratObj[[fn2]] <- toAdd
      
      if (!is.null(minFraction)){
        for (label in c(fn, fn2)) {
          print(paste0('Filtering ', label, ' below: ', minFraction))
          d <- data.frame(table(Label = unlist(seuratObj[[label]])))
          names(d) <- c('Label', 'Count')
          d$Fraction <- d$Count / sum(d$Count)
          
          d <- d %>% dplyr::arrange(dplyr::desc(Fraction))
          print(d)
          toRemove <- d$Label[d$Fraction < minFraction]
          if (length(toRemove) > 0) {
            print(paste0('Will remove: ', paste0(toRemove, collapse = ', ')))
          }
          
          if (length(toRemove) > 0) {
            l <- unlist(seuratObj[[label]])
            names(l) <- colnames(seuratObj)
            l[l %in% toRemove] <- 'Unknown'
            seuratObj[[label]] <- l
          }
          
          print('After filter:')
          d <- data.frame(table(Label = unlist(seuratObj[[label]])))
          names(d) <- c('Label', 'Count')
          print(d)
        }
      }
    }, error = function(e){
      print(paste0('Error running singleR for dataset: ', dataset))
      print(conditionMessage(e))
      traceback()
    })
  }
  
  print(paste0('Adding fields: ', paste0(allFields, collapse = ',')))
  df <- data.frame(CellBarcodes = colnames(seuratObj))
  for (fn in allFields){
    df[fn] <- seuratObj@meta.data[[fn]]
  }
  
  if (!is.null(rawDataFile)) {
    write.table(completeRawData, file = rawDataFile, sep = '\t', row.names = FALSE)
  }
  
  if (!is.null(resultTableFile)){
    write.table(file = resultTableFile, df, sep = '\t', row.names = F, quote = F)
  }
  
  if (length(allFields) == 0) {
    print('No singleR calls were added, this probably indicates there were errors with singleR')
  } else {
    if ('tsne' %in% names(seuratObj@reductions) || 'umap' %in% names(seuratObj@reductions)) {
      DimPlot_SingleR(seuratObj, plotIndividually = TRUE, datasets = datasets)
    }
    
    Tabulate_SingleR(seuratObj, plotIndividually = TRUE, datasets = datasets)
  }
  
  return(seuratObj)
}

#' @title Convert a Seurat object atlas into a Single Cell Experiment type Atlas and apply relevant metadata
#' 
#' @param seuratAtlas A Seurat object type Rhesus Atlas
#' @param label.main The Pheno.scope that is used for the highest singleR scope. Default = Pheno.scope1
#' @param label.fine The Pheno.scope that is used for the lower singleR scope. Default = Pheno.scope4
#' @return A SCE type atlas, ready to be used in SingleR functions. 
#' @examples
#' sceAtlas <- ConvertAtlasAndApplyMetadata(seuratAtlas)
#' sceAtlas <- ConvertAtlasAndApplyMetadata(seuratAtlas, label.main = "Pheno.scope2", label.fine = "Pheno.scope6")
#' @import SingleR
#' @import SingleCellExperiment
#' @import future
#' @import doFuture
#' @export
#' @importFrom scuttle logNormCounts
ConvertAtlasAndApplyMetadata <- function(seuratAtlas = NULL, label.main = "Pheno.scope1", label.fine = "Pheno.scope4"){
  require(doFuture)
  atlasObj <- GetAtlasMetadata()
  #Convert and lognorm the Atlas (seuratAtlas) 
  sceAtlas <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = GetAssayData(object = seuratAtlas, assay = "RNA", slot = "counts")))
  sceAtlas <- scuttle::logNormCounts(sceAtlas)
  
  #Mapping between the Atlas and metadata takes awhile, this parallelization scheme speeds it up a bit
  doFuture::registerDoFuture()
  future::plan(multiprocess)
  
  #apply label.main (Defaults to Pheno.scope1)
  labels <- foreach(cell_name = rownames(sceAtlas@colData), .combine = 'c') %dopar% {
    label.main <- atlasObj[cell_name,"Pheno.scope1"]
    label.main
  }
  sceAtlas@colData$label.main <- labels
  
  #apply label.fine (Defaults to Pheno.scope4)
  labels <- foreach(cell_name = rownames(sceAtlas@colData), .combine = 'c') %dopar% {
    label.fine <- atlasObj.All[cell_name,"Pheno.scope4"]
    label.fine
  }
  sceAtlas@colData$label.fine <- labels
  return(sceAtlas)
}

#' @title Helper function for ConvertAtlasAndApplyMetadata. 
#' 
#' @description This needs amending in the future, but for now it grabs the metadata on Box
#' @return A data frame with the current atlas metadata. 
#' @examples
#' atlasObj <- GetAtlasMetadata()
GetAtlasMetadata <- function(){
  metaBcell <- '~/Box/Bimber/BimberLab Shared/Experiments/RhesusAtlas/data/meta/June232021_freeze/RhesusAtlasV2.Meta.Bcell_June232021.rds'
  atlasObj.Bcell <- readRDS(file = metaBcell)
  metaCD8pos <- '~/Box/Bimber/BimberLab Shared/Experiments/RhesusAtlas/data/meta/June232021_freeze/RhesusAtlasV2.Meta.Tcell.CD8pos_June232021.rds'
  atlasObj.CD8pos <- readRDS(file = metaCD8pos)
  metaMcell <- '~/Box/Bimber/BimberLab Shared/Experiments/RhesusAtlas/data/meta/June232021_freeze/RhesusAtlasV2.Meta.Mcell_June232021.rds'
  atlasObj.Mcell <- readRDS(file = metaMcell)
  metaCD4pos <- '~/Box/Bimber/BimberLab Shared/Experiments/RhesusAtlas/data/meta/June232021_freeze/RhesusAtlasV2.Meta.Tcell.CD4pos_June232021.rds'
  atlasObj.CD4pos <- readRDS(file = metaCD4pos)
  
  metaDblNeg <- '~/Box/Bimber/BimberLab Shared/Experiments/RhesusAtlas/data/meta/June232021_freeze/RhesusAtlasV2.Meta.Tcell.DblNeg_June232021.rds'
  atlasObj.DblNeg <- readRDS(file = metaDblNeg)
  
  metaDblts <- '~/Box/Bimber/BimberLab Shared/Experiments/RhesusAtlas/data/meta/June232021_freeze/RhesusAtlasV2.Meta.Dblts_June232021.rds'
  atlasObj.Dblts <- readRDS(file = metaDblts)
  
  metaDblPos <- '~/Box/Bimber/BimberLab Shared/Experiments/RhesusAtlas/data/meta/June232021_freeze/RhesusAtlasV2.Meta.Tcell.DblPos_June232021.rds'
  atlasObj.DblPos <- readRDS(file = metaDblPos)
  
  metaTriNeg <- '~/Box/Bimber/BimberLab Shared/Experiments/RhesusAtlas/data/meta/June232021_freeze/RhesusAtlasV2.Meta.TriNeg_June232021.rds'
  atlasObj.TriNeg <- readRDS(file = metaTriNeg)
  
  atlasObj.All <- bind_rows(list(atlasObj.Bcell, atlasObj.CD4pos, atlasObj.CD8pos, atlasObj.Mcell,atlasObj.DblNeg ,atlasObj.DblPos,atlasObj.TriNeg, atlasObj.Dblts))
  
  return(atlasObj.All)
}
