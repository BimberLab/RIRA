


######### Seurat


#' @title integrate_reflist_CCA 
#'
#' @description integrate CCA given a list of ref Seurat object see SplitObject() of Seurat 
#' @param RefLS The seurat object ref list see SplitObject() 
#' @param dims CCA/PCA dims to integrate default 1:30 
#' @param dims_umap CCA/PCA dims to run UMAP default 1:20 
#' @return list of 2 times Ref.anchors and Ref.integrated
#' @export
integrate_reflist_CCA <- function(RefLS = NULL, dims = 1:30, dims_umap=1:20){
  
  Ref.anchors <- FindIntegrationAnchors(object.list = RefLS, dims = dims)
  Ref.integrated <- IntegrateData(anchorset = Ref.anchors, dims = dims)
  DefaultAssay(Ref.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  Ref.integrated <- ScaleData(Ref.integrated, verbose = FALSE)
  Ref.integrated <- RunPCA(Ref.integrated, npcs = 50, verbose = FALSE)
  ElbowPlot(Ref.integrated)
  Ref.integrated <- RunUMAP(Ref.integrated, reduction = "pca", dims = dims_umap, verbose = FALSE)
  return(list(Ref.anchors = Ref.anchors, Ref.integrated = Ref.integrated))
}





#' @title classifier_Seurat 
#'
#' @description calls seurat fx to classify a query dataset 
#' @param Ref.integrated The seurat object of the Ref.integrated
#' @param Query.Ser The seurat object of the Query
#' @param dims CCA/PCA dims to integrate default 1:30 
#' @param dims_umap CCA/PCA dims to run UMAP default 1:20 
#' @param Ref.meta metadata pertaining to Ref
#' @export
classifier_Seurat <- function(Ref.integrated = NULL, Query.Ser=NULL, dims = 1:30, dims_umap=1:20,
                              Ref.meta){
  Query_anchors <- FindTransferAnchors(reference = Ref.integrated, query = Query.Ser, dims = dims)
  predictions <- TransferData(anchorset = Query_anchors, refdata = Ref.meta, dims = dims)
  Query.Ser <- AddMetaData(Query.Ser, metadata = predictions)
  
  Ref.integrated <- RunUMAP(Ref.integrated, dims = dims_umap, reduction = "pca", return.model = TRUE)
  
  Query.Ser <- MapQuery(anchorset = Query_anchors,
                           reference = Ref.integrated,
                           query = Query.Ser,
                           refdata = Ref.integrated$Pheno.scope1,
                           reference.reduction = "pca", 
                          reduction.model = "umap")
  
  
 
  return(Query.Ser)
}





# classifier_SingleR : calls seurat fx to classify a query dataset 

classifier_SingleR <- function(){
  
  
  
  
  # outLS must contain appropriate output for classify_cells()
  outLS = list()
  return(outLS)
}







#### cFIT
# devtools::install_github("pengminshi/cFIT")
# require(cFIT)

# =============================================================================
# classifier_cFIT : calls cFIT fx to classify a query dataset 
# split reference data by batch (ref_meta_batchID = column in metadata with batch IDs--e.g., experimentID)
# also can use this function for changing labels for transfer (prior to asign_labels)

#' @title classifier_Seurat 
#'
#' @description format_refdata takes seurat objects and returns appropriate list format
#' @param ref_seurat A Seurat object as reference
#' @param ref_meta_batchID a vector that defines splits/batches in the ref_seurat e.g. ref_seurat$myID
#' @param ref_meta_labels(optional) charateristic vector of cell labels (default NULL)
#' @param ref_meta the metadata of the ref object
#' @param dataset.name A name defaults 'atlas:'
#' @export
format_refdata <- function(ref_seurat, ref_meta_batchID, ref_meta = NULL, ref_meta_labels = NULL, dataset.name = 'atlas:') {
  
  if(is.null(ref_meta)) ref_meta = ref_seurat@meta.data
  if(is.null(ref_meta_labels)) stop("Need some labels! ad a list or DF to ref_meta_labels")
  
  data.list = cFIT::split_dataset_by_batch(X=t(as.matrix(ref_seurat@assays$RNA@counts)),
                                           batch = ref_meta_batchID, 
                                           labels = ref_meta_labels, 
                                           metadata = ref_meta, 
                                           dataset.name = dataset.name)
  return(data.list)
}





#' @title integrate_ref  perform integration on reference
#'
#' @description perform integration on reference cFIT
#' @param genes A set of genes usually top variable genes
#' @param data.list the output of format_refdata()
#' @export
integrate_ref <- function(genes, data.list) {
  # reference data preprocessing
  exprs.list = cFIT::preprocess_for_integration(data.list$X.list, genes, scale.factor=10^4, scale=T, center=F)
  # perform reference data integration across batches (longest step)
  int.out = cFIT::CFITIntegrate(X.list=exprs.list, r=15, verbose=F, max.niter = 100, seed=0)
  return(int.out)
}

#' @title transfer_sigs  perform CFITTransfer
#'
#' @description perform CFITTransfer
#' @param genes A set of genes usually top variable genes
#' @param target_seurat A query seurat object
#' @param int.out the output of integrate_ref()
#' @export
transfer_sigs <- function(genes, target_seurat, int.out) {
  #prepare target data for transfer
  target.X <- t(as.matrix(target_seurat@assays$RNA@counts))
  genes_target = intersect(genes, colnames(target.X))
  x = target.X[, genes_target]
  scale.factor=10^4
  x = log(x/rowSums(x) * scale.factor + 1)
  x = t(Seurat::ScaleData(t(x), do.center = F, do.scale = T, verbose = F))
  target.exprs.list <- x
  names(target.exprs.list) = names(target.X)
  #transfer learned gene expression signatures from integration onto target dataset 
  tf.out = cFIT::CFITTransfer(Xtarget=target.exprs.list, Wref=int.out$W, max.niter=100, seed=0, verbose=T)
  return(tf.out)
}

