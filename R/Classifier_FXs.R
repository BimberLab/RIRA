
# classifier_Seurat : calls seurat fx to classify a query dataset 

getSeuratAnchors <- function(ref_seurat, target_seurat) {
  anchors <- FindTransferAnchors(reference = ref_seurat, query = target_seurat, dims = 1:30)
  return(anchors)
}

classifier_Seurat <- function(anchors, ref_metalabels){
  predictions <- TransferData(anchorset = anchors, refdata = ref_metalabels, 
                               dims = 1:30)
  
  return(predictions)
}



# =============================================================================
# classifier_cFIT : calls cFIT fx to classify a query dataset 

# split reference data by batch (ref_meta_batchID = column in metadata with batch IDs--e.g., experimentID)
# also can use this function for changing labels for transfer (prior to asign_labels)
format_refdata <- function(ref_seurat, ref_meta_batchID, ref_meta) {
  data.list = cFIT::split_dataset_by_batch(X=t(as.matrix(ref_seurat@assays$RNA@counts)),
                                           batch = ref_meta_batchID, 
                                           labels = ref_meta_labels, 
                                           metadata = ref_meta, 
                                           dataset.name = 'atlas:')
  return(data.list)
}

# select highly variable genes in reference
genes = cFIT::select_genes(data.list$X.list, ngenes=2000, verbose=F)

# perform integration on reference
integrate_ref <- function(genes, data.list) {
  # reference data preprocessing
  exprs.list = cFIT::preprocess_for_integration(data.list$X.list, genes, scale.factor=10^4, scale=T, center=F)
  # perform reference data integration across batches (longest step)
  int.out = cFIT::CFITIntegrate(X.list=exprs.list, r=15, verbose=F, max.niter = 100, seed=0)
  return(int.out)
}

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

#assign labels from reference to the target data
est.labels = cFIT::asign_labels(exprs.source=do.call(rbind, int.out$H.list), 
                                exprs.target=tf.out$H, 
                                labels.source=do.call(c, data.list$labels.list))

# =============================================================================




# classifier_SingleR : calls seurat fx to classify a query dataset 

classifier_SingleR <- function(){
  
  
  
  
  # outLS must contain appropriate output for classify_cells()
  outLS = list()
  return(outLS)
}

