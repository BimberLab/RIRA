# DietSeurat_Cust : a custom version of the DietSeurat() to really do the job

DietSeurat_Cust <- function(serobj = NULL, addCITE = F){
  serobj_temp = CreateSeuratObject(counts = serobj@assays$RNA@counts, 
                                   meta.data = serobj@meta.data)
  serobj_temp = SetAssayData(
    object = serobj_temp,
    assay.type = 'RNA',
    new.data = serobj@assays$RNA@data,
    slot = 'data'
  )
  
  if(addCITE){
    
    serobj_temp = SetAssayData(
      object = serobj_temp,
      assay.type = 'ADT',
      new.data = serobj@assays$ADT@counts,
      slot = 'counts'
    )
    serobj_temp = SetAssayData(
      object = serobj_temp,
      assay.type = 'ADT',
      new.data = serobj@assays$ADT@data,
      slot = 'data'
    )
    
    # can implement adding ADT tsne and umap
    
  }
  
  serobj_temp$rnaTSNE_1 = serobj@reductions$tsne@cell.embeddings[,1]
  serobj_temp$rnaTSNE_2 = serobj@reductions$tsne@cell.embeddings[,2]
  
  serobj_temp$rnaUMAP_1 = serobj@reductions$umap@cell.embeddings[,1]
  serobj_temp$rnaUMAP_2 = serobj@reductions$umap@cell.embeddings[,2]
  
  # can also implement adding PCAs
  
  
  return(serobj_temp)
}