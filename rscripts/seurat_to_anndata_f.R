#' @title
#' Seurat to Anndata (R part)
#' 
#' @description
#' Convert a .h5Seurat Seurat object into anndata. Intermediate 
#' results are dumped to temporary directory `/tmp`.
#' 
#' @details
#' 
#' @param f_name File name of the Seurat object. The .h5Seurat extension should be omitted.
#' @param f_path File path to the Seurat object.
#' @param out_dir Directory for output of the intermediate files. A folder named `/tmp` will be created.
#' @param embeddings Character vector specifying the embeddings of non-linear dimension reductions to 
#' save, i.e. c("tsne","umap").
#' @param slot Slot to save. By default the count matrix, i.e. "counts".
#' @param dr Character vector specifying the dimension reduction to save, by default "pca".
#' 
#' @author Sam Morabito
#' @references https://smorabit.github.io/tutorials/8_velocyto/
#' 
SeuratToAnndata <- function(
    f_name,
    f_path,
    out_dir,
    embeddings = NULL,
    slot = "counts",
    dr = "pca"
){
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
  
  # assuming that you have some Seurat object called seurat_obj:
  seurat_obj <- LoadH5Seurat(paste0(f_path, f_name, ".h5Seurat"))
  
  # save metadata table:
  seurat_obj[['cell_id']] <- colnames(seurat_obj)
  
  seurat_obj[[""]] <- seurat_obj@reductions$umap@cell.embeddings[,1]
  seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
  seurat_obj$UMAPint_1 <- seurat_obj@reductions$umap.int@cell.embeddings[,1]
  seurat_obj$UMAPint_2 <- seurat_obj@reductions$umap.int@cell.embeddings[,2]
  write.csv(seurat_obj@meta.data, file=paste0(out_dir,f_name,'_metadata.csv'), quote=F, row.names=F)
  
  
}