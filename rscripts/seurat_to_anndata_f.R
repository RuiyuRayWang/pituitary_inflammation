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
#' @param assay Assay to get expression data, by default "RNA".
#' @param slot Slot to save. By default the count matrix, i.e. "counts".
#' @param dr Character vector specifying the dimension reduction to save, by default "pca".
#' 
#' @author Sam Morabito
#' @references https://smorabit.github.io/tutorials/8_velocyto/
#' 
#' @importFrom SeuratDisk Loadh5Seurat
#' @importFrom Seurat Embeddings GetAssayData
#' @importFrom Matrix writeMM
#' 
SeuratToAnndata <- function(
    f_name,
    f_path,
    out_dir,
    embeddings = NULL,
    assay = "RNA",
    slot = "counts",
    dim.red = "pca"
){
  out_dir = file.path(out_dir,"tmp")
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
  embeddings = embeddings %||% "umap"
  
  # assuming that you have some Seurat object called seurat_obj:
  seurat_obj <- LoadH5Seurat(paste0(f_path, f_name, ".h5Seurat"))
  
  # save metadata table
  seurat_obj[['cell_id']] <- colnames(seurat_obj)
  write.csv(seurat_obj@meta.data, file=file.path(out_dir,paste0(f_name,'_metadata.csv')), quote=F, row.names=F)
  
  # save embeddings
  lapply(X = embeddings, FUN = function(x){
    write.csv(
      x = Embeddings(object = seurat_obj, reduction = x), 
      file = file.path(out_dir,paste0('embeddings_',x,'.csv')),
      quote = F, row.names = F
      )
  })
  # emb <- do.call(what = cbind, args = emb)
  # write.csv(emb, file=file.path(out_dir,paste0(f_name,'_embeddings.csv')), quote=F, row.names=F)
  
  # write expression counts matrix
  counts_mtx <- GetAssayData(object = seurat_obj, assay = assay, slot = slot)
  writeMM(obj = counts_mtx, file = file.path(out_dir, paste0(f_name,'_counts.mtx')))
  
  # write dimensionality reduction matrix
  lapply(X = dim.red, FUN = function(x){
    write.csv(
      x = Embeddings(object = seurat_obj, reduction = x), 
      file = file.path(out_dir,paste0('dr_',x,'.csv')),
      quote = F, row.names = F
    )
  })
  # dr <- do.call(what = cbind, args = dr)
  # write.csv(dr, file=file.path(out_dir,paste0(f_name,'_dr.csv')), quote=F, row.names=F)
}