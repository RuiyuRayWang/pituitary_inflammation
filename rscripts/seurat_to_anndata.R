## https://smorabit.github.io/tutorials/8_velocyto/

library(Seurat)
library(SeuratDisk)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

f_name = "hpcs_lps_state_marked"
f_path = "../data/"
out_dir = paste0(f_path,'tmp/')
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

# assuming that you have some Seurat object called seurat_obj:
seurat_obj <- LoadH5Seurat(paste0(f_path, f_name, ".h5Seurat"))

# save metadata table:
seurat_obj$cell_id <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
seurat_obj$UMAPint_1 <- seurat_obj@reductions$umap.int@cell.embeddings[,1]
seurat_obj$UMAPint_2 <- seurat_obj@reductions$umap.int@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file=paste0(out_dir,f_name,'_metadata.csv'), quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
writeMM(counts_matrix, file = paste0(out_dir,f_name,'_counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
for(pca in grep("pca", names(seurat_obj@reductions), value = T)){
  write.csv(seurat_obj@reductions[[pca]]@cell.embeddings, file=paste0(out_dir,f_name,'_',sub("\\.","",pca),'.csv'), quote=F, row.names=F)
}

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file=paste0(out_dir,f_name,'_gene_names.csv'),
  quote=F,row.names=F,col.names=F
)

