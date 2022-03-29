## https://smorabit.github.io/tutorials/8_velocyto/

library(SeuratDisk)

# assuming that you have some Seurat object called seurat_obj:
seurat_obj <- LoadH5Seurat('data/cells_postprocessed.h5Seurat')

# save metadata table:
seurat_obj$cell_id <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='data/scvelo/metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
writeMM(counts_matrix, file = 'data/scvelo/counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file='data/scvelo/pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='data/scvelo/gene_names.csv',
  quote=F,row.names=F,col.names=F
)

