library(Seurat)
library(SeuratDisk)

cells <- LoadH5Seurat("data/cells_postprocessed.h5Seurat")

SSW <- function(object, def_assay = "RNA", n_feat = 2000, npcs = 50, dims_use = 1:50, n.neighbors = 30,
                res = 0.8){
  DefaultAssay(object = object) <- def_assay
  if (def_assay == "RNA"){
    object <- NormalizeData(object = object)
    object <- FindVariableFeatures(object = object, nfeatures = n_feat)
    object <- ScaleData(object = object)
    object <- RunPCA(object = object, npcs = npcs, reduction.name = "pca", reduction.key = "PCA_")
    object <- FindNeighbors(object = object, dims = dims_use, reduction = "pca")
    object <- FindClusters(object = object, resolution = res)
    object <- RunUMAP(object = object, n.neighbors = n.neighbors, dims = dims_use, reduction = "pca", reduction.name = "umap", reduction.key = "UMAP_")
  } else if (def_assay == "integrated"){
    object <- ScaleData(object = object)
    object <- RunPCA(object = object, npcs = npcs, reduction.name = "pca.int", reduction.key = "PCAint_")
    object <- FindNeighbors(object = object, dims = dims_use, reduction = "pca.int")
    object <- FindClusters(object = object, resolution = res)
    object <- RunUMAP(object = object, n.neighbors = n.neighbors, dims = dims_use, reduction = "pca.int", reduction.name = "umap.int", reduction.key = "UMAPint_")
  }
  return(object)
}

cells.lps <- subset(cells, subset = treat %in% c("Saline","LPS"))

hpcs.lps <- subset(cells.lps, subset = cell_type_brief %in% c("Som","Lac","Cort","Gonad","Mel","Thyro"))
hpcs.lps <- SSW(hpcs.lps, def_assay = "RNA")
hpcs.lps <- SSW(hpcs.lps, def_assay = "integrated")

resolutions <- seq(0.1, 1.5, 0.1)

som.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
# som.lps <- SSW(som.lps, def_assay = "RNA")
# som.lps <- SSW(som.lps, def_assay = "integrated")
DefaultAssay(som.lps) <- "RNA"
som.lps <- som.lps %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50)
res = 0.02; som.lps <- FindClusters(som.lps, resolution = res)  # a range of resolutions tested manually
nrow(unique(som.lps[[paste0("RNA_snn_res.",res)]]))



lac.lps <- subset(hpcs.lps, subset = cell_type_brief == "Lac")
lac.lps <- SSW(lac.lps, def_assay = "RNA")
lac.lps <- SSW(lac.lps, def_assay = "integrated")

cort.lps <- subset(hpcs.lps, subset = cell_type_brief == "Cort")
cort.lps <- SSW(cort.lps, def_assay = "RNA")
cort.lps <- SSW(cort.lps, def_assay = "integrated")

gonad.lps <- subset(hpcs.lps, subset = cell_type_brief == "Gonad")
gonad.lps <- SSW(gonad.lps, def_assay = "RNA")
gonad.lps <- SSW(gonad.lps, def_assay = "integrated")

mel.lps <- subset(hpcs.lps, subset = cell_type_brief == "Mel")
mel.lps <- SSW(mel.lps, def_assay = "RNA")
mel.lps <- SSW(mel.lps, def_assay = "integrated")

thyro.lps <- subset(hpcs.lps, subset = cell_type_brief == "Thyro")
thyro.lps <- SSW(thyro.lps, def_assay = "RNA")
thyro.lps <- SSW(thyro.lps, def_assay = "integrated")

