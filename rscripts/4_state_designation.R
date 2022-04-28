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

Idents(hpcs.lps) <- "Undetermined"

## Iteratively determine cell states
## Repeat same procedures for all major hormonal cell types
som.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
# som.lps <- SSW(som.lps, def_assay = "RNA")
# som.lps <- SSW(som.lps, def_assay = "integrated")
DefaultAssay(som.lps) <- "RNA"
som.lps <- som.lps %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50)
res = 0.02; som.lps <- FindClusters(som.lps, resolution = res)  ## A range of resolutions tested manually
nrow(unique(som.lps[[paste0("RNA_snn_res.",res)]])) == 2
Idents(som.lps) <- paste0("RNA_snn_res.",res)
som.lps <- RenameIdents(som.lps, `0` = "Inflammation", `1` = "Healthy")
som.lps$state <- Idents(som.lps)
for (s in levels(som.lps$state)){
  hpcs.lps <- SetIdent(hpcs.lps, cells = WhichCells(som.lps, expression = state == s), value = s)  ## Copy labels to master object (hpcs.lps)
}

lac.lps <- subset(hpcs.lps, subset = cell_type_brief == "Lac")
DefaultAssay(lac.lps) <- "RNA"
lac.lps <- lac.lps %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50)
res = 0.1; lac.lps <- FindClusters(lac.lps, resolution = res)  # a range of resolutions tested manually
nrow(unique(lac.lps[[paste0("RNA_snn_res.",res)]])) == 2
Idents(lac.lps) <- paste0("RNA_snn_res.",res)
# DimPlot(lac.lps) | DimPlot(cort.lps, group.by = "stim")
lac.lps <- RenameIdents(lac.lps, `0` = "Healthy", `1` = "Inflammation")
lac.lps$state <- Idents(lac.lps)
for (s in levels(lac.lps$state)){
  hpcs.lps <- SetIdent(hpcs.lps, cells = WhichCells(lac.lps, expression = state == s), value = s)
}

cort.lps <- subset(hpcs.lps, subset = cell_type_brief == "Cort")
DefaultAssay(cort.lps) <- "RNA"
cort.lps <- cort.lps %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50)
res = 0.1; cort.lps <- FindClusters(cort.lps, resolution = res)  # a range of resolutions tested manually
nrow(unique(cort.lps[[paste0("RNA_snn_res.",res)]])) == 2
Idents(cort.lps) <- paste0("RNA_snn_res.",res)
# DimPlot(cort.lps) | DimPlot(cort.lps, group.by = "stim")
cort.lps <- RenameIdents(cort.lps, `0` = "Healthy", `1` = "Inflammation")
cort.lps$state <- Idents(cort.lps)
for (s in levels(cort.lps$state)){
  hpcs.lps <- SetIdent(hpcs.lps, cells = WhichCells(cort.lps, expression = state == s), value = s)
}

gonad.lps <- subset(hpcs.lps, subset = cell_type_brief == "Gonad")
DefaultAssay(gonad.lps) <- "RNA"
gonad.lps <- gonad.lps %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50)
res = 0.1; gonad.lps <- FindClusters(gonad.lps, resolution = res)  # a range of resolutions tested manually
nrow(unique(gonad.lps[[paste0("RNA_snn_res.",res)]])) == 2
Idents(gonad.lps) <- paste0("RNA_snn_res.",res)
# DimPlot(gonad.lps) | DimPlot(gonad.lps, group.by = "stim")
gonad.lps <- RenameIdents(gonad.lps, `0` = "Healthy", `1` = "Inflammation")
gonad.lps$state <- Idents(gonad.lps)
for (s in levels(gonad.lps$state)){
  hpcs.lps <- SetIdent(hpcs.lps, cells = WhichCells(gonad.lps, expression = state == s), value = s)
}

mel.lps <- subset(hpcs.lps, subset = cell_type_brief == "Mel")
DefaultAssay(mel.lps) <- "RNA"
mel.lps <- mel.lps %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50)
res = 0.1; mel.lps <- FindClusters(mel.lps, resolution = res)  # a range of resolutions tested manually
nrow(unique(mel.lps[[paste0("RNA_snn_res.",res)]])) == 2
Idents(mel.lps) <- paste0("RNA_snn_res.",res)
# DimPlot(mel.lps) | DimPlot(mel.lps, group.by = "stim")
mel.lps <- RenameIdents(mel.lps, `1` = "Healthy", `0` = "Inflammation")
mel.lps$state <- Idents(mel.lps)
for (s in levels(mel.lps$state)){
  hpcs.lps <- SetIdent(hpcs.lps, cells = WhichCells(mel.lps, expression = state == s), value = s)
}

thyro.lps <- subset(hpcs.lps, subset = cell_type_brief == "Thyro")
DefaultAssay(thyro.lps) <- "RNA"
thyro.lps <- thyro.lps %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:50) %>% RunUMAP(dims = 1:50)
res = 0.8; thyro.lps <- FindClusters(thyro.lps, resolution = res)  # a range of resolutions tested manually
nrow(unique(thyro.lps[[paste0("RNA_snn_res.",res)]])) == 2
Idents(thyro.lps) <- paste0("RNA_snn_res.",res)
# DimPlot(thyro.lps) | DimPlot(thyro.lps, group.by = "stim")
thyro.lps <- RenameIdents(thyro.lps, `0` = "Healthy", `1` = "Inflammation")
thyro.lps$state <- Idents(thyro.lps)
for (s in levels(thyro.lps$state)){
  hpcs.lps <- SetIdent(hpcs.lps, cells = WhichCells(thyro.lps, expression = state == s), value = s)
}

SaveH5Seurat(hpcs.lps, "data/hpcs_state_marked.h5Seurat", overwrite = T, verbose = F)
Convert(source = "data/hpcs_state_marked.h5Seurat", dest = "h5ad", overwrite = T)
