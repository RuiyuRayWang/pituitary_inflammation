library(Seurat)
library(SeuratDisk)
library(tidyverse)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

cells <- LoadH5Seurat("../data/processed/cells_postprocessed.h5Seurat")
hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")

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

nonhpcs.lps <- subset(cells.lps, subset = cell_type_brief %in% c("Endo","Peri","Pou1f1","Stem","WBCs","RBCs"))
nonhpcs.lps <- SSW(nonhpcs.lps, def_assay = "RNA")
nonhpcs.lps <- SSW(nonhpcs.lps, def_assay = "integrated")

Idents(nonhpcs.lps) <- "Undetermined"

## Iteratively determine cell states
## Repeat same procedures for all non-hormonal cell types
endo.lps <- subset(nonhpcs.lps, subset = cell_type_brief == "Endo")
# endo.lps <- SSW(endo.lps, def_assay = "RNA")
# endo.lps <- SSW(endo.lps, def_assay = "integrated")
DefaultAssay(endo.lps) <- "RNA"
endo.lps <- endo.lps %>%
  NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% FindNeighbors(dims = 1:20) %>% RunUMAP(dims = 1:20)
res = 0.4; endo.lps <- FindClusters(endo.lps, resolution = res)  ## A range of resolutions tested manually
nrow(unique(endo.lps[[paste0("RNA_snn_res.",res)]])) == 2
Idents(endo.lps) <- paste0("RNA_snn_res.",res)
endo.lps <- RenameIdents(endo.lps, `0` = "Inflammation", `1` = "Healthy")
endo.lps$state <- Idents(endo.lps)
for (s in levels(endo.lps$state)){
  nonhpcs.lps <- SetIdent(nonhpcs.lps, cells = WhichCells(endo.lps, expression = state == s), value = s)  ## Copy labels to master object (nonhpcs.lps)
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

hpcs.lps$state <- Idents(hpcs.lps)

SaveH5Seurat(hpcs.lps, "../data/processed/hpcs_lps_state_marked.h5Seurat", overwrite = T, verbose = F)
# Convert(source = "../data/processed/hpcs_lps_state_marked.h5Seurat", dest = "h5ad", overwrite = T, verbose = F)  ## DO NOT RUN, use seurat_to_anndata.R

## Copy this file to 'scenic_protocol/files' for SCENIC analyses
dst_dir = "scenic_protocol/files/"
if(!dir.exists(dst_dir)){dir.create(dst_dir, recursive = T)}
f_loom = c('data/hpcs_lps.loom', file.path(dst_dir, 'hpcs_lps.loom'))
if(any(file.exists(f_loom))){
  ## Overwriting files causes pyscenic to fail
  file.remove(f_loom)
}
SaveLoom(hpcs.lps, filename = f_loom[1])
file.copy(from = f_loom[1],
          to = f_loom[2])
