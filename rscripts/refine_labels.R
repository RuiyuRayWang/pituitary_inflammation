library(Seurat)
library(SeuratDisk)
library(tidyverse)

cells <- LoadH5Seurat('data/cells_postprocessed.h5Seurat')

# Subsetting out hormone producing cells (HPCs)
cells.hpcs <- subset(cells, subset = cell_type_brief %in% c("Som","Cort","Gonad","Lac","Thyro","Mel"))

DefaultAssay(cells.hpcs) <- "RNA"
cells.hpcs <- cells.hpcs %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50)

d <- DimPlot(cells.hpcs, reduction = "umap")

DefaultAssay(cells.hpcs) <- "integrated"
cells.hpcs <- cells.hpcs %>%
  ScaleData() %>%
  RunPCA(reduction.name = "pca.int", reduction.key = "PCAint_") %>%
  RunUMAP(dims = 1:50, reduction = "pca.int", reduction.name = "umap.int", reduction.key = "UMAPint_")

d_int <- DimPlot(cells.hpcs, reduction = "umap.int")

## Iteratively query cell identities in integrated/un-integrated umaps, and refine their labels
DefaultAssay(cells.hpcs) <- "RNA"
HoverLocator(d, information = FetchData(cells.hpcs, vars = c("cell_type_brief","stim","Ghrhr","Pomc","Prl","Fshb","Lhb","Tshb")))

cell.id <- "CAAGGAGC_YT19070303"
DimPlot(cells.hpcs, cells.highlight = cell.id, reduction = "umap.int")

# Idents(cells.hpcs) <- "cell_type_brief"
cells.hpcs <- SetIdent(cells.hpcs, cells = cell.id, value = "Som")
d <- DimPlot(cells.hpcs, reduction = "umap")

## Store refined labels
cells.hpcs$cell_type_brief <- Idents(cells.hpcs)
x <- cells.hpcs$cell_type_brief
save(x, file = "misc/cells_hpcs_manual_refined_labels.rda")

SaveH5Seurat(cells.hpcs, filename = "data/cells_hpcs.h5Seurat", overwrite = T)

# LPS groups only
hpcs.lps <- subset(cells.hpcs, subset = treat %in% c("Saline","LPS"))

# Redo integration on subsetted object
DefaultAssay(hpcs.lps) <- "RNA"
hpcs.lps.list <- SplitObject(hpcs.lps, split.by = "stim")
hpcs.lps.list <- lapply(hpcs.lps.list, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
})
features <- SelectIntegrationFeatures(hpcs.lps.list)

anchors <- FindIntegrationAnchors(object.list = hpcs.lps.list, anchor.features = features, k.filter = 180)
hpcs.lps <- IntegrateData(anchors)

DefaultAssay(hpcs.lps) <- "integrated"
hpcs.lps <- hpcs.lps %>%
  ScaleData() %>%
  RunPCA(reduction.name = "pca.int", reduction.key = "PCAint_") %>%
  RunUMAP(dims = 1:40, reduction = "pca.int", reduction.name = "umap.int", reduction.key = "UMAPint_")

DefaultAssay(hpcs.lps) <- "RNA"
hpcs.lps <- hpcs.lps %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:40)

d <- DimPlot(hpcs.lps, reduction = "umap")
d_int <- DimPlot(hpcs.lps, reduction = "umap.int")
d | d_int

HoverLocator(d_int, information = FetchData(cells.hpcs, vars = c("cell_type_brief","stim","Gh","Ghrhr","Pomc","Prl","Fshb","Lhb","Tshb")))
cell.id <- "CACCTTAC_YT19070312"
DimPlot(hpcs.lps, cells.highlight = cell.id)

# Idents(hpcs.lps) <- "cell_type_brief"
hpcs.lps <- SetIdent(hpcs.lps, cells = cell.id, value = "Ambig")
d <- DimPlot(hpcs.lps, reduction = "umap")