#!/usr/bin/env Rscript

## Post-processing workflow
### Seurat integration
### Integrated analysis
### Cell type annotation and refinement
### Standard workflow without integration

library(Seurat)
library(SeuratDisk)
library(tidyverse)

cells <- LoadH5Seurat('data/cells_assigned.h5Seurat')

cells.list <- SplitObject(cells, split.by = "stim")

# Normalize, feature selection by each stim
for (i in 1:length(cells.list)){
  cells.list[[i]] <- NormalizeData(cells.list[[i]], verbose = TRUE)
  cells.list[[i]] <- FindVariableFeatures(cells.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}
features <- SelectIntegrationFeatures(object.list = cells.list)

# Perform integration
cells.anchor <- FindIntegrationAnchors(object.list = cells.list, anchor.features = features, dims = 1:50, k.filter = 150)
cells.integrated <- IntegrateData(anchorset = cells.anchor)

# Perform an integrated analysis
DefaultAssay(cells.integrated) <- "integrated"

cells.integrated <- ScaleData(cells.integrated)
cells.integrated <- RunPCA(cells.integrated, npcs = 50, reduction.name = "pca.int", reduction.key = "PCAint_")
cells.integrated <- RunUMAP(cells.integrated, dims = 1:50, reduction = "pca.int", reduction.name = "umap.int", reduction.key = "UMAPint_")
cells.integrated <- FindNeighbors(cells.integrated, reduction = "pca.int", dims = 1:50)
cells.integrated <- FindClusters(cells.integrated, resolution = 0.8)  


# Also proceed through un-integrated analysis
DefaultAssay(cells.integrated) <- "RNA"
cells.integrated <- NormalizeData(cells.integrated)
cells.integrated <- FindVariableFeatures(cells.integrated)
cells.integrated <- ScaleData(cells.integrated)
cells.integrated <- RunPCA(cells.integrated, npcs = 50)

# cells.integrated <- JackStraw(cells.integrated, dims = 50)
# cells.integrated <- ScoreJackStraw(cells.integrated, dims = 1:50)
# dims_use <- JS(cells.integrated[["pca"]], slot = "overall") %>% as.data.frame %>% filter(Score < 0.05) %>% pull(PC)
dims_use <- 1:50

cells.integrated <- FindNeighbors(cells.integrated, dims = dims_use)
cells.integrated <- FindClusters(cells.integrated)
cells.integrated <- RunUMAP(cells.integrated, dims = dims_use)

SaveH5Seurat(cells.integrated, filename = 'data/cells_integrated.h5Seurat', overwrite = T, verbose = F)

# Visualize Louvain clustering and cellassign results side by side
d1 <- DimPlot(cells.integrated, reduction = "umap.int", group.by = "integrated_snn_res.0.8", label = T) + NoLegend()
d2 <- DimPlot(cells.integrated, group.by = "cell_type_cellassign", reduction = "umap.int", label = T) + NoLegend()
d3 <- DimPlot(cells.integrated, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T) + NoLegend()
d4 <- DimPlot(cells.integrated, group.by = "cell_type_cellassign", reduction = "umap", label = T) + NoLegend()
(d1 | d2) / (d3 | d4)

# For some clear cut cell identities, directly assign labels to cells
# Generally, Louvain clustering results are more reliable, cellassign results serves as auxiliary evidence that guides annotation
# For mixed clusters, select cells manually
# Idents(cells.integrated) <- "cell_type_cellassign"
Idents(cells.integrated) <- "integrated_snn_res.0.8"
## RBCs
f1 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 13), reduction = "umap.int") + NoLegend() + ggtitle("14")
f2 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Red Blood Cells"), reduction = "umap.int") + NoLegend() + ggtitle("RBCs")
f3 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 13), reduction = "umap") + NoLegend() + ggtitle("14")
f4 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Red Blood Cells"), reduction = "umap") + NoLegend() + ggtitle("RBCs")
(f1 | f2) / (f3 | f4)
cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Red Blood Cells")

# rbcs_selected <- CellSelector(DimPlot(cells.integrated, reduction = "umap.int", cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 14)))
# rbcs_selected <- union(rbcs_selected, WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 14))
# cells.integrated <- SetIdent(
#   object = cells.integrated, 
#   cells = rbcs_selected, 
#   value = "Red Blood Cells"
#   )

## Pou1f1 Progenitors
f1 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 11), reduction = "umap.int") + NoLegend() + ggtitle("12")
f2 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Pou1f1 Progenitors"), reduction = "umap.int") + NoLegend() + ggtitle("Pou1f1")
f3 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 11), reduction = "umap") + NoLegend() + ggtitle("12")
f4 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Pou1f1 Progenitors"), reduction = "umap") + NoLegend() + ggtitle("Pou1f1")
(f1 | f2) / (f3 | f4)
cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Pou1f1 Progenitors")
# cells.integrated <- CellSelector(plot = f3, object = cells.integrated, ident = "Pou1f1 Progenitors")

## Melanotropes
f1 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 7), reduction = "umap.int") + NoLegend() + ggtitle("7")
f2 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Melanotropes"), reduction = "umap.int") + NoLegend() + ggtitle("Mel")
f3 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 7), reduction = "umap") + NoLegend() + ggtitle("7")
f4 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Melanotropes"), reduction = "umap") + NoLegend() + ggtitle("Mel")
(f1 | f2) / (f3 | f4)
cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Melanotropes")

## Corticotropes
f1 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 2), reduction = "umap.int") + NoLegend() + ggtitle("1")
f2 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Corticotropes"), reduction = "umap.int") + NoLegend() + ggtitle("Cort")
f3 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 2), reduction = "umap") + NoLegend() + ggtitle("1")
f4 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Corticotropes"), reduction = "umap") + NoLegend() + ggtitle("Cort")
(f1 | f2) / (f3 | f4)
cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Corticotropes")

## Gonadotropes
f1 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 5), reduction = "umap.int") + NoLegend() + ggtitle("5")
f2 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Gonadotropes"), reduction = "umap.int") + NoLegend() + ggtitle("Gonad")
f3 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 5), reduction = "umap") + NoLegend() + ggtitle("5")
f4 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Gonadotropes"), reduction = "umap") + NoLegend() + ggtitle("Gonad")
(f1 | f2) / (f3 | f4)
cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Gonadotropes")

## Thyrotropes
f1 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 10), reduction = "umap.int") + NoLegend() + ggtitle("11")
f2 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Thyrotropes"), reduction = "umap.int") + NoLegend() + ggtitle("Thyro")
f3 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 10), reduction = "umap") + NoLegend() + ggtitle("11")
f4 <- DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Thyrotropes"), reduction = "umap") + NoLegend() + ggtitle("Thyro")
(f1 | f2) / (f3 | f4)
cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Thyrotropes")

cells.integrated$cell_type_anno <- Idents(cells.integrated)

# Iterative clustering: Subsetting un-annotated populations
`%ni%` <- Negate(`%in%`)
cells.sub <- subset(cells.integrated, subset = cell_type_anno %ni% c("Corticotropes", "Melanotropes", "Gonadotropes",
                                                                     "Red Blood Cells", "Thyrotropes", "Pou1f1 Progenitors"))
# Rerun Seurat Standard Workflow on subsetted population
DefaultAssay(cells.sub) <- "integrated"
cells.sub <- ScaleData(cells.sub)
cells.sub <- RunPCA(cells.sub, npcs = 50, reduction.name = "pca.int", reduction.key = "PCAint_")
cells.sub <- RunUMAP(cells.sub, dims = 1:35, reduction = "pca.int", reduction.name = "umap.int", reduction.key = "UMAPint_")
cells.sub <- FindNeighbors(cells.sub, reduction = "pca.int", dims = 1:35)
cells.sub <- FindClusters(cells.sub, resolution = 0.8)

DefaultAssay(cells.sub) <- "RNA"
cells.sub <- NormalizeData(cells.sub)
cells.sub <- FindVariableFeatures(cells.sub)
cells.sub <- ScaleData(cells.sub)
cells.sub <- RunPCA(cells.sub, npcs = 50)
cells.sub <- RunUMAP(cells.sub, dims = 1:35)
cells.sub <- FindNeighbors(cells.sub, reduction = "pca", dims = 1:35)
cells.sub <- FindClusters(cells.sub, resolution = 0.8)

d1 <- DimPlot(cells.sub, reduction = "umap.int", group.by = "integrated_snn_res.0.8", label = T) + NoLegend()
d2 <- DimPlot(cells.sub, group.by = "cell_type_cellassign", reduction = "umap.int", label = T) + NoLegend()
d3 <- DimPlot(cells.sub, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T) + NoLegend()
d4 <- DimPlot(cells.sub, group.by = "cell_type_cellassign", reduction = "umap", label = T) + NoLegend()
(d1 | d2) / (d3 | d4)

Idents(cells.sub) <- "integrated_snn_res.0.8"
## Pericytes
f1 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 13), reduction = "umap.int") + NoLegend() + ggtitle("13")
f2 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Pericytes"), reduction = "umap.int") + NoLegend() + ggtitle("Peri")
f3 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 13), reduction = "umap") + NoLegend() + ggtitle("13")
f4 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Pericytes"), reduction = "umap") + NoLegend() + ggtitle("Peri")
(f1 | f2) / (f3 | f4)
cells.sub <- CellSelector(plot = f1, object = cells.sub, ident = "Pericytes")
# cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Pericytes")
cells.integrated <- SetIdent(cells.integrated, 
                             cells = WhichCells(cells.sub, idents = "Pericytes"), 
                             value = "Pericytes")

## Pituicytes
f1 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 11), reduction = "umap.int") + NoLegend() + ggtitle("11")
f2 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Pituicyte"), reduction = "umap.int") + NoLegend() + ggtitle("Pitui")
f3 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 11), reduction = "umap") + NoLegend() + ggtitle("1")
f4 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Pituicyte"), reduction = "umap") + NoLegend() + ggtitle("Pitui")
(f1 | f2) / (f3 | f4)
cells.sub <- CellSelector(plot = f1, object = cells.sub, ident = "Pituicytes")
# cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Pituicytes")
cells.integrated <- SetIdent(cells.integrated, 
                             cells = WhichCells(cells.sub, idents = "Pituicytes"), 
                             value = "Pituicytes")

## Endo
f1 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 9), reduction = "umap.int") + NoLegend() + ggtitle("9")
f2 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Endothelial Cells"), reduction = "umap.int") + NoLegend() + ggtitle("Endo")
f3 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 9), reduction = "umap") + NoLegend() + ggtitle("9")
f4 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Endothelial Cells"), reduction = "umap") + NoLegend() + ggtitle("Endo")
(f1 | f2) / (f3 | f4)
cells.sub <- CellSelector(plot = f1, object = cells.sub, ident = "Endothelial Cells")
# cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Endothelial Cells")
cells.integrated <- SetIdent(cells.integrated, 
                             cells = WhichCells(cells.sub, idents = "Endothelial Cells"), 
                             value = "Endothelial Cells")

## WBCs
f1 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 6), reduction = "umap.int") + NoLegend() + ggtitle("6")
f2 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "White Blood Cells"), reduction = "umap.int") + NoLegend() + ggtitle("WBCs")
f3 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 6), reduction = "umap") + NoLegend() + ggtitle("6")
f4 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "White Blood Cells"), reduction = "umap") + NoLegend() + ggtitle("WBCs")
(f1 | f2) / (f3 | f4)
cells.sub <- CellSelector(plot = f1, object = cells.sub, ident = "White Blood Cells")
# cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "White Blood Cells")
cells.integrated <- SetIdent(cells.integrated, 
                             cells = WhichCells(cells.sub, idents = "White Blood Cells"), 
                             value = "White Blood Cells")

## Stem cells
f1 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 10), reduction = "umap.int") + NoLegend() + ggtitle("10")
f2 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Stem Cells"), reduction = "umap.int") + NoLegend() + ggtitle("Stem")
f3 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 10), reduction = "umap") + NoLegend() + ggtitle("10")
f4 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Stem Cells"), reduction = "umap") + NoLegend() + ggtitle("Stem")
(f1 | f2) / (f3 | f4)
cells.sub <- CellSelector(plot = f1, object = cells.sub, ident = "Stem Cells")
# cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Stem Cells")
cells.integrated <- SetIdent(cells.integrated, 
                             cells = WhichCells(cells.sub, idents = "Stem Cells"), 
                             value = "Stem Cells")

## Lactotropes
f1 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 2), reduction = "umap.int") + NoLegend() + ggtitle("2")
f2 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Lactotropes"), reduction = "umap.int") + NoLegend() + ggtitle("Lac")
f3 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 2), reduction = "umap") + NoLegend() + ggtitle("2")
f4 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Lactotropes"), reduction = "umap") + NoLegend() + ggtitle("Lac")
(f1 | f2) / (f3 | f4)
cells.sub <- CellSelector(plot = f1, object = cells.sub, ident = "Lactotropes")
# cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Lactotropes")
cells.integrated <- SetIdent(cells.integrated, 
                             cells = WhichCells(cells.sub, idents = "Lactotropes"), 
                             value = "Lactotropes")

## Somatotropes
f1 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 %in% c(0,1,4,5,7,8,12)), reduction = "umap.int") + NoLegend() + ggtitle("")
f2 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Lactotropes"), reduction = "umap.int") + NoLegend() + ggtitle("Som")
f3 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 %in% c(0,1,4,5,7,8,12)), reduction = "umap") + NoLegend() + ggtitle("")
f4 <- DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Lactotropes"), reduction = "umap") + NoLegend() + ggtitle("Som")
(f1 | f2) / (f3 | f4)
cells.sub <- CellSelector(plot = f1, object = cells.sub, ident = "Somatotropes")
# cells.integrated <- CellSelector(plot = f1, object = cells.integrated, ident = "Somatotropes")
cells.integrated <- SetIdent(cells.integrated, 
                             cells = WhichCells(cells.sub, idents = "Somatotropes"), 
                             value = "Somatotropes")

cells.integrated$cell_type_anno <- Idents(cells.integrated)
ambiguous <- WhichCells(cells.integrated, expression = cell_type_anno %ni% c("Somatotropes", "Lactotropes", "Corticotropes", "Melanotropes", "Gonadotropes", 
                                                                             "Thyrotropes", "Stem Cells", "Endothelial Cells", "Red Blood Cells", "White Blood Cells", 
                                                                             "Pou1f1 Progenitors", "Pituicytes", "Pericytes"))
cells.integrated <- SetIdent(object = cells.integrated, cells = ambiguous, value = "Ambiguous")

cells.integrated$cell_type_refined <- Idents(cells.integrated)


cells.integrated$cell_type_refined <- factor(cells.integrated$cell_type_refined, 
                                             levels = c("Somatotropes", "Corticotropes", "Lactotropes", "Melanotropes",
                                                        "Gonadotropes", "Thyrotropes", "Pou1f1 Progenitors", "Stem Cells",
                                                        "White Blood Cells", "Red Blood Cells", "Endothelial Cells",
                                                        "Pericytes", "Pituicytes", "Ambiguous"))
Idents(cells.integrated) <- "cell_type_refined"
cells.integrated <- RenameIdents(cells.integrated, `Somatotropes` = "Som", `Corticotropes` = "Cort", `Lactotropes` = "Lac", 
                                 `Melanotropes` = "Mel", `Gonadotropes` = "Gonad", `Thyrotropes` = "Thyro",
                                 `Pou1f1 Progenitors` = "Pou1f1", `Stem Cells` = "Stem", `White Blood Cells` = "WBCs",
                                 `Red Blood Cells` = "RBCs", `Endothelial Cells` = "Endo", `Pericytes` = "Peri",
                                 `Pituicytes` = "Pitui", `Ambiguous` = "Ambig")
cells.integrated$cell_type_brief <- Idents(cells.integrated)
cells.integrated$cell_type_brief <- factor(cells.integrated$cell_type_brief,
                                           levels = c("Som", "Cort", "Lac", "Mel",
                                                      "Gonad", "Thyro", "Pou1f1", "Stem",
                                                      "WBCs", "RBCs", "Endo",
                                                      "Peri", "Pitui", "Ambig"))
cells.integrated$cell_type_anno <- NULL

cells.integrated$cell_type_refined <- as.character(cells.integrated$cell_type_refined)
cells.integrated$cell_type_brief <- as.character(cells.integrated$cell_type_brief)

SaveH5Seurat(object = cells.integrated, filename = "data/cells_postprocessed.h5Seurat", overwrite = T)
Convert("data/cells_postprocessed.h5Seurat", dest = "h5ad", overwrite = T)
