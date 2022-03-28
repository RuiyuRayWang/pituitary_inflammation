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

# Visualize Louvain clustering and cellassign results side by side
DimPlot(cells.integrated, reduction = "umap.int", group.by = "integrated_snn_res.0.8", label = T) /
  (DimPlot(cells.integrated, group.by = "cell_type_cellassign", reduction = "umap.int", label = T) + NoLegend())

# For some clear cut cell identities, directly assign labels to cells
# Generally, Louvain clustering results are more reliable, cellassign results serves as auxiliary evidence that guides annotation
# For mixed clusters, select cells manually
Idents(cells.integrated) <- "integrated_snn_res.0.8"
## RBCs
DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 14), reduction = "umap.int") /
  DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Red blood cells"), reduction = "umap.int")
rbcs_selected <- CellSelector(DimPlot(cells.integrated, reduction = "umap.int", cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 14)))
rbcs_selected <- union(rbcs_selected, WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 14))
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = rbcs_selected, 
  value = "Red Blood Cells"
  )
## Pou1f1 Progenitors
DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 12), reduction = "umap.int") /
  DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Pou1f1 Progenitors"), reduction = "umap.int")
pou1f1_selected <- CellSelector(DimPlot(cells.integrated, reduction = "umap.int", cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 12)))
pou1f1_selected <- union(pou1f1_selected, WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 12))
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = pou1f1_selected, 
  value = "Pou1f1 Progenitors"
)
## Melanotropes
DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 7), reduction = "umap.int") /
  DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Melanotropes"), reduction = "umap.int")
melanotropes_selected <- CellSelector(DimPlot(cells.integrated, reduction = "umap.int", group.by = "integrated_snn_res.0.8", label = T))
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = melanotropes_selected, 
  value = "Melanotropes"
)
## Corticotropes
DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 1), reduction = "umap.int") /
  DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Corticotropes"), reduction = "umap.int")
corticotropes_selected <- CellSelector(DimPlot(cells.integrated, reduction = "umap.int", group.by = "integrated_snn_res.0.8", label = T))
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = corticotropes_selected, 
  value = "Corticotropes"
)
## Gonadotropes
DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 5), reduction = "umap.int") /
  DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Gonadotropes"), reduction = "umap.int")
gonadotropes_selected <- CellSelector(DimPlot(cells.integrated, reduction = "umap.int", cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 5)))
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = gonadotropes_selected, 
  value = "Gonadotropes"
)
## Thyrotropes
DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = integrated_snn_res.0.8 == 11), reduction = "umap.int") /
  DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Thyrotropes"), reduction = "umap.int")
thyrotropes_selected <- CellSelector(DimPlot(cells.integrated, reduction = "umap.int", cells.highlight = gonadotropes_selected))
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = thyrotropes_selected, 
  value = "Thyrotropes"
)

cells.integrated$cell_type_anno <- Idents(cells.integrated)

# Iterative clustering: Subsetting un-annotated populations
`%ni%` <- Negate(`%in%`)
cells.sub <- subset(cells.integrated, subset = cell_type_anno %ni% c("Corticotropes", "Melanotropes", "Gonadotropes",
                                                                     "Red Blood Cells", "Thyrotropes", "Pou1f1 Progenitors"))
# Rerun Seurat Standard Workflow on subsetted population
DefaultAssay(cells.sub) <- "integrated"
cells.sub <- ScaleData(cells.sub)
cells.sub <- RunPCA(cells.sub, npcs = 50)
cells.sub <- RunUMAP(cells.sub, dims = 1:35)
cells.sub <- FindNeighbors(cells.sub, reduction = "pca", dims = 1:35)
cells.sub <- FindClusters(cells.sub, resolution = 0.8)

DimPlot(cells.sub, reduction = "umap", group.by = "integrated_snn_res.0.8", label = T) / 
  (DimPlot(cells.sub, group.by = "cell_type_cellassign", reduction = "umap", label = T) + NoLegend())

## Pericytes
DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 11)) /
  DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Pericytes"))
pericytes_selected <- CellSelector(DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 11)))
cells.sub <- SetIdent(
  object = cells.sub, 
  cells = pericytes_selected, 
  value = "Pericytes"
)
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = pericytes_selected, 
  value = "Pericytes"
)
## Pituicytes
DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 3)) /
  DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Pituicyte"))
DimPlot(cells.integrated, cells.highlight = WhichCells(cells.integrated, expression = cell_type_cellassign == "Pituicyte"), reduction = "umap") /
  FeaturePlot(cells.integrated, features = "Scn7a", reduction = "umap")
pituicytes_selected <- CellSelector(FeaturePlot(cells.integrated, features = "Scn7a", reduction = "umap"))
cells.sub <- SetIdent(
  object = cells.sub, 
  cells = pituicytes_selected, 
  value = "Pituicytes"
)
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = pituicytes_selected, 
  value = "Pituicytes"
)
## Endo
DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 7)) /
  DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Endothelial cells"))
endo_selected <- CellSelector(DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 7)))
cells.sub <- SetIdent(
  object = cells.sub, 
  cells = endo_selected, 
  value = "Endothelial Cells"
)
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = endo_selected, 
  value = "Endothelial Cells"
)
## WBCs
DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 6)) /
  DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "White blood cells"))
wbcs_selected <- CellSelector(DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 6)))
cells.sub <- SetIdent(
  object = cells.sub, 
  cells = wbcs_selected, 
  value = "White Blood Cells"
)
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = wbcs_selected, 
  value = "White Blood Cells"
)
## Stem cells
DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 8)) /
  DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Stem cells"))
stem_selected <- CellSelector(DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 8)))
cells.sub <- SetIdent(
  object = cells.sub, 
  cells = stem_selected, 
  value = "Stem Cells"
)
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = stem_selected, 
  value = "Stem Cells"
)
## Lactotropes
DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 2)) /
  DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Lactotropes"))
lactotropes_selected <- CellSelector(DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 == 2)))
cells.sub <- SetIdent(
  object = cells.sub, 
  cells = lactotropes_selected, 
  value = "Lactotropes"
)
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = lactotropes_selected, 
  value = "Lactotropes"
)
## Somatotropes
DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 %in% c(0,1,4,5,8,9,10))) /
  DimPlot(cells.sub, reduction = "umap", cells.highlight = WhichCells(cells.sub, expression = cell_type_cellassign == "Lactotropes"))
somatotropes_selected <- CellSelector(DimPlot(cells.sub, cells.highlight = WhichCells(cells.sub, expression = integrated_snn_res.0.8 %in% c(0,1,4,5,8,9,10))))
cells.sub <- SetIdent(
  object = cells.sub, 
  cells = somatotropes_selected, 
  value = "Somatotropes"
)
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = somatotropes_selected, 
  value = "Somatotropes"
)

ambiguous <- WhichCells(cells.integrated, idents = c(0,1,3,4,6,10,13))
cells.integrated <- SetIdent(
  object = cells.integrated, 
  cells = ambiguous, 
  value = "Ambiguous"
)

cells.integrated$cell_type_refined <- Idents(cells.integrated)
cells.integrated$cell_type_refined <- factor(cells.integrated$cell_type_refined, 
                                             levels = c("Somatotropes", "Corticotropes", "Lactotropes", "Melanotropes",
                                                        "Gonadotropes", "Thyrotropes", "Pou1f1 Progenitors", "Stem Cells",
                                                        "White Blood Cells", "Red Blood Cells", "Endothelial Cells",
                                                        "Pericytes", "Pituicytes", "Ambiguous"))
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



save(somatotropes_selected, corticotropes_selected, lactotropes_selected, melanotropes_selected, gonadotropes_selected,
     thyrotropes_selected, pou1f1_selected, stem_selected, wbcs_selected, rbcs_selected, endo_selected, pericytes_selected,
     pituicytes_selected, ambiguous, 
     file = paste0("data/selected_cells_",Sys.Date(),".rda"))
SaveH5Seurat(object = cells.integrated, filename = "data/cells_postprocessed.h5Seurat", overwrite = T)
