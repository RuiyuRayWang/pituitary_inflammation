library(Seurat)
library(SeuratDisk)
library(tidyverse)

cells <- LoadH5Seurat(file = "data/cells_postprocessed.h5Seurat")

hpcs <- subset(cells, subset = cell_type_brief %in% c("Som","Cort","Lac","Mel","Gonad","Thyro","Pitui"))

hpcs.lps <- subset(hpcs, subset = treat %in% c("Saline","LPS"))


SSW <- function(object, def_assay = "RNA", n_feat = 2000, npcs = 50, dims_use = 1:50, n.neighbors = 30,
                res = 0.8){
  DefaultAssay(object = object) <- def_assay
  object <- NormalizeData(object = object)
  object <- FindVariableFeatures(object = object, nfeatures = n_feat)
  object <- ScaleData(object = object)
  object <- RunPCA(object = object, npcs = npcs)
  object <- FindNeighbors(object = object, dims = dims_use)
  object <- FindClusters(object = object, resolution = res)
  object <- RunUMAP(object = object, n.neighbors = n.neighbors, dims = dims_use)
}

hpcs.lps <- SSW(hpcs.lps)

som.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
cort.lps <- subset(hpcs.lps, subset = cell_type_brief == "Cort")
lac.lps <- subset(hpcs.lps, subset = cell_type_brief == "Lac")
mel.lps <- subset(hpcs.lps, subset = cell_type_brief == "Mel")
gonad.lps <- subset(hpcs.lps, subset = cell_type_brief == "Gonad")
thyro.lps <- subset(hpcs.lps, subset = cell_type_brief == "Thyro")
pitui.lps <- subset(hpcs.lps, subset = cell_type_brief == "Pitui")

som.lps <- SSW(som.lps, n_feat = 1000, dims_use = 1:15)
cort.lps <- SSW(cort.lps, n_feat = 1500, dims_use = 1:10)
lac.lps <- SSW(lac.lps, dims_use = 1:10)
mel.lps <- SSW(mel.lps, dims_use = 1:10, res = 0.1)
gonad.lps <- SSW(gonad.lps, dims_use = 1:10)
thyro.lps <- SSW(thyro.lps, dims_use = 1:8)
pitui.lps <- SSW(pitui.lps, npcs = 10, dims_use = 1:10,  res = 1)

Idents(mel.lps) <- "RNA_snn_res.0.1"
mel.lps <- RenameIdents(mel.lps, `0` = "Inflammation", `1` = "Healthy")
mel.mks.infl_vs_ctrl <- FindMarkers(mel.lps, ident.1 = "Inflammation", only.pos = T)


Idents(pitui.lps) <- "RNA_snn_res.1"
# pitui.lps <- RenameIdents(pitui.lps, `0` = "Healthy", `1` = "Inflammation")
pitui.lps <- CellSelector(plot = DimPlot(pitui.lps, group.by = "treat"), object = pitui.lps, ident = "Healthy")
pitui.lps <- CellSelector(plot = DimPlot(pitui.lps, group.by = "treat"), object = pitui.lps, ident = "Inflammation")
pitui.mks.infl_vs_ctrl <- FindMarkers(pitui.lps, ident.1 = "Inflammation", only.pos = T)

write.csv(mel.mks.infl_vs_ctrl, file = "outs/melanotropes.mks.infl_vs_ctrl.csv", quote = F)
write.csv(mel.mks.infl_vs_ctrl, file = "outs/pituicytes.mks.infl_vs_ctrl.csv", quote = F)