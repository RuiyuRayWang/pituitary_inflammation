library(SCENIC)
library(SCopeLoomR)
library(tidyverse)
library(Seurat)
library(SeuratDisk)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat", verbose = F)
metadata <- hpcs.lps@meta.data
write.csv(metadata, "../data/processed/hpcs_lps_metadata.csv")

pyScenicDir <- '../data/scenic_protocol'
pyScenicLoomFile <- file.path(pyScenicDir, "files", "hpcs_lps_pyscenic_output.loom")
loom <- open_loom(pyScenicLoomFile, mode="r")

## Read information from loom file: Regulon AUC
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

close_loom(loom)

## Manually explore and select binarization thresholds for AUCell scores, using shiny app.
thresholds_csv <- read.csv("../data/scenic_protocol/files/thresholds.csv", row.names = 1)
thresholds <- thresholds_csv$threshold; names(thresholds) <- rownames(thresholds_csv)
DefaultAssay(hpcs.lps) <- "RNA"
tsne_coords <- hpcs.lps %>% 
  RunTSNE(reduction = "pca", dims = 1:50) %>% 
  Embeddings(reduction = "tsne") %>% 
  as.data.frame() %>% 
  rename("_X" = "tSNE_1", "_Y" = "tSNE_2") %>%
  as.matrix()
aucellApp <- AUCell_createViewerApp(
  auc = regulonsAUC,
  thresholds = thresholds,
  tSNE = tsne_coords
  )
savedSelections <- shiny::runApp(aucellApp)
thresholds <- savedSelections$thresholds
saveRDS(savedSelections, '../data/scenic_protocol/files/hpcslps_savedSelections.rds')

## Update and export
thresholds_csv <- as.data.frame(thresholds, row.names = names(thresholds)); colnames(thresholds_csv) <- "threshold"
write.csv(thresholds_csv, "../data/scenic_protocol/files/thresholds.csv")

auc_mtx <- getAUC(regulonsAUC)
auc_mtx <- auc_mtx[names(thresholds),]

## Binarize SCENIC AUC using manually curated thresholds
bin_mtx <- apply(X = auc_mtx, MARGIN = 2, FUN = function(x){
  ifelse(x > thresholds, 1, 0)
})

write.csv(x = bin_mtx, file = "../data/scenic_protocol/files/bin_mtx.csv")

hpcs.lps[["AUC"]] <- CreateAssayObject(
  data = getAUC(regulonsAUC)
)

hpcs.lps[["BIN"]] <- CreateAssayObject(
  count = bin_mtx
)

## Perform dimension reduction on AUC matrix of all HPCs
DefaultAssay(hpcs.lps) <- "AUC"
hpcs.lps <- RunUMAP(
  hpcs.lps, 
  umap.method = "umap-learn",
  features = rownames(hpcs.lps), 
  reduction.name = "umap.scenic", 
  reduction.key = "UMAPscenic_"
)


## Som
som.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
DefaultAssay(som.lps) <- "AUC"
som.lps <- FindVariableFeatures(som.lps, assay = "AUC", selection.method = "disp", nfeatures = 200)
# som.lps@assays$AUC@scale.data <- GetAssayData(som.lps, slot = "data")[VariableFeatures(som.lps),]
som.lps <- ScaleData(som.lps, do.scale = FALSE, do.center = TRUE)
som.lps <- RunPCA(som.lps, reduction.name = "pca.scenic", reduction.key = "PCAscenic_")
som.lps <- FindNeighbors(som.lps, reduction = "pca.scenic")
res = 0.1
som.lps <- FindClusters(som.lps, graph.name = "AUC_snn", resolution = res)
som.lps <- RunUMAP(som.lps, features = rownames(som.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
# DimPlot(som.lps, reduction = "umap.scenic", group.by = paste0("AUC_snn_res.",res)) | DimPlot(som.lps, reduction = "umap.scenic", group.by = "state")
som.lps <- RenameIdents(som.lps, `0` = "Inflammation", `1` = "Healthy")
som.lps$scenic_state <- factor(Idents(som.lps), levels = c("Healthy","Inflammation"))
for (s in levels(som.lps$scenic_state)){
  hpcs.lps <- SetIdent(object = hpcs.lps, cells = WhichCells(som.lps, expression = scenic_state == s), value = s)
}

## Lac
lac.lps <- subset(hpcs.lps, subset = cell_type_brief == "Lac")
DefaultAssay(lac.lps) <- "AUC"
lac.lps <- FindVariableFeatures(lac.lps, assay = "AUC", selection.method = "disp", nfeatures = 200)
# lac.lps@assays$AUC@scale.data <- GetAssayData(lac.lps, slot = "data")[VariableFeatures(lac.lps),]
lac.lps <- ScaleData(lac.lps, do.scale = FALSE, do.center = TRUE)
lac.lps <- RunPCA(lac.lps, reduction.name = "pca.scenic", reduction.key = "PCAscenic_")
lac.lps <- FindNeighbors(lac.lps, reduction = "pca.scenic")
res = 0.3
lac.lps <- FindClusters(lac.lps, graph.name = "AUC_snn", resolution = res)
lac.lps <- RunUMAP(lac.lps, features = rownames(lac.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(lac.lps, reduction = "umap.scenic", group.by = paste0("AUC_snn_res.",res)) | DimPlot(lac.lps, reduction = "umap.scenic", group.by = "state")
lac.lps <- RenameIdents(lac.lps, `0` = "Healthy", `1` = "Inflammation")
lac.lps$scenic_state <- factor(Idents(lac.lps), levels = c("Healthy","Inflammation"))
for (s in levels(lac.lps$scenic_state)){
  hpcs.lps <- SetIdent(object = hpcs.lps, cells = WhichCells(lac.lps, expression = scenic_state == s), value = s)
}

## Cort
cort.lps <- subset(hpcs.lps, subset = cell_type_brief == "Cort")
DefaultAssay(cort.lps) <- "AUC"
cort.lps <- FindVariableFeatures(cort.lps, assay = "AUC", selection.method = "disp", nfeatures = 200)
# cort.lps@assays$AUC@scale.data <- GetAssayData(cort.lps, slot = "data")[VariableFeatures(cort.lps),]
cort.lps <- ScaleData(cort.lps, do.scale = FALSE, do.center = TRUE)
cort.lps <- RunPCA(cort.lps, reduction.name = "pca.scenic", reduction.key = "PCAscenic_")
cort.lps <- FindNeighbors(cort.lps, reduction = "pca.scenic")
res = 0.3
cort.lps <- FindClusters(cort.lps, graph.name = "AUC_snn", resolution = res)
cort.lps <- RunUMAP(cort.lps, features = rownames(cort.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(cort.lps, reduction = "umap.scenic", group.by = paste0("AUC_snn_res.",res)) | DimPlot(cort.lps, reduction = "umap.scenic", group.by = "state")
cort.lps <- RenameIdents(cort.lps, `0` = "Healthy", `1` = "Inflammation")
cort.lps$scenic_state <- factor(Idents(cort.lps), levels = c("Healthy","Inflammation"))
for (s in levels(cort.lps$scenic_state)){
  hpcs.lps <- SetIdent(object = hpcs.lps, cells = WhichCells(cort.lps, expression = scenic_state == s), value = s)
}

## Mel
mel.lps <- subset(hpcs.lps, subset = cell_type_brief == "Mel")
DefaultAssay(mel.lps) <- "AUC"
mel.lps <- FindVariableFeatures(mel.lps, assay = "AUC", selection.method = "disp", nfeatures = 250)
# mel.lps@assays$AUC@scale.data <- GetAssayData(mel.lps, slot = "data")[VariableFeatures(mel.lps),]
mel.lps <- ScaleData(mel.lps, do.scale = FALSE, do.center = TRUE)
mel.lps <- RunPCA(mel.lps, reduction.name = "pca.scenic", reduction.key = "PCAscenic_")
mel.lps <- FindNeighbors(mel.lps, reduction = "pca.scenic")
res = 0.3
mel.lps <- FindClusters(mel.lps, graph.name = "AUC_snn", resolution = res)
mel.lps <- RunUMAP(mel.lps, features = rownames(mel.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(mel.lps, reduction = "umap.scenic", group.by = paste0("AUC_snn_res.",res)) | DimPlot(mel.lps, reduction = "umap.scenic", group.by = "state")
mel.lps <- RenameIdents(mel.lps, `0` = "Healthy", `1` = "Inflammation")
mel.lps$scenic_state <- factor(Idents(mel.lps), levels = c("Healthy","Inflammation"))
for (s in levels(mel.lps$scenic_state)){
  hpcs.lps <- SetIdent(object = hpcs.lps, cells = WhichCells(mel.lps, expression = scenic_state == s), value = s)
}

## Gonad
gonad.lps <- subset(hpcs.lps, subset = cell_type_brief == "Gonad")
DefaultAssay(gonad.lps) <- "AUC"
gonad.lps <- FindVariableFeatures(gonad.lps, assay = "AUC", selection.method = "disp", nfeatures = 250)
# gonad.lps@assays$AUC@scale.data <- GetAssayData(gonad.lps, slot = "data")[VariableFeatures(gonad.lps),]
gonad.lps <- ScaleData(gonad.lps, do.scale = FALSE, do.center = TRUE)
gonad.lps <- RunPCA(gonad.lps, reduction.name = "pca.scenic", reduction.key = "PCAscenic_")
gonad.lps <- FindNeighbors(gonad.lps, reduction = "pca.scenic")
res = 0.5
gonad.lps <- FindClusters(gonad.lps, graph.name = "AUC_snn", resolution = res)
gonad.lps <- RunUMAP(gonad.lps, features = rownames(gonad.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(gonad.lps, reduction = "umap.scenic", group.by = paste0("AUC_snn_res.",res)) | DimPlot(gonad.lps, reduction = "umap.scenic", group.by = "state")
gonad.lps <- RenameIdents(gonad.lps, `0` = "Healthy", `1` = "Inflammation")
gonad.lps$scenic_state <- factor(Idents(gonad.lps), levels = c("Healthy","Inflammation"))
for (s in levels(gonad.lps$scenic_state)){
  hpcs.lps <- SetIdent(object = hpcs.lps, cells = WhichCells(gonad.lps, expression = scenic_state == s), value = s)
}

## Thyro
thyro.lps <- subset(hpcs.lps, subset = cell_type_brief == "Thyro")
DefaultAssay(thyro.lps) <- "AUC"
thyro.lps <- FindVariableFeatures(thyro.lps, assay = "AUC", selection.method = "disp", nfeatures = 300)
# thyro.lps@assays$AUC@scale.data <- GetAssayData(thyro.lps, slot = "data")[VariableFeatures(thyro.lps),]
thyro.lps <- ScaleData(thyro.lps, do.scale = FALSE, do.center = TRUE)
thyro.lps <- RunPCA(thyro.lps, reduction.name = "pca.scenic", reduction.key = "PCAscenic_")
thyro.lps <- FindNeighbors(thyro.lps, reduction = "pca.scenic")
res = 0.8
thyro.lps <- FindClusters(thyro.lps, graph.name = "AUC_snn", resolution = res)
thyro.lps <- RunUMAP(thyro.lps, features = rownames(thyro.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(thyro.lps, reduction = "umap.scenic", group.by = paste0("AUC_snn_res.",res)) | DimPlot(thyro.lps, reduction = "umap.scenic", group.by = "state")
thyro.lps <- RenameIdents(thyro.lps, `0` = "Healthy", `1` = "Inflammation")
thyro.lps$scenic_state <- factor(Idents(thyro.lps), levels = c("Healthy","Inflammation"))
for (s in levels(thyro.lps$scenic_state)){
  hpcs.lps <- SetIdent(object = hpcs.lps, cells = WhichCells(thyro.lps, expression = scenic_state == s), value = s)
}

hpcs.lps$scenic_state <- Idents(hpcs.lps)
hpcs.lps$scenic_state <- as.character(hpcs.lps$scenic_state)

SaveH5Seurat(object = hpcs.lps, filename = "../data/processed/hpcs_lps_scenic_state_marked.h5Seurat", overwrite = T, verbose = F)
