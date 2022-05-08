library(SCENIC)
library(SCopeLoomR)
library(tidyverse)
library(scales)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_scenic_state_marked.h5Seurat", verbose = F)

# Fig1g
DimPlot(hpcs.lps, group.by = "state", reduction = "umap.scenic", cols = brewer.pal(n = 12, name = "Paired")[c(8,2)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_scenic_hpcs_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/", 
  width = 4, height = 4,
  dpi = 300,
)

# Fig1h
if(!dir.exists('../figures/Fig1/F1h'))dir.create('../figures/Fig1/F1h')
DefaultAssay(hpcs.lps) <- "AUC"
hpcs.lps$scenic_state <- factor(hpcs.lps$scenic_state, levels = c("Healthy","Inflammation"))
## Som
som.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
som.lps <- RunUMAP(som.lps, features = rownames(som.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(som.lps, group.by = "scenic_state", cols = brewer.pal(n = 12, name = "Paired")[c(5,6)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_scenic_som_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1h", 
  width = 3, height = 3,
  dpi = 300,
)

## Lac
lac.lps <- subset(hpcs.lps, subset = cell_type_brief == "Lac")
lac.lps <- RunUMAP(lac.lps, features = rownames(lac.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(lac.lps, group.by = "scenic_state", cols = brewer.pal(n = 12, name = "Paired")[c(7,8)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_scenic_lac_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1h", 
  width = 3, height = 3,
  dpi = 300,
)

## Cort
cort.lps <- subset(hpcs.lps, subset = cell_type_brief == "Cort")
cort.lps <- RunUMAP(cort.lps, features = rownames(cort.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(cort.lps, group.by = "scenic_state", cols = brewer.pal(n = 12, name = "Paired")[c(1,2)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_scenic_cort_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1h", 
  width = 3, height = 2,
  dpi = 300,
)

## Mel
mel.lps <- subset(hpcs.lps, subset = cell_type_brief == "Mel")
mel.lps <- RunUMAP(mel.lps, features = rownames(mel.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(mel.lps, group.by = "scenic_state", cols = brewer.pal(n = 12, name = "Paired")[c(3,4)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_scenic_mel_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1h", 
  width = 2, height = 3,
  dpi = 300,
)

## Gonad
gonad.lps <- subset(hpcs.lps, subset = cell_type_brief == "Gonad")
gonad.lps <- RunUMAP(gonad.lps, features = rownames(gonad.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(gonad.lps, group.by = "scenic_state", cols = c("#7570B3","#BEBADA")) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_scenic_gonad_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1h", 
  width = 3, height = 3,
  dpi = 300,
)

## Thyro
thyro.lps <- subset(hpcs.lps, subset = cell_type_brief == "Thyro")
thyro.lps <- RunUMAP(thyro.lps, features = rownames(thyro.lps), reduction = "pca.scenic", umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(thyro.lps, group.by = "scenic_state", cols = c("#F781BF","#E7298A")) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_scenic_thyro_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1h", 
  width = 3, height = 3,
  dpi = 300,
)




