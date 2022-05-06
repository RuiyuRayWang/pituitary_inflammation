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

hpcs.lps <- LoadH5Seurat("../data/hpcs_lps_state_marked.h5Seurat", verbose = F)
metadata <- hpcs.lps@meta.data

pyScenicDir <- '../scenic_protocol'
pyScenicLoomFile <- file.path(pyScenicDir, "files", "hpcs_lps_pyscenic_output.loom")
loom <- open_loom(pyScenicLoomFile, mode="r")

## Read information from loom file: Regulon AUC
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

close_loom(loom)

bin_mtx <- t(read.csv("../scenic_protocol/files/bin_mtx.csv", row.names = 1))

hpcs.lps[["AUC"]] <- CreateAssayObject(
  data = getAUC(regulonsAUC)
)

DefaultAssay(hpcs.lps) <- "AUC"
hpcs.lps <- RunUMAP(
  hpcs.lps, 
  umap.method = "umap-learn",
  features = rownames(hpcs.lps), 
  reduction.name = "umap.scenic", 
  reduction.key = "UMAPscenic_",
  # a = 1.2, b = 1.5,
  # n.epochs = 1000
  )

# Fig1g
DimPlot(hpcs.lps, group.by = "state", reduction = "umap.scenic", cols = brewer.pal(n = 12, name = "Paired")[c(8,2)], pt.size = .6) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
# FetchData(
#   object = hpcs.lps, 
#   vars = c("UMAPscenic_1", "UMAPscenic_2", "state", "cell_type_brief")
#   ) %>%
#   mutate(state = factor(state, levels = c("Healthy","Inflammation")),
#          cell_type_brief = factor(cell_type_brief, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))) %>%
#   ggplot(aes(x = UMAPscenic_1, y = UMAPscenic_2)) +
#   geom_point(aes(color = state, fill = cell_type_brief), shape = 21, size = 1.2, stroke = .8) +
#   # geom_mark_hull(aes(color = cell_type_brief), concavity = 5) +
#   scale_fill_manual(values = hue_pal()(13)[c(1,2,9,6,10,11)]) +
#   scale_color_manual(values = brewer.pal(n = 8, name = "Paired")[c(2,8)]) +
#   # scale_shape_manual(values = c(20,5)) +
#   theme(
#     axis.line = element_blank(),
#     axis.title = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     panel.background = element_blank(),
#     legend.position = "none",
#   )
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
## Som
som.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
som.lps <- RunUMAP(som.lps, features = rownames(som.lps), umap.method = "umap-learn", reduction.name = "umap.scenic", reduction.key = "UMAPscenic_")
DimPlot(som.lps, group.by = "cell_type_brief", cols = "#F8766D") +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_som_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1h", 
  width = 3, height = 4,
  dpi = 300,
)
DefaultAssay(som.lps) <- "AUC"

auc_mtx <- t(GetAssayData(som.lps, slot = "data", assay = "AUC"))

som_auc_mtx <- auc_mtx[Cells(som.lps),]
cm <- e1071::cmeans(som_auc_mtx, centers = 2)
method_use <- c()
d <- dist(som_auc_mtx, method = "euclidean")
hclust()
km <- kmeans(x = som_auc_mtx, centers = 2, iter.max = 20)
DimPlot(som.lps, reduction = "umap.scenic", cells.highlight = names(km$cluster[km$cluster==1]))

res_dbscan <- dbscan::dbscan(som_auc_mtx, eps = 5, minPts = 4)

DefaultAssay(som.lps) <- "AUC"
som.lps <- FindVariableFeatures(som.lps, assay = "AUC", selection.method = "disp", nfeatures = 200)
# som.lps@assays$AUC@scale.data <- GetAssayData(som.lps, slot = "data")[VariableFeatures(som.lps),]
som.lps <- ScaleData(som.lps, do.scale = FALSE, do.center = TRUE)
som.lps <- RunPCA(som.lps, reduction.name = "pca.scenic", reduction.key = "PCAscenic_")
som.lps <- FindNeighbors(som.lps, reduction = "pca.scenic")
som.lps <- FindClusters(som.lps, graph.name = "AUC_snn", resolution = 0.1)
DimPlot(som.lps, reduction = "umap.scenic", group.by = "AUC_snn_res.0.1") | DimPlot(som.lps, reduction = "umap.scenic", group.by = "state")

clust_res <- NbClust::NbClust(som_auc_mtx, distance = "euclidean", method = "ward.D", min.nc = 2, max.nc = 4)

ht <- Heatmap(
  matrix = t(auc_mtx)[,Cells(som.lps)], 
  cluster_rows = T, 
  cluster_columns = T, 
  show_column_names = F, 
  show_row_names = F, 
  column_km = 2
  )
ht <- draw(ht)
c.dend <- column_dend(ht)
cclust.list <- column_order(ht)

pheatmap::pheatmap(
  som_auc_mtx, 
  
  show_rownames = F, 
  show_colnames = F)

# res <- mixtools::mvnormalmixEM(
#   x = som_auc_mtx,
#   epsilon = 1e-02,
#   verb = TRUE
# )
km <- kmeans(
  x = som_auc_mtx,
  centers = 2, nstart = 1, algorithm = "MacQueen"
  )

## Lac
lac.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
lac.lps %>%
  RunUMAP(features = rownames(hpcs.lps), umap.method = "umap-learn") %>%
  DimPlot(group.by = "cell_type_brief", cols = "#E18A00") +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_lac_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1h", 
  width = 3, height = 3,
  dpi = 300,
)