library(Seurat)
library(ggplot2)

cells <- readRDS('data/cells.rds')
fit <- readRDS('data/fit_cellassign_20211205.rds')

cells[["cellassign_celltype"]] <- fit$cell



# p1 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Somatotropes")) + NoLegend()
p2 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Lactotropes")) + NoLegend() + ggtitle("Lactotropes")
p3 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Corticotropes")) + NoLegend() + ggtitle("Corticotropes")
p4 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Melanotropes")) + NoLegend() + ggtitle("Melanotropes")
p5 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Gonadotropes")) + NoLegend() + ggtitle("Gonadotropes")
p6 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Thyrotropes")) + NoLegend() + ggtitle("Thyrotropes")
p7 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Pou1f1 Progenitors")) + NoLegend() + ggtitle("Pou1f1 Progenitors")
p8 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Endothelial cells")) + NoLegend() + ggtitle("Endothelial cells")
p9 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Folliculostellate cells")) + NoLegend() + ggtitle("Folliculostellate cells")
p10 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Pericytes")) + NoLegend() + ggtitle("Pericytes")
p11 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Pituicyte")) + NoLegend() + ggtitle("Pituicyte")
p12 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Stem cells")) + NoLegend() + ggtitle("Stem cells")
p13 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "White blood cells")) + NoLegend() + ggtitle("White blood cells")
p14 <- DimPlot(cells, cells.highlight = WhichCells(cells, expression = cell_type == "Red blood cells")) + NoLegend() + ggtitle("Red blood cells")

(p2 | p3 | p4 | p5 | p6) / (p7 | p8 | p9 | p10 | p11) / (p12 | p13 | p14 | patchwork::plot_spacer() | patchwork::plot_spacer())


cells.list <- SplitObject(cells, split.by = "stim")

cells.list.new <- lapply(X = cells.list, FUN = function(obj){
  return(obj %>% 
           NormalizeData() %>%
           FindVariableFeatures() %>%
           ScaleData() %>%
           RunPCA() %>%
           FindNeighbors(dims = 1:50) %>%
           FindClusters() %>%
           RunUMAP(dims = 1:50))
})

## Manual curation is too cumbersome. Use cellassign instead.
d <- "Ctrl"
f <- "Emcn"
c <- CellSelector(FeaturePlot(cells.list.new[[d]], features = f))
HoverLocator(FeaturePlot(cells.list.new[[d]], features = f), information = FetchData(cells.list.new[[d]], vars = c(f)))

cells.list.new[[d]] <- SetIdent(cells.list.new[[d]], cells = c, value = "Pericytes")
table(Idents(cells.list.new[[d]]))
