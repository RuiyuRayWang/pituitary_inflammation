library(Seurat)

cells <- readRDS('data/cells_pre.rds')

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
