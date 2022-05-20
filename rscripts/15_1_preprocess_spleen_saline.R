library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(scater)
library(GEOquery)
library(org.Mm.eg.db)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("renameRows.R")

if(!file.exists('../data/GSE132642/GSE132642_RAW.tar')){
  options(timeout = 6000)
  getGEOSuppFiles(GEO = "GSE132642", makeDirectory = T, baseDir = "../data/")
}

untar(tarfile = '../data/GSE132642/GSE132642_RAW.tar', 
      files = c("GSM3885341_mouse.spleen.saline_1_1.csv.gz",
                "GSM3885342_mouse.spleen.saline_1_2.csv.gz",
                "GSM3885343_mouse.spleen.saline_2_1.csv.gz",
                "GSM3885344_mouse.spleen.saline_2_2.csv.gz"),
      exdir = "../data/GSE132642/", 
)

## Load data
data <- rownames_to_column(read.csv(file = "../data/GSE132642/GSM3885341_mouse.spleen.saline_1_1.csv.gz", row.names = 1), "ens_id") %>%
  left_join(rownames_to_column(read.csv(file = "../data/GSE132642/GSM3885342_mouse.spleen.saline_1_2.csv.gz", row.names = 1), "ens_id")) %>%
  left_join(rownames_to_column(read.csv(file = "../data/GSE132642/GSM3885343_mouse.spleen.saline_2_1.csv.gz", row.names = 1), "ens_id")) %>%
  left_join(rownames_to_column(read.csv(file = "../data/GSE132642/GSM3885344_mouse.spleen.saline_2_2.csv.gz", row.names = 1), "ens_id")) %>%
  column_to_rownames("ens_id")
gc()

## Rename rows
anno <- data.frame(
  ensembl_id = rownames(data), 
  symbol = mapIds(
    x = org.Mm.eg.db,
    keys = rownames(data), 
    column = "SYMBOL", 
    keytype = "ENSEMBL",
    multiVals = "first")
)
anno <- anno[!is.na(anno$symbol),]
data <- renameRows(df = data, anno = anno)

## Parse metadata
meta <- read.csv("../data/GSE132642/pbio.3000528.s030.csv")
meta$cell <- sub(pattern = "-", replacement = ".", meta$cell)
meta$cell <- sub(pattern = "^([^\\.]*\\.[^\\.]*\\.[^\\.]*\\.[^\\.]*)\\.", replacement = "\\1_", meta$cell)

## Create SCE object
cells.sce <- SingleCellExperiment(list(counts=as.matrix(data)), colData = meta)  # Appending metadata

## QC
## Rename rows to gene symbols
# rowData(cells.sce) <- tab[match(rownames(cells.sce), tab$ens_id),"symbol"]
# rownames(cells.sce) <- rowData(cells.sce)$value 

## Cell level QC
per.cell <- perCellQCMetrics(cells.sce,
                             subsets = list(Mito = grep("^mt-", rownames(cells.sce)),
                                            Ribo = grep("\\bRp[sl]\\d+[^k]*\\b", rownames(cells.sce))))
qc.stats <- quickPerCellQC(per.cell, percent_subsets=c("subsets_Mito_percent","subsets_Ribo_percent"))
colSums(as.matrix(qc.stats))
colData(cells.sce) <- cbind(colData(cells.sce), cbind(per.cell, qc.stats))
cells.filtered <- cells.sce[,!cells.sce$discard]

## Feature level QC
keep_feature <- nexprs(cells.filtered, byrow=TRUE) > 0
## Remove mitochondrial and ribosomal genes which does not contribute to analysis (i.e. cell clustering)
keep_feature[grep("\\bRp[sl]\\d+[^k]*\\b", rownames(cells.sce))] <- FALSE
keep_feature[grep("mt-", rownames(cells.sce))] <- FALSE
cells.filtered <- cells.filtered[keep_feature,]
dim(cells.filtered)

cells <- as.Seurat(cells.sce, data = NULL)
cells[["nCount_originalexp"]] <- NULL; cells[["nFeature_originalexp"]] <- NULL
cells <- RenameAssays(cells, originalexp = "RNA")

# Standard Seurat Workflow
dim_use = 1:50
cells <- NormalizeData(cells) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = dim_use) %>%
  FindClusters() %>%
  RunUMAP(dims = dim_use)

# DimPlot(cells, group.by = "converged.cell.type")
Idents(cells) <- "converged.cell.type"
cells <- RenameIdents(cells, `T-cell.1` = "T cells", `T-cell.2` = "T cells", `Plasmacytoid-dendritic-cell` = "pDCs", `NK-cell` = "NK cells",
                      `B-cell` = "B cells", `Macrophage` = "Macrophages", `Monocyte` = "Monocytes", `Neutrophil` = "Neutrophils")
cells[["cell_type"]] <- Idents(cells)
cells[["treat"]] <- "Saline"

SaveH5Seurat(object = cells, filename = "../data/GSE132642/hilton_mouse_spleen_saline.h5Seurat", overwrite = T)
