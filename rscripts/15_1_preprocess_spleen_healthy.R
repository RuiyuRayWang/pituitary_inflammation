library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(scater)
library(GEOquery)
# library(org.Mm.eg.db)
# library(awesomeSC)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

if(!file.exists('../data/GSE132642/GSE132642_RAW.tar')){
  options(timeout = 6000)
  getGEOSuppFiles(GEO = "GSE132642", makeDirectory = T, baseDir = "../data/")
}

### TODO: untar GSE, get meta

meta <- read.csv(file.path(wd,"Hilton_spleen","pbio.3000528.s031.csv"))
meta$cell <- sub(pattern = "-", replacement = ".", meta$cell)
meta$cell <- sub(pattern = "^([^\\.]*\\.[^\\.]*\\.[^\\.]*\\.[^\\.]*)\\.", replacement = "\\1_", meta$cell)
tab <- read.table('/mnt/ZA1BT1ER/yanting/pituitary_scRNAseq_analysis/QC/gencode.vM21.annotation.tab', col.names = c("ens_id","symbol"))
tab$ens_id <- sub(pattern = "\\.\\d+", "", x = tab$ens_id)

mt_genes_ens_id <- tab[grep(pattern = "mt-", tab$symbol),"ens_id"]
rb_genes_ens_id <- tab[c(grep(pattern = "^Rp[sl]\\d+[^k]",tab$symbol),
                         grep(pattern = "^Rp[sl]p",tab$symbol)),"ens_id"]


# Hilton et al. Mouse. LPS
## Parse data
dat_1 <- read.csv(file = file.path(wd,"Hilton_spleen","GSM3885345_mouse.spleen.LPS_1_1.csv.gz"), row.names = 1)
dat_2 <- read.csv(file = file.path(wd,"Hilton_spleen","GSM3885346_mouse.spleen.LPS_1_2.csv.gz"), row.names = 1)
dat_3 <- read.csv(file = file.path(wd,"Hilton_spleen","GSM3885347_mouse.spleen.LPS_2_1.csv.gz"), row.names = 1)
dat_4 <- read.csv(file = file.path(wd,"Hilton_spleen","GSM3885348_mouse.spleen.LPS_2_2.csv.gz"), row.names = 1)

dat_1 <- dat_1 %>% rownames_to_column(var = "ensembl_id")
dat_2 <- dat_2 %>% rownames_to_column(var = "ensembl_id")
dat_3 <- dat_3 %>% rownames_to_column(var = "ensembl_id")
dat_4 <- dat_4 %>% rownames_to_column(var = "ensembl_id")
dat <- dat_1 %>%
  left_join(dat_2, by = "ensembl_id") %>%
  left_join(dat_3, by = "ensembl_id") %>%
  left_join(dat_4, by = "ensembl_id") %>%
  column_to_rownames(var = "ensembl_id")
rm(list = ls(pattern = "dat_.*\\d"))

## Create SCE object
cells.sce <- SingleCellExperiment(list(counts=as.matrix(dat)), colData = meta)  # Appending metadata
all(rownames(colData(cells.sce)) == colData(cells.sce)$cell)  # Check that cell names match

## QC
## Remove features with NA matches in gencode gtf annotation
keep_features = !is.na(tab[match(rownames(cells.sce), tab$ens_id),"symbol"])
cells.sce <- cells.sce[keep_features,]

## Rename rows to gene symbols
rowData(cells.sce) <- tab[match(rownames(cells.sce), tab$ens_id),"symbol"]
rownames(cells.sce) <- rowData(cells.sce)$value
# per.cell <- perCellQCMetrics(cells.sce,
#                              subsets = list(Mito = match(mt_genes_ens_id, rownames(cells.sce))))
# cells.sce <- addPerCellQC(cells.sce,
#                                         subsets = list(
#                                           Mito=match(mt_genes_ens_id[mt_genes_ens_id %in% rownames(cells.sce)], 
#                                                      rownames(cells.sce)),
#                                           Ribo=match(rb_genes_ens_id[rb_genes_ens_id %in% rownames(cells.sce)], 
#                                                      rownames(cells.sce))))
cells.sce <- addPerCellQC(cells.sce,
                          subsets = list(
                            Mito=grep(pattern = "mt-", rownames(cells.sce)),
                            Ribo=c(grep(pattern = "^Rp[sl]\\d+[^k]",rownames(cells.sce)),
                                   grep(pattern = "^Rp[sl]p",rownames(cells.sce)))))
qc.stats <- quickPerCellQC(colData(cells.sce), percent_subsets=c("subsets_Mito_percent","subsets_Ribo_percent"))
colSums(as.matrix(qc.stats))
cells.sce <- cells.sce[,!qc.stats$discard]
# cells.sce <- logNormCounts(cells.sce)

cells.seurat <- as.Seurat(cells.sce, data = NULL)

# Standard Seurat Workflow
cells.seurat <- NormalizeData(cells.seurat)
cells.seurat <- FindVariableFeatures(cells.seurat)
cells.seurat <- ScaleData(cells.seurat)
cells.seurat <- RunPCA(cells.seurat, npcs = 50)
cells.seurat <- JackStraw(cells.seurat, dims = 50)
cells.seurat <- ScoreJackStraw(cells.seurat, dims = 1:50)
JackStrawPlot(cells.seurat, dims = 1:50)
plotly::ggplotly()
dim_use <- 1:50
cells.seurat <- FindNeighbors(cells.seurat, dims = dim_use)
cells.seurat <- FindClusters(cells.seurat)
cells.seurat <- RunUMAP(cells.seurat, dims = dim_use)
cells.seurat <- RunTSNE(cells.seurat, dims = dim_use)

# DimPlot(cells.seurat, group.by = "converged.cell.type")
Idents(cells.seurat) <- "converged.cell.type"
cells.seurat <- RenameIdents(cells.seurat, `T-cell` = "T cells", `Plasmacytoid-dendritic-cell` = "pDCs", `NK-cell` = "NK cells",
                             `B-cell` = "B cells", `Macrophage` = "Macrophages", `Monocyte` = "Monocytes", `Neutrophil` = "Neutrophils", 
                             `Plasma-cell` = "B cells")
cells.seurat[["cell_type"]] <- Idents(cells.seurat)

saveRDS(cells.seurat, file.path(wd,"Hilton_spleen","hilton.mouse.lps.seurat.rds"))
