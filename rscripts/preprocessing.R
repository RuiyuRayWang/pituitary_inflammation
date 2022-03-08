library(tidyverse)
library(Seurat)
library(scater)
library(scWidgets)

data_long <- read.table(file = 'data/Pituitary_counts_all.tsv.gz', sep = "\t", header = TRUE)
data <- pivot_wider(data_long, names_from = "cell", values_from = "count", values_fill = 0)
ens_id <- data %>% pull("gene")
data <- data %>% column_to_rownames(var = "gene") %>% as.data.frame()

## Rename rows
anno <- read.table(file = 'misc/MM.GRCm38.102.annotation.tab', sep = "\t", col.names = c("ensembl_id", "symbol"))
data <- renameRows(data, anno)

## Attach Feature metadata
cells <- CreateSeuratObject(data, project = "PIT_INFL")
cells@assays$RNA@meta.features[["ens_id"]] <- ens_id

cells[["barcode"]] <- str_split(string = rownames(cells[[]]), pattern = "_", simplify = TRUE)[,1]
cells[["library"]] <- str_split(string = rownames(cells[[]]), pattern = "_", simplify = TRUE)[,2]

cells <- DietSeurat(cells)
cells.sce <- as.SingleCellExperiment(cells)
rowData(cells.sce)$ens_id <- ens_id  ## Preserve Ens IDs in feature metadata
cells.sce <- logNormCounts(cells.sce)
dim(cells.sce)

## Cell level QC
per.cell <- perCellQCMetrics(cells.sce,
                             subsets = list(Mito = grep("^mt-", rownames(cells.sce)),
                                            Ribo = grep("\\bRp[sl]\\d+[^k]*\\b", rownames(cells.sce))))
colData(cells.sce) <- cbind(colData(cells.sce), per.cell)

qc.stats <- quickPerCellQC(per.cell, percent_subsets = c("subsets_Mito_percent","subsets_Ribo_percent"))
colData(cells.sce) <- cbind(colData(cells.sce), qc.stats)
colSums(as.matrix(qc.stats))
## low_lib_size = 3; low_n_features = 9; high_subsets_Mito_percent = 14; high_subsets_Ribo_percent = 4; discard = 22
cells.filtered <- cells.sce[,!qc.stats$discard]

## Feature level QC
keep_feature <- nexprs(cells.filtered, byrow=TRUE) > 0
## Remove mitochondrial and ribosomal genes which does not contribute to analysis (i.e. cell clustering)
keep_feature[grep("\\bRp[sl]\\d+[^k]*\\b", rownames(cells.sce))] <- FALSE
keep_feature[grep("mt-", rownames(cells.sce))] <- FALSE
cells.filtered <- cells.filtered[keep_feature,]

cells.new <- as.Seurat(cells.filtered)
cells.new[["nCount_originalexp"]] <- NULL; cells.new[["nFeature_originalexp"]] <- NULL
cells.new <- RenameAssays(cells.new, originalexp = "RNA")

## Parse metadata
meta <- read.csv('misc/metadata.csv')
for (ident in names(meta)[2:ncol(meta)]){
  print(paste0("Parsing ",ident,"..."))
  Idents(cells.new) <- "library"
  
  for (lib in meta[["library"]]){
    if(lib %in% Idents(cells.new)){
      cells.new <- SetIdent(cells.new, 
                            cells = WhichCells(cells.new, idents = lib),
                            value = meta %>% filter(library==lib) %>% pull(ident))
    }
  }
  cells.new[[ident]] <- Idents(cells.new)
}
cells.new[["ident"]] <- NULL

cells.new[["treat_dose_duration"]] <- paste(cells.new$treat, paste(cells.new$dose, cells.new$duration, sep = " "), sep = " ")
Idents(cells.new) <- "treat_dose_duration"
cells.new <- RenameIdents(cells.new, 
                          'Saline 0g 6h' = 'Ctrl',
                          'Saline 0g 3h' = 'Ctrl',
                          'LPS 500ug 3w' = 'LPS >3w',
                          'LPS 500ug 4w' = 'LPS >3w',
                          'LPS 500ug 5w' = 'LPS >3w',
                          'LPS 1mg 3w' = 'LPS >3w',
                          'Poly(i:c) 20mg 3w' = 'Poly(i:c) >3w',
                          'Poly(i:c) 20mg 5w' = 'Poly(i:c) >3w',
                          'Poly(i:c) 20mg 8w' = 'Poly(i:c) >3w')
cells.new[["stim"]] <- Idents(cells.new)

saveRDS(cells.new, 'data/cells.rds')
