library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(scater)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data_long <- read.table(file = '../data/processed/Pituitary_counts_all.tsv.gz', sep = "\t", header = TRUE)
data <- pivot_wider(data_long, names_from = "cell", values_from = "count", values_fill = 0)
ens_id <- data %>% pull("gene")
data <- data %>% column_to_rownames(var = "gene") %>% as.data.frame()

## Rename rows
# anno <- read.table(file = 'misc/RESOURCES/ensembl/release_102/gtf/mus_musculus/MM.GRCm38.102.annotation.tab', sep = "\t", col.names = c("ensembl_id", "symbol"))
# renameRows <- function(df, anno){
#   stopifnot(is.data.frame(df))
#   stopifnot("ensembl_id" %in% colnames(anno) | "symbol" %in% colnames(anno))
#   
#   if(any(!rownames(df) %in% anno$ensembl_id)){
#     df <- dplyr::slice(df, -which(!rownames(df) %in% anno$ensembl_id))  # df may contain rows not in the annotation table. Remove these rows.
#   }
#   rn <- rownames(df)
#   rn <- anno[match(rn,anno$ensembl_id),"symbol"]
#   dup <- rn[duplicated(rn)]
#   print(paste0("The following genes are duplicated: ", paste0(dup, collapse = ", "), ". Making unique."))
#   for (d in dup){
#     count = 0
#     for (i in 1:length(rn)){
#       if (rn[i] == d) {
#         count = count+1
#         if (count>1){
#           rn[i] = paste0(rn[i],"-",count)
#         }
#       }
#     }
#   }
#   rownames(df) <- rn
#   return(df)
# }
# data <- renameRows(data, anno)

## Attach Feature metadata
cells <- CreateSeuratObject(data, project = "PIT_INFL")
# cells@assays$RNA@meta.features[["ens_id"]] <- ens_id

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
saveRDS(cells.filtered, file = "../data/processed/cells.filtered.sce.rds")  ## for cellassign

cells.new <- as.Seurat(cells.filtered)
cells.new[["nCount_originalexp"]] <- NULL; cells.new[["nFeature_originalexp"]] <- NULL
cells.new <- RenameAssays(cells.new, originalexp = "RNA")
cells.new[["ident"]] <- NULL

## Parse metadata
meta <- read.csv('miscs/metadata.csv')
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

# saveRDS(cells.new, '../data/processed/cells.rds')
SaveH5Seurat(object = cells.new, filename = "../data/processed/cells_pre.h5Seurat", overwrite = T)
