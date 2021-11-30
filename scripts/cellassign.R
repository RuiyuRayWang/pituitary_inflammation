library(Seurat)
library(SingleCellExperiment)


Sys.setenv(RETICULATE_PYTHON="/home/luolab/miniconda3/envs/r-reticulate/bin/python")
Sys.setenv(LD_LIBRARY_PATH=
             paste0("/usr/local/cuda/extras/CUPTI:/usr/local/lib:/usr/local/cuda/lib64:/home/luolab/u-net/lib:/home/unetuser/u-net/extlib:",
                    Sys.getenv("LD_LIBRARY_PATH")))
library(tensorflow)
library(reticulate)
# Sys.getenv()
library(cellassign)
tensorflow::tf_config()


library(SingleCellExperiment)
data(example_sce)
data(example_marker_mat)
s <- sizeFactors(example_sce)
fit <- cellassign(exprs_obj = example_sce[rownames(example_marker_mat),],
                  marker_gene_info = example_marker_mat,
                  s = s,
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)


cells <- readRDS("data/cells_pre.rds")
cells.sce <- as.SingleCellExperiment(cells)

# Parse marker list
pituitary_marker_list <- list(
  Somatotropes = c("Gh","Ghrhr","Chgb","Scg2"),
  Lactotropes = c("Prl","Angpt1","Chgb","Scg2"),
  Corticotropes = c("Pomc","Crhr1","Tbx19","Chgb","Scg2"),
  Melanotropes = c("Pomc","Tbx19","Pax7","Pcsk2","Rbfox3","Chgb","Scg2"),
  Gonadotropes = c("Fshb","Lhb","Gnrhr","Cga","Chgb","Scg2"),
  Thyrotropes = c("Tshb","Trhr","Cga","Chgb","Scg2"),
  `Pou1f1 Progenitors` = c("Pbk","Top2a","Mki67"),
  `Red blood cells` = c("Hbb-bt","Hbb-bs"),
  `White blood cells` = c("C1qa","Ctss"), 
  `Folliculostellate cells` = c("S100b","Fxyd1"),
  `Endothelial cells` = c("Pecam1","Emcn","Plvap"),
  `Pituicyte` = c("Gja1","Scn7a","Col25a1"),
  `Pericytes` = c("Col1a1","Dcn","Ogn","Lum","Pdgfrb"),
  `Stem cells` = c("Sox2","Aldh3a1","Aldh1a2","Cyp2f2")
)
mks <- marker_list_to_mat(pituitary_marker_list, include_other = FALSE)
# pheatmap::pheatmap(mks)

# Making sure that the markers exist in the SingleCellExperiment object.
rowData(cells.sce)$Symbol <- rownames(rowData(cells.sce))
marker_in_sce <- match(rownames(mks), rowData(cells.sce)$Symbol)
stopifnot(all(!is.na(marker_in_sce)))

# Subset `sce` object to retain only marker genes.
sce_marker <- cells.sce[marker_in_sce, ]
stopifnot(all.equal(rownames(mks), rowData(sce_marker)$Symbol))

# Call `cellassign` passing in the `SingleCellExperiment`, marker info, the size factors we’ve calculated, as well as 
# various other options
counts(sce_marker) <- as.matrix(counts(sce_marker))
print(dim(sce_marker))
print(dim(mks))

fit <- cellassign(
  exprs_obj = sce_marker,
  marker_gene_info = mks,
  s = sizeFactors(cells.sce),
  shrinkage = TRUE,
  max_iter_adam = 50,
  min_delta = 2,
  verbose = TRUE
)
