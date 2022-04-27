# ## Installation instructions
# 
# ## Build multiple R versions: https://ruiyuraywang.github.io/post/mult-r/ and switch to R-3.6.2
# 
# ## Install some package dependencies
# install.packages("reticulate")
# install.packages("devtools")
# install.packages("BiocManager")
# install.packages("Seurat")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("splatter")
# devtools::install_github("erichson/rsvd@v1.0.0")
# BiocManager::install("BiocSingular")
# devtools::install_version("RcppAnnoy", "0.0.16", repos="http://cran.us.r-project.org")
# BiocManager::install("BiocNeighbors", version="3.10")
# BiocManager::install("scater")
# BiocManager::install("scran")
# 
# ## Create dedicated environment in shell
# # conda create -n cellassign python=3.7
# 
# ## Setup environment to use
# Sys.setenv(RETICULATE_PYTHON="/home/luolab/miniconda3/envs/cellassign/bin/python")
# library(reticulate)
# reticulate::use_condaenv("cellassign")
# 
# ## Build core packages for cellassign
# 
# ## Install tensorflow
# ## NOT DO! Use devtools and build tensorflow from shell!
# # install.packages("tensorflow")
# # tensorflow::install_tensorflow(extra_packages="tensorflow-probability", version = "2.1.0")
# devtools::install_github("rstudio/tensorflow@v2.4.0")
# 
# ## In shell, activate conda environment `cellassign` and do
# # $ pip install tensorflow==2.4.0
# # $ pip install tensorflow-probability==0.12.0
# # $ pip install tensorflow-gpu==2.4.0
# 
# ## Install cellassign
# devtools::install_github("Irrationone/cellassign")

## Run pipe
library(SingleCellExperiment)

library(tensorflow)
tensorflow::tf_config()  # Make sure tensorflow is installed properly

library(cellassign)

cells.sce <- readRDS("data/cells.filtered.sce.rds")
s <- cells.sce$sizeFactor

# Parse marker list
pituitary_marker_list <- list(
  Somatotropes = c("Gh","Ghrhr","Chgb","Scg2"),
  Lactotropes = c("Prl","Angpt1","Chgb","Scg2"),
  Corticotropes = c("Pomc","Crhr1","Tbx19","Chgb","Scg2"),
  Melanotropes = c("Pomc","Tbx19","Pax7","Pcsk2","Rbfox3","Chgb","Scg2"),
  Gonadotropes = c("Fshb","Lhb","Gnrhr","Cga","Chgb","Scg2"),
  Thyrotropes = c("Tshb","Trhr","Cga","Chgb","Scg2"),
  `Pou1f1 Progenitors` = c("Pbk","Top2a","Mki67"),
  `Red Blood Cells` = c("Hbb-bt","Hbb-bs"),
  `White Blood Cells` = c("C1qa","Ctss"), 
  `Folliculostellate Cells` = c("S100b","Fxyd1"),
  `Endothelial Cells` = c("Pecam1","Emcn","Plvap"),
  `Pituicyte` = c("Gja1","Scn7a","Col25a1"),
  `Pericytes` = c("Col1a1","Dcn","Ogn","Lum","Pdgfrb"),
  `Stem Cells` = c("Sox2","Aldh3a1","Aldh1a2","Cyp2f2")
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
  s = s,
  shrinkage = TRUE,
  max_iter_adam = 50,
  min_delta = 2,
  verbose = TRUE
)
saveRDS(fit, file = paste0("data/fit_cellassign_",Sys.Date(),".rds"))

## Switch R version back
## Assign to objects
library(SeuratDisk)
# cells.seurat <- readRDS('data/cells_pre.rds')
cells.seurat <- LoadH5Seurat("data/cells_pre.h5Seurat")
fit <- readRDS(file = "data/fit_cellassign_2022-04-27.rds")
cells.seurat$cell_type_cellassign <- fit$cell_type
# saveRDS(cells.seurat, file = "data/cells.rds")
SaveH5Seurat(cells.seurat, filename = "data/cells_assigned.h5Seurat", overwrite = T)  ## use h5 format whenever possible

