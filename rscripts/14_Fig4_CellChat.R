library(Seurat)
library(SeuratDisk)
library(CellChat)
library(tidyverse)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

# Part I: Data input & processing and initialization of CellChat object

## Load data
hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")

### Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) is required as input for CellChat analysis.
data.input = GetAssayData(hpcs.lps, slot = "data", assay = "RNA")
meta <- hpcs.lps@meta.data

## Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type_brief")

## Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

### Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

### use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

### set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

